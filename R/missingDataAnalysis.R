#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param long_data PARAM_DESCRIPTION
#' @param plot.dir PARAM_DESCRIPTION, Default: 'Missing_Data'
#' @param status.config PARAM_DESCRIPTION, Default: NULL
#' @param growth.config PARAM_DESCRIPTION, Default: NULL
#' @param focus.variable PARAM_DESCRIPTION, Default: 'SCALE_SCORE'
#' @param student.factors PARAM_DESCRIPTION, Default: NULL
#' @param factor.labels PARAM_DESCRIPTION, Default: NULL
#' @param group PARAM_DESCRIPTION, Default: NULL
#' @param output.plots PARAM_DESCRIPTION, Default: FALSE
#' @param output.type PARAM_DESCRIPTION, Default: 'SVG'
#' @param ... PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[data.table]{J}}, \code{\link[data.table]{setattr}}, \code{\link[data.table]{setkey}}, \code{\link[data.table]{dcast.data.table}}, \code{\link[data.table]{data.table-package}}
#'  \code{\link[stats]{setNames}}
#'  \code{\link[naniar]{mcar_test}}, \code{\link[naniar]{gg_miss_upset}}, \code{\link[naniar]{gg_miss_fct}}
#'  \code{\link[mice]{md.pattern}}
#'  \code{\link[finalfit]{missing_plot}}, \code{\link[finalfit]{missing_pairs}}, \code{\link[finalfit]{extract_variable_label}}, \code{\link[finalfit]{summary_factorlist}}, \code{\link[finalfit]{ff_merge}}
#'  \code{\link[ggplot2]{scale_colour_gradient}}, \code{\link[ggplot2]{scale_manual}}, \code{\link[ggplot2]{labs}}
#'  \code{\link[tools]{toTitleCase}}
#'  \code{\link[fixest]{fixef_reexported}}
#'  \code{\link[utils]{combn}}
#'  \code{\link[grDevices]{pdf}}, \code{\link[grDevices]{dev2}}, \code{\link[grDevices]{recordPlot}}
#'  \code{\link[VIM]{mosaicMiss}}, \code{\link[VIM]{spineMiss}}
#'  \code{\link[svglite]{svglite}}
#' @rdname missingDataAnalysis
#' @export 
#' @author AUTHOR [AUTHOR_2]
#' @keywords KEYWORD_TERM
#' @author AUTHOR [AUTHOR_2]
#' @keywords KEYWORD_TERM
#' @importFrom data.table SJ setnames setkey setkeyv dcast setattr data.table
#' @importFrom stats setNames
#' @importFrom naniar mcar_test gg_miss_upset gg_miss_fct
#' @importFrom mice md.pattern
#' @importFrom finalfit missing_plot missing_pairs extract_variable_label summary_factorlist ff_merge
#' @importFrom ggplot2 scale_fill_gradient scale_fill_manual ylab xlab
#' @importFrom tools toTitleCase
#' @importFrom fixest fixef
#' @importFrom utils combn
#' @importFrom grDevices pdf dev.control recordPlot
#' @importFrom VIM mosaicMiss spineMiss
#' @importFrom svglite svglite
missingDataAnalysis =
  function(
    long_data,
    plot.dir = "Missing_Data",
    status.config = NULL,
    growth.config = NULL,
    focus.variable = "SCALE_SCORE",
    student.factors = NULL,
    factor.labels = NULL,
    group = NULL,
    output.plots = FALSE,
    output.type = "SVG", # "PDF"
    ...
  ) {

  ###  Avoid "Undefined global functions or variables:" from R CMD check
# GRADE <- VALID_CASE <- Z_SCORE <- TMP_MEAN_Z <-
  ID <- YEAR <- CONTENT_AREA <- TMP_STFAC <- Missing <- Avg_Focus_Var_Z <- NULL


  ###  Initial checks
  if (output.plots) {
    diagn.dir <-
      file.path(plot.dir, "missingDataAnalysis", "missing_data_plots")
    if (!dir.exists(diagn.dir)) dir.create(diagn.dir, recursive = TRUE)
  }

  if (length(group) > 1) {
    stop("Only one 'group' variable can be used as an impution factor.")
  }

  ###  Combine and augment config lists
  configs <- status.config

  if (!is.null(status.config)) {
    for (f in seq(configs)) configs[[f]][["analysis.type"]] <- "STATUS"
  }
  if (!is.null(growth.config)) {
    for (f in seq(growth.config)) growth.config[[f]][["analysis.type"]] <- "GROWTH"
    configs <- c(configs, growth.config)
  }

  if (is.null(configs)) {
    stop("Either a 'growth.config' or 'status.config' must be supplied")
  }

  long.to.wide.vars <- c(focus.variable, group, student.factors) # "CONTENT_AREA", "GRADE",

  ###  Cycle through configs to get results by cohort
  res.list <- vector(mode = "list", length = length(configs))
  names(res.list) <- names(configs)

  for (K in seq(configs)) {
    cohort.iter <- res.list[[K]][["config"]] <- configs[[K]]
    names(cohort.iter) <- gsub("^sgp[.]", "", names(cohort.iter))

    prior.years <- utils::head(unique(cohort.iter[["panel.years"]]), -1)
    current.year <- utils::tail(unique(cohort.iter[["panel.years"]]), 1)
    current.grade <- utils::tail(unique(cohort.iter[["grade.sequences"]]), 1)
    if (length(subject <- unique(cohort.iter[["content.areas"]])) > 1) subject <- NULL
    grade.length <- length(cohort.iter[["grade.sequences"]])

    cohort.lookup <-
      data.table::SJ(
        "VALID_CASE",
        utils::tail(cohort.iter[["content.areas"]], grade.length),
        utils::tail(cohort.iter[["panel.years"]], grade.length),
        cohort.iter[["grade.sequences"]]
      ) |>
      data.table::setnames(getKey(long_data) %w/o% "ID") |>
        data.table::setkey(YEAR) |> # ensure lookup table is ordered by years.
          data.table::setkey(NULL)  # NULL out key so doesn't corrupt the join in dcast.

    data.table::setkeyv(long_data, getKey(long_data))

    tmp_long <-
      long_data[cohort.lookup][, c(getKey(long_data), long.to.wide.vars), with = FALSE]

    if (cohort.iter[["analysis.type"]] == "GROWTH") {
      ###  convert long to wide
      wide_data <-
        data.table::dcast(data = tmp_long,
                          formula = ID ~ YEAR + CONTENT_AREA,
                          sep = ".",
                          drop = FALSE,
                          value.var = c("VALID_CASE", long.to.wide.vars))

      all.years <- c(prior.years, current.year)
      tmp.vcase <- paste0(paste0(focus.variable, ".", all.years), ".",
                          rep(unique(cohort.iter[["content.areas"]]), each = length(all.years)))
      excl.idx <- Reduce(intersect,
                         lapply(tmp.vcase, \(f) which(is.na(wide_data[, get(f)]))))
      if (length(excl.idx)) {
        wide_data <- wide_data[-excl.idx]
      }

      subset.vars <- focus.vars <-
        grep(paste0(focus.variable, "[.]", collapse = "|"), names(wide_data), value = TRUE)
      if (length(group)) {
        group.vars <- grep(group, names(wide_data), value = TRUE)
        subset.vars <- c(subset.vars, group.vars)
      } else {
        group.vars <- NULL
      }
      if (length(student.factors)) {
        for (sf in student.factors) {
          sf.wide.nms <- rev(grep(sf, names(wide_data), value = TRUE))
          wide_data[, eval(sf) := get(sf.wide.nms[1])]
          for (r in seq(sf.wide.nms)[-1]) {
            wide_data[is.na(get(sf)), eval(sf) := get(sf.wide.nms[r])]
          }
        }
        subset.vars <- c(subset.vars, paste0(student.factors, "$"))
      }

      #  create subset of wide data with only variables to be used in analysis
      wide_data <-
        wide_data[,
          grep(
            paste0("^ID$|VALID_CASE|", paste(subset.vars, collapse = "|")),
            names(wide_data)
          ),
          with = FALSE
        ]

      #  remove columns that are all NA (e.g., SGP for 3rd grade priors)
      na.vars <-
        unlist(lapply(names(wide_data),
                      \(f) all(is.na(wide_data[, get(f)]))))
      if (any(na.vars)) {
        na.var.nms <- names(wide_data)[na.vars]
        focus.vars <- focus.vars %w/o% na.var.nms
        student.factors <- student.factors %w/o% na.var.nms
        group.vars <- group.vars %w/o% na.var.nms
        wide_data <- wide_data[, names(wide_data)[!na.vars], with = FALSE]
      }

      #  assign labels to variables for use in the finalfit package functions
      for (fv in focus.vars) {
        fv.lab <- paste(rev(strsplit(fv, "[.]")[[1]][-1]), collapse = " ")
        # finalfit::ff_label doesn't work...
        data.table::setattr(wide_data[[fv]], "label", fv.lab)
      }

      if (length(factor.labels)) {
        for (sfac in names(factor.labels)) {
          data.table::setattr(wide_data[[sfac]], "label", factor.labels[[sfac]])
        }
      } # finalfit::extract_variable_label(wide_data)

      group.mean.vars <- NULL # Add these???

    } else {  #  END "GROWTH"  --  Begin "STATUS"
      ##  Create wide_data and tmp_long_priors
      status.lookup <- cohort.lookup[YEAR == current.year]
      data.table::setkeyv(tmp_long, getKey(tmp_long))

      wide.fmla <- stats::as.formula(paste("ID",
                      ifelse(is.null(group), "", paste("+", group)),
                      "~ YEAR + CONTENT_AREA"))
      wide_data <-
        data.table::dcast(
          data = tmp_long[status.lookup],
          formula = wide.fmla,
          sep = ".", drop = FALSE,
          value.var = c("VALID_CASE", focus.variable, student.factors)
        )

      if (!is.null(group)) {
        group.lookup <- unique(tmp_long[status.lookup, c("ID", group), with = FALSE])
        wide_data <- wide_data[group.lookup]
      }

      focus.vars <-
        grep(paste0(focus.variable, "[.]", collapse = "|"), names(wide_data), value = TRUE)
      if (length(student.factors)) {
        for (sf in student.factors) {
          sf.wide.nms <- rev(grep(sf, names(wide_data), value = TRUE))
          wide_data[, eval(sf) := get(sf.wide.nms[1])]
          for (r in seq(sf.wide.nms)[-1]) {
            wide_data[is.na(get(sf)), eval(sf) := get(sf.wide.nms[r])]
          }
        }
      }
      wide_data <-
        wide_data[,
          c("ID", group, student.factors, focus.vars,
            grep("VALID_CASE", names(wide_data), value = TRUE)
          ),
          with = FALSE
        ]
      #  remove columns that are all NA (e.g., SGP for 3rd grade priors)
      na.vars <-
        unlist(lapply(names(wide_data),
                      \(f) all(is.na(wide_data[, get(f)]))))
      if (any(na.vars)) {
        na.var.nms <- names(wide_data)[na.vars]
        focus.vars <- focus.vars %w/o% na.var.nms
        student.factors <- student.factors %w/o% na.var.nms
        group.vars <- group.vars %w/o% na.var.nms
        wide_data <- wide_data[, names(wide_data)[!na.vars], with = FALSE]
      }

      if (length(group)) {
        group.vars <- group
        priors.lookup <- cohort.lookup[YEAR != current.year]
        tmp_long_priors <-
          tmp_long[priors.lookup][,
            c("YEAR", "CONTENT_AREA", "GRADE", focus.variable, group),
            with = FALSE
          ]

        ##  Prior year(s) summaries to use
        smry.eval.expression <-
          paste0("MEAN_DIFF__", focus.variable,
                 " = ", "mean(", focus.variable, ", na.rm=TRUE)")
        smry.eval.expression <-
          stats::setNames(smry.eval.expression,
                          sub("^(.*) = .*", "\\1", smry.eval.expression))

        tmp_grp_smry <- tmp_long_priors[!is.na(get(group)),
                          lapply(smry.eval.expression, \(f) eval(parse(text = f))),
                            keyby = c("YEAR", "CONTENT_AREA", group)]

        tmp_grp_smry <- data.table::dcast(tmp_grp_smry, get(group) ~ YEAR + CONTENT_AREA,
                                          sep = ".", value.var = names(smry.eval.expression))

        data.table::setnames(tmp_grp_smry,
                             c(group,
                               paste0(names(smry.eval.expression), ".", names(tmp_grp_smry)[-1])))

        wide_data <- merge(wide_data, tmp_grp_smry, by = group, all.x = TRUE)

        ##  Put in cross school mean for schools with no students in prior years
        group.mean.vars <- grep("MEAN_DIFF__", names(wide_data), value = TRUE)
        for (gmv in group.mean.vars) {
          tmp.grp.mean <- mean(wide_data[, get(gmv)], na.rm = TRUE)
          wide_data[is.na(get(gmv)), eval(gmv) := tmp.grp.mean]
        }

        for (ca in unique(cohort.iter[["content.areas"]])) {
          tmp.m.vars <- grep(paste0("MEAN_DIFF__.+.", ca), group.mean.vars, value = TRUE)
          fv <- paste(focus.variable, current.year, ca, sep = ".")
          for (M in tmp.m.vars) {
            wide_data[, eval(M) := get(M) * -1]
            wide_data[, eval(M) := rowSums(.SD), .SDcols = c(M, fv)]
            wide_data[is.na(get(M)), eval(M) := 0]
            # wide_data[, as.list(summary(get(M))), keyby = group]
          }
        }
      } else {
        group.vars <- group.mean.vars <- NULL
      }

      #  assign labels to variables for use in the finalfit package functions
      for (fv in focus.vars) {
        fv.lab <- paste(rev(strsplit(fv, "[.]")[[1]][-1]), collapse = " ")
        data.table::setattr(wide_data[[fv]], "label", fv.lab)
      }

      if (length(factor.labels)) {
        for (sfac in names(factor.labels)) {
          data.table::setattr(wide_data[[sfac]], "label", factor.labels[[sfac]])
        }
      } # finalfit::extract_variable_label(wide_data)
    } ###  END "STATUS"

    #####
    ###   Analyze Missing Data/Patterns
    #####

    res.list[[K]][["mcar_test"]] <-
      naniar::mcar_test(data = wide_data[, ..focus.vars])

    res.list[[K]][["patterns"]] <-
      mice::md.pattern(wide_data[, ..focus.vars], plot = FALSE)
    non.miss <- which(res.list[[K]][["patterns"]][, ncol(res.list[[K]][["patterns"]])] == 0L)
    md.freqs <- (as.numeric(rownames(res.list[[K]][["patterns"]])) / nrow(wide_data))[-(non.miss)]
    md.freqs <- round(md.freqs / sum(md.freqs, na.rm = TRUE), 3) # re-weight to sum to 1
    patts.tk <- which(md.freqs > 0.005)
    res.list[[K]][["upset_plot"]] <- try(
      naniar::gg_miss_upset(
        data = wide_data[, ..focus.vars],
        nsets = length(focus.vars), nintersects = length(patts.tk)
      )
    )
    res.list[[K]][["missing_plot"]] <- try(
      finalfit::missing_plot(
        .data = wide_data[, ..focus.vars],
        title = "", use_labels = TRUE,
        plot_opts = # match `VIM` colors
          ggplot2::scale_fill_gradient(low = "#87ceeb", high = "#ff3333") # low = "#80d9f7", high = "#fc3f3f"
      )
    )
    res.list[[K]][["missing_pairs"]] <- suppressMessages(
      finalfit::missing_pairs(
        .data = wide_data[, ..focus.vars],
        title = "", position = "fill", use_labels = TRUE
      ) + ggplot2::scale_fill_manual(values = c("#87ceeb", "#ff3333")) # c("#80d9f7", "#fc3f3f")
    )

    if (length(student.factors)) {
      tmp_data <- wide_data[, c(focus.vars, student.factors), with = FALSE]
      data.table::setnames(tmp_data, finalfit::extract_variable_label(tmp_data))
      tmp.sfacs <-
        if (is.null(factor.labels)) {
          student.factors
        } else {
          finalfit::extract_variable_label(wide_data)[student.factors]
        }

      for (sfac in tmp.sfacs) {
        tmp_data[, TMP_STFAC := factor(get(sfac))]
        data.table::setattr(tmp_data[["TMP_STFAC"]], "levels",
                            gsub(".*: ", "", levels(tmp_data[["TMP_STFAC"]])))
        res.list[[K]][[names(tmp.sfacs)[which(tmp.sfacs == sfac)]]][["factor_plot"]] <-
          naniar::gg_miss_fct(x = tmp_data[, !..tmp.sfacs], fct = TMP_STFAC) +
            ggplot2::ylab(gsub("_|[.]", " ", focus.variable)) +
            ggplot2::xlab(gsub("_|[.]", " ", sfac))
      }
      tmp_data[, TMP_STFAC := NULL]
    }

    for (fv in focus.vars) {
      if (cohort.iter[["analysis.type"]] == "GROWTH") {
        gv <- paste0(group, ".", sub(".+?[.]", "", fv))
      } else {
        gv <- group
      }
      vc <- paste0("VALID_CASE.", sub(".+?[.]", "", fv))
      excluded <- c("ID", group.vars,
                    grep(sub("(.+)[.].*", "\\1", fv), focus.vars, value = TRUE),
                    grep("VALID_CASE", names(wide_data), value = TRUE),
                    grep(sub(paste0(".*", current.year, "[.]"), "", fv),
                         group.mean.vars, invert = TRUE, value = TRUE))
      iv <- names(wide_data) %w/o% excluded
      dv <- "Missing"

      wide_data[, Missing :=
          factor(ifelse(is.na(get(fv)), 1, 0),
                 levels = c(0, 1), labels = c("Observed", "Missing"))]

      tmp_wide <- wide_data[, c(..vc, ..dv, ..iv, ..gv)][get(vc) == "VALID_CASE", ]
      fct.iv <- iv %w/o% student.factors
      stu.iv <- iv %w/o% fct.iv

      ##  Create a single mean `focus.variable` - too many missing in some years.
      ##  Also tried stepwise regressions with each year, but those got pretty messy.
      ##  This provides the (?) best/most complete cases for regressions.
      for (f in fct.iv) {
        tmp_wide[, eval(gsub(focus.variable, "Z_SCORE", f)) := scale(get(f))]
      }
      tmp_wide[, Avg_Focus_Var_Z := rowMeans(.SD, na.rm = TRUE),
                 .SDcols = gsub(focus.variable, "Z_SCORE", fct.iv)]
      tmp_wide[, gsub(focus.variable, "Z_SCORE", fct.iv) := NULL]
      tmp.label <- tools::toTitleCase(tolower(strsplit(focus.variable, "_|[.]")[[1]]))

      data.table::setattr(tmp_wide[["Avg_Focus_Var_Z"]], "label",
        paste(c("Mean", tmp.label, "(Standardized)"), collapse = " "))

      all.iv <- c(fct.iv, "Avg_Focus_Var_Z", stu.iv)
      mod.iv <- c(stu.iv, "csw0(Avg_Focus_Var_Z)")

      if (length(group)) {
        res.list[[K]][[fv]][["FE"]][["model"]] <-
          glmfixed(.data = tmp_wide,
                   dependent = dv,
                   explanatory = mod.iv,
                   fixed_effect = gv,
                   lean = TRUE
          )
        fxd.tab <-
          fit2df.fixest(utils::tail(res.list[[K]][[fv]][["FE"]][["model"]], 1)[[1]],
                        estimate_suffix = " (fixed effects)")
        res.list[[K]][[fv]][["FE"]][["fixed_effects"]] <-
          fixest::fixef(
            glmfixed(
              .data = tmp_wide,
              dependent = dv,
              explanatory = c(stu.iv, "Avg_Focus_Var_Z"),
              fixed_effect = gv
            )
          )
      } else {
        fxd.tab <- NULL
      }
      res.list[[K]][[fv]][["GLM"]][["model"]] <-
        glmfixed(
          .data = tmp_wide,
          dependent = dv,
          explanatory = mod.iv,
          lean = TRUE
        )
      glm.tab <-
        fit2df.fixest(
          utils::tail(res.list[[K]][[fv]][["GLM"]][["model"]], 1)[[1]],
          estimate_suffix = " (multivariable)"
        )

      res.list[[K]][[fv]][["table"]] <-
        finalfit::summary_factorlist(
          .data = tmp_wide,
          dependent = dv,
          explanatory = all.iv,
          na_include = TRUE,
          p = TRUE, # p_cat = "fisher", p_cont_para = "t.test",
          fit_id = TRUE
        ) |>
          # https://ivelasq.rbind.io/blog/understanding-the-r-pipe/
          # https://stackoverflow.com/questions/30604107/r-conditional-evaluation-when-using-the-pipe-operator
          (\(.) {
            if (is.null(fxd.tab)) return(.)
            finalfit::ff_merge(factorlist = .,
                               fit2df_df = fxd.tab)
          })() |>
            finalfit::ff_merge(fit2df_df = glm.tab,
                              last_merge = TRUE) |>
              data.table::data.table()

      if (length(student.factors)) {
        if (length(student.factors > 1)) {
          perms <- utils::combn(student.factors, 2, simplify = FALSE)
          for (cmb in seq(perms)) {
            grDevices::pdf(NULL)
            grDevices::dev.control(displaylist = "enable")
            VIM::mosaicMiss(
              x = wide_data[, mget(c(fv, perms[[cmb]]))], highlight = fv,
              plotvars = perms[[cmb]], miss.labels = FALSE, only.miss = FALSE,
              col = c("#87ceeb", "#ff3333"), # skyblue and slightly muted 'red'
              labels = list(
                set_varnames =
                  finalfit::extract_variable_label(wide_data)[perms[[cmb]]]
              )
            )
            res.list[[K]][[fv]][[paste(perms[[cmb]], collapse = "__x__")]] <-
              grDevices::recordPlot()
            invisible(grDevices::dev.off())
          }
        } else {
          grDevices::pdf(NULL)
          grDevices::dev.control(displaylist = "enable")
          VIM::spineMiss(
            x = as.data.frame(wide_data[, mget(c(student.factors, fv))]),
            col = c("#87ceeb", "#ff3333"),
            miss.labels = FALSE, xlab = factor.labels[[student.factors]]
          )
          res.list[[K]][[fv]][[student.factors]] <-
            grDevices::recordPlot(attach = "VIM")
          invisible(grDevices::dev.off())
        }
      }

      wide_data[, Missing := NULL]
    }

    if (output.plots) {
      switch(output.type,
        SVG = outDev <-
          function(plot.path, w = 7, h = 7) {
            svglite::svglite(
              filename = plot.path,
              width = w, height = h,
              bg = "transparent"
            )
          },
        PDF = outDev <-
          function(plot.path, w = 7, h = 7) {
            grDevices::pdf(
              file = plot.path,
              width = w, height = h,
              bg = "transparent"
            )
          }
      )

      cohort.plot.path <-
        file.path(
          diagn.dir,
          paste("Grade", current.grade, subject, current.year, sep = "_")
        )
      if (!dir.exists(cohort.plot.path)) dir.create(cohort.plot.path)

      for (P in c("upset_plot", "missing_plot", "missing_pairs")) {
        plot.nm <-
          file.path(cohort.plot.path, paste0(P, ".", tolower(output.type)))

        outDev(plot.path = plot.nm)
        invisible(print(res.list[[K]][[P]]))
        grDevices::dev.off()
      }

      if (length(student.factors)) {
        factor.paths <- file.path(cohort.plot.path, "factor_plots")
        if (!dir.exists(factor.paths)) dir.create(factor.paths)

        for (sf in student.factors) {
          outDev(
            plot.path =
              file.path(factor.paths, paste0(sf, ".", tolower(output.type)))
          )
          invisible(print(res.list[[K]][[sf]][["factor_plot"]]))
          grDevices::dev.off()
        }

        for (fv in focus.vars) {
          if (!dir.exists(file.path(factor.paths, fv)))
            dir.create(file.path(factor.paths, fv))
          tmp.nms <- grep(paste0("^", student.factors, collapse = "|"),
                          names(res.list[[K]][[fv]]), value = TRUE)
          for (mosaic in tmp.nms) {
            outDev(
              plot.path =
                file.path(factor.paths, fv, paste0(mosaic, ".", tolower(output.type)))
            )
            invisible(print(res.list[[K]][[fv]][[mosaic]]))
            grDevices::dev.off()
          }
        }
      }
    }
  }  ###  END K
  res.list
}
