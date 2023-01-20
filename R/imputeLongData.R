#' @title Multiple imputation of missing values in longitudinal data
#' @description Function produces multiple imputations for missing data
#'  via methods available for the "mice" package.
#' @param long_data The incomplete dataset in which to impute student
#'  (scale score) values.
#' @param return.all.data Return all records from the long data used to
#'  create the imputed scale scores? If TRUE, only records from the current
#'  year are returned.
#'  Default: TRUE
#' @param plot.dir Directory path for diagnostic plots to be placed. Default is
#'  the working directory. Additional subdirectories are created internally as needed.
#'  Default: 'Missing_Data'
#' @param growth.config An elongated `SGP` style config script with an entry for
#'  each grade/content_area/year cohort definition. Configs will be used to subset
#'  the `long_data` provided as required for cohort level data imputation.
#'  Default: NULL
#' @param status.config An elongated `SGP` style config script with an entry for
#'  each grade/content_area/year cohort definition. Configs will be used to subset
#'  the `long_data` provided as required for cohort level data imputation. Unlike
#'  a `growth.config` entry, `status.config` entries use data from the same grade,
#'  but from the prior year(s) (i.e. not individual variables). For example you
#'  might impute missing 3rd grade ELA scores based on a previous year's 3rd grade
#'  school mean scale score, FRL status, etc.
#'  Default: NULL
#' @param focus.variable The variable to be imputed. Currently only a single
#'  variable is allowed. Default: `SCALE_SCORE`
#' @param student.factors Demographic or other student level background information
#'  that will be used in the imputation calculations. The default (NULL) means
#'  that no additional student level factors (beyond observed values of the
#'  specified `focus.variable` - i.e., `SCALE_SCORE`) are included.
#'  Default: NULL
#' @param group Grouping indicator (e.g. institution ID) used to construct group
#'  level means of the `focus.variable` and any `student.factors` provided.
#'  Default: NULL
#' @param missing.group.id If `group` is provided, a placeholder value for
#'  students with no group affiliation (possibly internally after the `long_data`
#'  has been widened and extended).
#'  Default: -99999L
#' @param impute.long Some imputation methods (e.g., "2l.pan") can be used to
#'  impute the data in a longitudinal form (clustered by student over time/grade),
#'  in which case this argument should be set to TRUE. Defaults to FALSE, meaning
#'  scores are imputed in a "wide", or cross-sectional format.
#'  Default: FALSE
#' @param impute.method The name of the method from the `mice` package or other
#'  add-on package functions used to impute the `focus.variable`. The default is
#'  `NULL`, which translates to the default method in the `mice` package - "pmm"
#'  for predictive mean matching. The only other tested method is currently the
#'  `panImpute` functionality from the `mice` and `mitml` packages.
#'  Default: NULL
#' @param parallel.config The default is NULL meaning imputations will be calculated
#'  sequentially in a call to `mice::mice`. Alternatively, a list with named elements
#'  cores, packages and/or cluster.type can be provided. The "cores" element should
#'  be a single numeric value for the number of CPU cores/hyperthreads to use. The
#'  "packages" element is a character string of package name(s) for the requested
#'  imputation methods. "cluster.type" specifies the parallel backend to use (typically
#'  FORK for Linux/Unix/Mac and "PSOCK" for Windows). An example list might look
#'  like: `list(packages = c('mice', 'miceadds'), cores=10)`.
#'  Default: NULL
#' @param seed A random seed set for the imputation process to allow for
#'  replication of results, or for alternative results using the same code.
#'  Default: 4224
#' @param M The number of imputed values to compute.
#'  Default: 15
#' @param maxit The number of iterations allowed for each imputation process.
#'  See the `mice` package documentation for details.
#'  Default: 30
#' @param verbose Defaults to FALSE, meaning progress information from the
#'  `mice` package is not printed out to the console.
#'  Default: FALSE
#' @param ... Additional arguments for the `mice::mice` function and any particular
#'  imputation method/function. See each function/package documentation for details.
#' @return Function returns `long_data` with additional columns of imputed values.
#' @details The function returns the dataset supplied to the `long_data` argument
#'  along with `M` additional columns populated with the imputed values.
#' @examples 
#' 	\dontrun{
#'     data_to_impute <- SGPdata::sgpData_LONG_COVID
#'
#'     ###   Read in STEP 0 SGP configuration scripts
#'     source("SGP_CONFIG/STEP_0/Ampute_2021/Growth.R")
#'     source("SGP_CONFIG/STEP_0/Ampute_2021/Status.R")
#'     Test_Data_LONG <-
#'         imputeScaleScore(
#'             long_data = data_to_impute,
#'             growth.config = growth_config_2021,
#'             status.config = status_config_2021,
#'             M = 1)
#'  }
#' @seealso 
#'  \code{\link[mice]{mice}}, \code{\link[mitml]{panImpute}}
#' @rdname imputeLongData
#' @export 
#' @author Adam R. VanIwaarden \email{avaniwaarden@nciea.org}
#' @keywords misc
#' @keywords models
#' @importFrom parallel detectCores
#' @importFrom data.table SJ setnames setkey setkeyv dcast melt data.table rbindlist as.data.table key
#' @importFrom stats setNames
#' @importFrom mice make.blocks make.predictorMatrix mice densityplot complete

imputeLongData =
  function(
    long_data,
    return.all.data = TRUE,
    plot.dir = "Missing_Data",
    status.config = NULL,
    growth.config = NULL,
    focus.variable = "SCALE_SCORE",
    student.factors = NULL,
    group = NULL,
    missing.group.id = -99999L,
    impute.long = FALSE,
    impute.method = NULL,
    parallel.config = NULL, # define cores, packages, cluster.type
    seed = 4224L,
    M = 15,
    maxit = 30,
    verbose = FALSE,
    ...
  ) {

  call <- match.call()

  ###  Avoid "Undefined global functions or variables:" from R CMD check
  ID <- VALID_CASE <- YEAR <- CONTENT_AREA <- GRADE <- .imp <- ALL_NA <- NULL

  ###  Initial checks
  if (!dir.exists(diagn.dir <- file.path(plot.dir, "imputeLongData", "diagnostic_plots"))) {
      dir.create(diagn.dir, recursive = TRUE)
  }

  if (length(group) > 1) {
    stop("Only one 'group' variable can be used as an imputation factor.")
  }

  if (is.null(impute.method)) {
    impute.method <- "pmm"
  }

  if (is.null(parallel.config)) {
    parexecute <- "SEQ"
  } else {
    parexecute <- "PAR"
    if (is.logical(parallel.config)) {
      parallel.config <- list()
    }
    if (is.null(parallel.config$packages)) {
      parallel.config$packages <- "mice"
    }
    if (is.null(parallel.config$cores)) {
      parallel.config$cores <- parallel::detectCores() - 1L
    }
    if (is.null(parallel.config$cluster.type)) {
      parallel.config$cluster.type <-
        ifelse(.Platform$OS.type == "unix", "FORK", "PSOCK")
    }
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

  long.to.wide.vars <-
    unique(c(focus.variable, group, student.factors)) # "CONTENT_AREA", "GRADE",

  ###  Cycle through configs to get results by cohort
  res.list <- vector(mode = "list", length = length(configs))

  for (K in seq(configs)) {
    cohort.iter <- configs[[K]]
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
      long_data[cohort.lookup][,
        c(getKey(long_data), long.to.wide.vars), with = FALSE]

    prior.years <- utils::head(unique(cohort.iter[["panel.years"]]), -1)
    current.year <- utils::tail(unique(cohort.iter[["panel.years"]]), 1)
    current.grade <- utils::tail(unique(cohort.iter[["grade.sequences"]]), 1)
    prior.grades <- utils::head(unique(cohort.iter[["grade.sequences"]]), -1)

    if (cohort.iter$analysis.type == "GROWTH") {
      ###  convert long to wide
      wide_data <-
        data.table::dcast(
          data = tmp_long,
          formula = ID ~ GRADE + CONTENT_AREA,
          sep = ".", drop = FALSE,
          value.var = c("VALID_CASE", long.to.wide.vars)
        )

      ###  Exclude kids missing IN current AND most recent year/grade's data
      ###  (keeps them if in data with missing score, but set to "VALID_CASE")
      tmp.vcase <- paste0(paste0("VALID_CASE.", c(utils::tail(prior.grades, 1), current.grade)),
                          ".", rep(unique(cohort.iter$content.areas), each = 2))
      excl.idx <-
        Reduce(intersect,
               lapply(tmp.vcase, \(f) which(is.na(wide_data[, get(f)])))
        )
      if (length(excl.idx)) {
        wide_data <- wide_data[-excl.idx]
      }

      #  remove columns that are all NA (e.g., SGP for 3rd grade priors)
      na.vars <-
        unlist(
          lapply(
            names(wide_data),
            \(f) all(is.na(wide_data[, get(f)]))
          )
        )
      if (any(na.vars)) {
        wide_data <- wide_data[, names(wide_data)[!na.vars], with = FALSE]
      }

      meas.list <- vector(mode = "list", length = length(long.to.wide.vars))
      meas.list <-
        lapply(
          long.to.wide.vars,
          \(f) {meas.list[[f]] <- grep(paste0("^", f, "[.]"), names(wide_data))}
        )
      names(meas.list) <- long.to.wide.vars

      long_final <- data.table::melt(data = wide_data,
                                     id.vars = "ID",
                                     variable.name = "GRADE",
                                     measure.vars = meas.list)
      long_final[,
        CONTENT_AREA :=
          as.character(factor(GRADE, labels = cohort.lookup[["CONTENT_AREA"]]))
      ][,
        YEAR :=
          as.character(factor(GRADE, labels = cohort.lookup[["YEAR"]]))
      ][,
        GRADE :=
          as.character(factor(GRADE, labels = cohort.lookup[["GRADE"]]))
      ]

      if (impute.long) {
        tmp.focus.variable <- focus.variable
        ##  Convert character/factor to numeric (0/1) data
        student.vars <-
          grep(
            paste(student.factors, collapse = "|"),
            focus.variable,
            value = TRUE
          )
        if (length(student.vars)) {
          for (s.var in student.vars) {
            long_final[, eval(s.var) := as.integer(factor(get(s.var))) - 1L]
            if (length(group)) {
              tmp.imv <- paste0("IMV___", s.var, "___", group)
              long_final[,
                eval(tmp.imv) := mean(get(s.var), na.rm = TRUE), by = c("GRADE", group)
              ]
              long_final[,
                eval(tmp.imv) := scale(get(tmp.imv)), keyby = "GRADE"
              ]
              long_final[is.na(get(tmp.imv)), eval(tmp.imv) := 0]
              tmp.focus.variable <- c(tmp.focus.variable, tmp.imv)
            }
          }
        }

        if (length(group)) {
          tmp.imv <- paste0("IMV___", focus.variable, "___", group)
          long_final[, eval(tmp.imv) :=
            mean(get(focus.variable), na.rm = TRUE), by = c("GRADE", group)]
          long_final[, eval(tmp.imv) :=
            scale(get(tmp.imv)), keyby = "GRADE"]
          long_final[is.na(get(tmp.imv)), eval(tmp.imv) := 0]
          tmp.focus.variable <- c(tmp.focus.variable, tmp.imv)
        }

        if (impute.method %in% c("2l.pan", "2l.lmer", "2l.norm")) {
          long_final[, ID := as.integer(ID)]
        } else {
          long_final[, ID := as.character(ID)]
        }
        long_final[, GRADE := as.numeric(GRADE)]
      }
      long_final[, VALID_CASE := "VALID_CASE"]

    } else {  #  END "GROWTH"  --  Begin "STATUS"
      ##  Create wide_data and tmp_long_priors
      status.lookup <- cohort.lookup[YEAR == current.year]
      data.table::setkeyv(tmp_long, getKey(tmp_long))

      wide.fmla <-
        stats::as.formula(
          paste("ID", ifelse(is.null(group), "", paste("+", group)),
                "~ GRADE + CONTENT_AREA"
          ))
      wide_data <-
        data.table::dcast(
          data = tmp_long[status.lookup],
          formula = wide.fmla,
          sep = ".", drop = FALSE,
          value.var = c(focus.variable, student.factors)
        )

      if (!is.null(group)) {
        group.lookup <-
          unique(tmp_long[status.lookup, c("ID", group), with = FALSE])
        wide_data <- wide_data[group.lookup]
      }

      if (length(student.factors)) {
        subset.vars <- names(wide_data) %w/o% "ID" # group?
        focus.vars <- grep(focus.variable, subset.vars, value = TRUE)
        student.vars <- subset.vars %w/o% c(group, focus.vars)
        for (sv in student.vars) {
          wide_data[is.na(get(sv)), eval(sv) := 0]
        }
      } else {
        tmp.wide.names <- grep(current.grade, names(wide_data), value = TRUE)
        subset.vars <- focus.vars <- paste0(focus.variable, ".", tmp.wide.names)
        data.table::setnames(wide_data, tmp.wide.names, focus.vars)
      }

      #  remove columns that are all NA (e.g., SGP for 3rd grade priors)
      na.vars <-
        unlist(
          lapply(
            names(wide_data),
            \(f) all(is.na(wide_data[, get(f)]))
          )
        )
      if (any(na.vars)) {
        na.var.nms <- names(wide_data)[na.vars]
        focus.vars <- focus.vars %w/o% na.var.nms
        student.factors <- student.factors %w/o% na.var.nms
        group.vars <- group.vars %w/o% na.var.nms
        wide_data <- wide_data[, names(wide_data)[!na.vars], with = FALSE]
      }

      if (length(group)) {
        group.vars <- paste0(group, ".", current.grade)
        priors.lookup <- cohort.lookup[YEAR != current.year]
        tmp_long_priors <-
          tmp_long[priors.lookup][,
            c("YEAR", "CONTENT_AREA", "GRADE", focus.variable, group),
            with = FALSE
          ]
        ##  Prior grade(s) summaries to use
        smry.eval.expression <-
          paste0("PRIOR_IMV__", focus.variable,
                 " = ", "mean(", focus.variable, ", na.rm=TRUE)")
        smry.eval.expression <-
          stats::setNames(
            smry.eval.expression,
            sub("^(.*) = .*", "\\1", smry.eval.expression)
          )

        tmp_grp_smry <-
          tmp_long_priors[
            !is.na(get(group)),
              lapply(smry.eval.expression, \(f) eval(parse(text = f))),
            keyby = c("GRADE", "CONTENT_AREA", group)
          ]

        tmp_grp_smry <-
          data.table::dcast(
            data = tmp_grp_smry,
            formula = get(group) ~ GRADE + CONTENT_AREA,
            sep = ".",
            value.var = names(smry.eval.expression)
          )
        data.table::setnames(
          tmp_grp_smry,
          c(group, paste0(names(smry.eval.expression), ".", names(tmp_grp_smry)[-1]))
        )

        wide_data <- merge(wide_data, tmp_grp_smry, by = group, all.x = TRUE)
        data.table::setnames(wide_data, group, group.vars)

        ##  Put in cross school mean for schools with no students in prior years
        for (prior.smry in grep("PRIOR_IMV__", names(wide_data), value = TRUE)) {
          tmp.grp.mean <- mean(wide_data[, get(prior.smry)], na.rm = TRUE)
          wide_data[is.na(get(prior.smry)), eval(prior.smry) := tmp.grp.mean]
        }
      }

      long_final <-
        long_data[status.lookup][,
          unique(c(key(long_data), long.to.wide.vars)), with = FALSE]
    } ###  END "STATUS"

    ##  Subset out scale scores and student background (factors)
    if (!impute.long && cohort.iter$analysis.type == "GROWTH") {  #  done above for "STATUS"
      subset.vars <- focus.vars <-
        grep(paste0(focus.variable, "[.]", collapse = "|"), names(wide_data), value = TRUE)
      if (length(group)) {
        group.vars <- grep(group, names(wide_data), value = TRUE)
        subset.vars <- c(subset.vars, group.vars)
      }
      if (length(student.factors)) {
        student.vars <- grep(paste0(student.factors, "[.]", collapse = "|"),
                             names(wide_data), value = TRUE)
        subset.vars <- c(subset.vars, student.vars)
        for (sv in student.vars) {
          wide_data[is.na(get(sv)), eval(sv) := 0]
        }
      }

      #  create subset of wide data with only variables to be used in imputation
      wide_data <-
        wide_data[,
          grep(
            paste0("^ID$|", paste(subset.vars, collapse = "|")),
            names(wide_data)
          ),
          with = FALSE
        ]

      # add in prior group level summaries if requested
      if (!is.null(cohort.iter[["status.years"]]) && length(group)) {
        yr.length <- length(cohort.iter[["status.years"]])
        unq.ca <- unique(cohort.iter[["content.areas"]])
        priors.lookup <-
          data.table::SJ(
            "VALID_CASE",
            rep(unq.ca, yr.length),
            rep(cohort.iter[["status.years"]], each = length(unq.ca)),
            rep(cohort.iter[["grade.sequences"]], each = length(unq.ca))
          ) |>
          data.table::setnames(getKey(long_data) %w/o% "ID") |>
            data.table::setkey(YEAR) |> # ensure lookup table is ordered by years.
              data.table::setkey(NULL)  # NULL out key so doesn't corrupt the join in dcast.

        data.table::setkeyv(long_data, getKey(long_data))

        tmp_long_priors <-
          long_data[priors.lookup][,
            c("YEAR", "CONTENT_AREA", "GRADE", focus.variable, group),
            with = FALSE
          ]
        ##  Prior grade(s) summaries to use
        smry.eval.expression <-
          paste0("PRIOR_IMV__", focus.variable,
                 " = ", "mean(", focus.variable, ", na.rm=TRUE)")
        smry.eval.expression <-
          stats::setNames(
            smry.eval.expression,
            sub("^(.*) = .*", "\\1", smry.eval.expression)
          )

        tmp_grp_smry <-
          tmp_long_priors[
            !is.na(get(group)),
              lapply(smry.eval.expression, \(f) eval(parse(text = f))),
            keyby = c("GRADE", "CONTENT_AREA", group) # "YEAR", 
          ]

        tmp_grp_smry <-
          data.table::dcast(
            data = tmp_grp_smry,
            formula = get(group) ~ GRADE + CONTENT_AREA, # YEAR + 
            sep = ".",
            value.var = names(smry.eval.expression)
          )

        tmp.grd.subj <- names(tmp_grp_smry)[-1]
        ##  Put in cross school mean for schools with no students in prior years
        for (prior.smry in tmp.grd.subj) {
          smry_names <-
            c(group, paste0("PRIOR_IMV__", focus.variable)) |>
              paste(prior.smry, sep = ".")
          if (!smry_names[1] %in% names(wide_data)) next
          sub_grp_smry <- tmp_grp_smry[, c("group", prior.smry), with = FALSE]
          setnames(sub_grp_smry, smry_names)

          wide_data <-
            merge(wide_data, sub_grp_smry, by = smry_names[1], all.x = TRUE)

          tmp.grp.mean <- mean(wide_data[, get(smry_names[2])], na.rm = TRUE)
          wide_data[
            is.na(get(smry_names[2])),
            eval(smry_names[2]) := tmp.grp.mean
          ]
        }
      }

      #  remove columns that are all NA (e.g., SGP for 3rd grade priors)
      na.vars <-
        lapply(
          names(wide_data),
          \(f) all(is.na(wide_data[, get(f)]))
        ) |> unlist()
      if (any(na.vars)) {
        na.var.nms <- names(wide_data)[na.vars]
        focus.vars <- focus.vars %w/o% na.var.nms
        student.factors <- student.factors %w/o% na.var.nms
        group.vars <- group.vars %w/o% na.var.nms
        wide_data <- wide_data[, names(wide_data)[!na.vars], with = FALSE]
      }
    }  ###  END "GROWTH"

    #####
    ###   IMPUTE
    #####

    if (!impute.long || cohort.iter$analysis.type == "STATUS") {
      ##  Final (if necessary) whittling down of the data
      if (length(group) & toupper(impute.method) == "PANIMPUTE") {
        # init.where <- tmp.where <- mice::make.where(data = wide_data) # init.where <-
        na_removed_records <- data.table::data.table()

        for (gv in group.vars) { # [3:6]
          # if (cluster.by.group) {
          if (cohort.iter$analysis.type == "GROWTH") {
            tmp.imp.var <- paste0(focus.variable, ".", sub(".+?[.]", "", gv))
            tmp.fixd.ef <- focus.vars %w/o% tmp.imp.var

            invisible(wide_data[is.na(get(gv)), eval(gv) := missing.group.id])

            ## Check for clusters with all NA values
            na_check <-
              wide_data[, list(ALL_NA = all(is.na(get(tmp.imp.var)))), keyby = gv]

            if (length(all.na.ids <- na_check[ALL_NA == TRUE, get(gv)])) {
              na_removed_records <-
                data.table::rbindlist(
                  list(
                    na_removed_records,
                    wide_data[get(gv) %in% all.na.ids, ]
                  )
                )
              wide_data <- wide_data[!(get(gv) %in% all.na.ids), ]
            }
          } else {
            invisible(wide_data[is.na(get(gv)), eval(gv) := missing.group.id])

            ## Check for clusters with all NA values
            na.check.expression <-
              paste0("ALL_NA_", focus.vars, # "_", group,
                    " = ", "all(is.na(", focus.vars, "))"
              )
            na.check.expression <-
              stats::setNames(
                na.check.expression,
                sub("^(.*) = .*", "\\1", na.check.expression)
              )
            na_check <-
              wide_data[,
                lapply(na.check.expression, \(f) eval(parse(text = f))),
                keyby = gv
              ]

            for (nac in names(na.check.expression)) {
              if (length(all.na.ids <- na_check[get(nac) == TRUE, get(gv)])) {
                na_removed_records <-
                  data.table::rbindlist(
                    list(
                      na_removed_records,
                      wide_data[get(gv) %in% all.na.ids, ]
                    )
                  )
                wide_data <- wide_data[!(get(gv) %in% all.na.ids), ]
              }
            }
          }
        }
      } else {
        na_removed_records <- NULL
      }

      #  re-check & remove columns that are all identical
      #  (causes error in pan::pan - 'NA/NaN/Inf in foreign function call (arg 20)')
      same.same <-
        lapply(
          names(wide_data),
          \(f) length(unique(wide_data[, get(f)])) == 1 # would work for NA's above too (`na.vars`)
        ) |> unlist()
      if (any(same.same)) {
        na.var.nms <- names(wide_data)[same.same]
        focus.vars <- focus.vars %w/o% na.var.nms
        # student.factors <- student.factors %w/o% na.var.nms
        group.vars <- group.vars %w/o% na.var.nms
        wide_data <- wide_data[, names(wide_data)[!same.same], with = FALSE]
        if (length(na_removed_records)) {
          na_removed_records <-
            na_removed_records[, names(na_removed_records)[!same.same], with = FALSE]
        }
      }

      ##  arguments/objects for `mice`
      tmp.meth <-
        rep(impute.method, length(focus.vars)) |>
          stats::setNames(focus.vars)
      tmp.pred <- names(wide_data) %w/o% "ID"

          ###   `mice`/... argument alignment
      if (!"blocks" %in% names(call)) {
        tmp.blok <-
          mice::make.blocks(
            data = focus.vars,
            calltype = "formula"
          )
      } else {
        tmp.blok <- call[["blocks"]]
      }

      if (!"formulas" %in% names(call)) {
        tmp.fmla <- list()
        for (fv in focus.vars) {
          tmp.iv <- focus.vars %w/o% fv# tmp.pred[fv, tmp.pred[fv, ] != 0]
          tmp.grd.subj <-
            gsub(paste0(focus.variable, "."), "", fv)

          if (length(student.factors)) {
            tmp.sf <-
              grep(
                paste0(student.factors, ".", tmp.grd.subj, collapse = "|"),
                tmp.pred,
                value = TRUE
              )
          } else {
            tmp.sf <- NULL
          }

          if (length(group)) {
            tmp.gv <- grep(tmp.grd.subj, group.vars, value = TRUE)
            if (!length(tmp.gv)) tmp.gv <- group.vars # for STATUS configs

            # SCALE_SCORE.X.Y.Z + I(clusterMeans(SCALE_SCORE.X.Y.Z, group)) + (1|group)
            rhs <- paste(tmp.iv, collapse = " + ")

            if (length(tmp.pr <- grep("PRIOR_IMV__", tmp.pred))) {
              rhs <-
                c(rhs,
                  paste0(tmp.pred[tmp.pr], collapse = " + ")
                )
            }
            if (length(tmp.sf)) {
              rhs <-
                c(rhs,
                  ifelse(toupper(impute.method) == "PANIMPUTE",
                    paste0(tmp.sf, ":I(clusterMeans(", tmp.sf, ", ", tmp.gv, "))", collapse = " + "),
                    paste0(tmp.sf, collapse = " + ")
                  )
                )
            }

            if (length(tmp.iv) && toupper(impute.method) == "PANIMPUTE") { #  && cohort.iter$analysis.type == "STATUS"
              rhs <-
                c(rhs,
                  paste0("I(clusterMeans(", tmp.iv, ", ", tmp.gv, "))", collapse = " + ")
                )
              rslope <- paste(" +", tmp.iv, collapse = "") # NULL # 
            } else {
              rslope <- NULL
            }

            if (toupper(impute.method) == "PANIMPUTE") { #  && cohort.iter$analysis.type == "STATUS"
              rhs <- c(rhs, paste0("(1", rslope, "|", tmp.gv, ")"))
            }

            tmp.fmla[[fv]] <-
              paste(fv, "~", paste(rhs, collapse = " + ")) |>
                stats::as.formula()
          } else {
            tmp.fmla[[fv]] <-
              paste(fv, "~", paste(c(tmp.iv, tmp.sf), collapse = " + ")) |> stats::as.formula()
          }
        }
      } else {
        tmp.fmla <- call[["formulas"]]
      }
    } else { # impute.long
      wide_data <-
        long_final[,
          c("ID", "GRADE", tmp.focus.variable %w/o% group),
          with = FALSE
        ]
      na_removed_records <- NULL

      if (is.null(impute.method)) {
        tmp.meth <- "2l.pan"
      } else {
        tmp.meth <- impute.method
      }
      tmp.meth <- stats::setNames(tmp.meth, focus.variable)

      tmp.pred <- mice::make.predictorMatrix(data = wide_data)
      tmp.pred[focus.variable, "ID"] <- -2
      tmp.pred[focus.variable, "GRADE"] <- 2 # random effect for GRADE (time)
    }

    if (length(na_removed_records)) {
      message(
        paste0("\n\t\t", nrow(na_removed_records),
              " student(s) will not be imputed with `", impute.method, "`.
                All records in their `group` are missing for one or more years.
                Imputations will be attempted without `group` via `pmm` separately.\n"
        )  #  Spacing above should align text to left (with tab padding).
      )
    }

    switch(parexecute,
      SEQ = imp <-
        suppressWarnings(
          mice::mice(
            data = wide_data, # where = tmp.where, ignore = tmp.ign, # can't use `where` or `ignore` with panImpute
            m = M, method = tmp.meth,
            visitSequence = focus.vars,
            blocks = tmp.blok, formulas = tmp.fmla,
            maxit = maxit, seed = seed, print = verbose,
            ...
        )),
      PAR = imp <-
        suppressWarnings(
          parMICE(
            data = wide_data, m = M, maxit = maxit, meth = tmp.meth,
            bloks = tmp.blok, frmlas = tmp.fmla, visit = focus.vars,
            seed = seed, nnodes = parallel.config$cores,
            cluster.type = parallel.config$cluster.type,
            packages = parallel.config$packages
        ))
      )

    ##  Save some diagnostic plots
    if (cohort.iter$analysis.type == "GROWTH") {
      imp.type <- "GROWTH_" # paste0(current.subject, "_")
    } else {
      imp.type <- "STATUS_"
    }

    grDevices::pdf(
      file = 
        file.path(
          diagn.dir,
          paste0("Grade_", current.grade, "_", current.year, "_", imp.type,
                 gsub("[.]", "", impute.method),
                 ifelse(impute.long, "_LONG", ""),
                 "_M_", M, "__maxit_", maxit,
                 ifelse(is.null(group), "", paste0("_x_", group)),
                 "__converge.pdf"
          )
        )
    )
    print(graphics::plot(imp))
    invisible(grDevices::dev.off())

    grDevices::pdf(
      file =
        file.path(
          diagn.dir,
          paste0("Grade_", current.grade, "_", current.year, "_", imp.type,
                 gsub("[.]", "", impute.method),
                 ifelse(impute.long, "_LONG", ""),
                 "_M_", M, "__maxit_", maxit,
                 ifelse(is.null(group), "", paste0("_x_", group)),
                 "__density.pdf"
          )
        ),
      width = 11, height = 8
    )
    if (!impute.long) { # && cohort.iter$analysis.type == "GROWTH"
      tryCatch(
        print(mice::densityplot(
          imp,
          eval(parse(text = paste0("~", paste(focus.vars, collapse = " + "))))
        )),
        error = function(e) TRUE
      ) -> err.tf
      if (is.logical(err.tf)) print(mice::densityplot(imp)) else rm(err.tf)
    } else {
      print(mice::densityplot(imp))
    }
    invisible(grDevices::dev.off())

    ###  Format and store results
    if (!impute.long || cohort.iter$analysis.type == "STATUS") {
      long_imputed <-
        data.table::as.data.table(
          mice::complete(imp, action = "long", include = TRUE)
        )
      if (cohort.iter$analysis.type == "STATUS") {
        long.ids <- c("ID", ".imp", gv)
      } else {
        long.ids <- c("ID", ".imp")
      }
      wide_imputed <-
        data.table::melt(
          data = long_imputed,
          id.vars = long.ids,
          value.name = focus.variable,
          variable.name = "GRADE",
          measure.vars = focus.vars
        )

      if (length(na_removed_records)) {
        tmp_imputed <-
          rbindlist(
            list(
              na_removed_records[, NA_RMVD := TRUE],
              data.table::as.data.table(
                mice::complete(imp, action = M, include = FALSE)
              )[, NA_RMVD := FALSE]
            )
          )
        if (cohort.iter$analysis.type != "STATUS") {
          tmp_imputed[, grep(group, names(tmp_imputed)) := NULL]
        }
        tmp.meth2 <-
        rep("", ncol(tmp_imputed)) |>
          stats::setNames(names(tmp_imputed))
        tmp.meth2[focus.vars] <- "pmm"
        imp2 <- #mice::mice(tmp_imputed)
          suppressWarnings(
            mice::mice(
              data = tmp_imputed,
              m = M, method = tmp.meth2,
              visitSequence = focus.vars,
              maxit = maxit, seed = seed, print = verbose
            )
          )
        long_impute2 <-
          data.table::as.data.table(
            mice::complete(imp2, action = "long", include = TRUE)
          )[NA_RMVD == TRUE][, NA_RMVD := NULL]
        wide_impute2 <-
          data.table::melt(
            data = long_impute2,
            id.vars = long.ids,
            value.name = focus.variable,
            variable.name = "GRADE",
            measure.vars = focus.vars
          )
        wide_imputed <-
          data.table::rbindlist(
            list(wide_impute2[, IMPUTE_METHOD := "pmm"],
                 wide_imputed[, IMPUTE_METHOD := impute.method]
            )
          )
      } else {
        invisible(wide_imputed[, IMPUTE_METHOD := impute.method])
      }

      if (cohort.iter$analysis.type == "STATUS") {
        wide_imputed[,
          CONTENT_AREA := gsub(paste0(".*", current.grade, "[.]"), "", GRADE)
        ][,
          GRADE := current.grade
        ][,
          YEAR := current.year
        ]
      } else {
        wide_imputed[, GRADE := gsub(paste0(focus.variable, "."), "", GRADE)]
        unq.ca <- unique(cohort.iter[["content.areas"]])

        wide_imputed[,
          CONTENT_AREA :=
            gsub(paste(paste0(".*", c(prior.grades, current.grade), "."), collapse = "|"), "", GRADE)
          # gsub(paste(paste0(".*.", c(prior.grades, current.grade), "."), collapse = "|"), "", GRADE)
        ][,
          GRADE := #as.character(factor(GRADE, labels = cohort.iter[["grade.sequences"]]))
            gsub(paste(paste0(".", unq.ca), collapse = "|"), "", GRADE)
        ][,
          YEAR :=
            as.character(
              factor(GRADE,
                levels = cohort.iter[["grade.sequences"]],
                labels = cohort.iter[["panel.years"]]
            ))
        ]
      }
        
      wide_imputed <-
        data.table::dcast(
          data = wide_imputed,
          formula = as.formula(
            paste(
              paste(long.ids %w/o% ".imp", collapse = " + "),
                    "+ CONTENT_AREA + GRADE + YEAR + IMPUTE_METHOD ~ .imp"
            )),
          value.var = focus.variable
        )
      if (cohort.iter$analysis.type == "STATUS") {
        rekey <- key(wide_imputed) %w/o% gv
        wide_imputed[, eval(gv) := NULL]
        setkeyv(wide_imputed, rekey)
      }
    } else {
      long_imputed <-
        data.table::as.data.table(
          mice::complete(imp, action = "long", include = TRUE)
        )[,
          c(".imp", "ID", "YEAR", "GRADE", focus.variable),
          with = FALSE
        ]

      wide_imputed <-
        data.table::dcast(
          data = long_imputed[GRADE == current.grade, ],
          formula = ID ~ .imp,
          value.var = focus.variable
        )
    }

    data.table::setnames(
      wide_imputed,
      as.character(0:M),
      c(focus.variable, paste0(focus.variable, ".IMP_", 1:M))
    )
    tmp.key <-
      c(data.table::key(wide_imputed), focus.variable) %w/o% "IMPUTE_METHOD" # "ID", ..., "YEAR"
    data.table::setkeyv(wide_imputed, tmp.key)
    data.table::setkeyv(long_final, tmp.key)
    if (nrow(wide_imputed) >= nrow(long_final)) {
      res.list[[K]] <- long_final[wide_imputed][, VALID_CASE := "VALID_CASE"]
    } else {
      res.list[[K]] <- wide_imputed[long_final]
    }
    # res.list[[K]] <- wide_imputed[long_final]
    rm(list = c("long_final", "wide_data", "wide_imputed", "imp"))
    invisible(gc())
  }  ###  END K

  ###   Compile and format final imputed data set
  final_imp_data <-
    data.table::rbindlist(
      res.list,
      use.names = TRUE,
      fill = TRUE
    )

  ##  remove dups for repeaters created from long to wide to long reshaping
  data.table::setkeyv(final_imp_data, getKey(final_imp_data))
  data.table::setkeyv(final_imp_data, data.table::key(final_imp_data) %w/o% "GRADE")
  dup.ids <-
    final_imp_data[
      which(duplicated(final_imp_data, by = data.table::key(final_imp_data))),
      ID
    ]
  if (length(dup.ids)) {
    final_imp_data[
      ID %in% dup.ids &
      is.na(get(focus.variable)),
        VALID_CASE := "NEW_DUP"
    ]

    final_imp_data <- final_imp_data[VALID_CASE != "NEW_DUP"]
  }

  if (return.all.data) {
    ##  Remove group, student.factors
    if (!is.null(group) || !is.null(student.factors)) {
      final_imp_data[, unique(c(group, student.factors)) := NULL]
    }

    fin.key <- c(getKey(final_imp_data), focus.variable)
    data.table::setkeyv(long_data, fin.key)
    data.table::setkeyv(final_imp_data, fin.key)

    final_imp_data <- final_imp_data[long_data]
  }
  return(final_imp_data)
}
