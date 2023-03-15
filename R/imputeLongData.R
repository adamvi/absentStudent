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
#' @param status.config An elongated `SGP` style config script with an entry for
#'  each grade/content_area/year cohort definition. Configs will be used to subset
#'  the `long_data` provided as required for cohort level data imputation. Unlike
#'  a `growth.config` entry, `status.config` entries use data from the same grade,
#'  but from the prior year(s) (i.e. not individual variables). For example you
#'  might impute missing 3rd grade ELA scores based on a previous year's 3rd grade
#'  school mean scale score, FRL status, etc.
#'  Default: NULL
#' @param growth.config An elongated `SGP` style config script with an entry for
#'  each grade/content_area/year cohort definition. Configs will be used to subset
#'  the `long_data` provided as required for cohort level data imputation.
#'  Default: NULL
#' @param focus.variable The variable to be imputed. Currently only a single
#'  variable is allowed. Default: `SCALE_SCORE`
#' @param focus.years Values of YEAR that should use the `impute.method`.
#'  If provided, `focus.variable` values missing in any excluded year(s) will be
#'  imputed using the (typically faster) `pmm` method.
#'  Default: NULL
#' @param student.factors Demographic or other student level background information
#'  that will be used in the imputation calculations. The default (NULL) means
#'  that no additional student level factors (beyond observed values of the
#'  specified `focus.variable` - i.e., `SCALE_SCORE`) are included.
#'  Default: NULL
#' @param group Grouping indicator (e.g. institution ID) used to construct group
#'  level means of the `focus.variable` and any `student.factors` provided.
#'  Default: NULL
#' @param impute.method The name of the method from the `mice` package or other
#'  add-on package functions used to impute the `focus.variable`. The default is
#'  `NULL`, which translates to the default method in the `mice` package - "pmm"
#'  for predictive mean matching. The only other tested method is currently the
#'  `panImpute` functionality from the `mice` and `mitml` packages.
#'  Default: NULL
#' @param formula.specs Imputation uses the the `formula` argument (and therefore
#'  only those imputation methods that use this argument can be used). This function
#'  determines those formula's dependent upon the supplied configs and data. This
#'  argument allows for the modification (simpler/more complex, whether to include
#'  values of the `focus.variable` formula as random slopes where applicable with
#'  `group`, etc.).
#'  Default: list(simple = TRUE, random.slope = FALSE)
#' @param parallel.config The default is NULL meaning imputations will be calculated
#'  sequentially in a call to `mice::mice`. Alternatively, a list with named elements
#'  cores, packages and/or cluster.type can be provided. The "cores" element should
#'  be a single numeric value for the number of CPU cores/hyperthreads to use. The
#'  "packages" element is a character string of package name(s) for the requested
#'  imputation methods. "cluster.type" specifies the parallel backend to use
#'  (typically FORK for Linux/Unix/Mac and "PSOCK" for Windows). An example list
#'  might look like: `list(packages = c('mice', 'miceadds'), cores=10)`.
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
#'     source("SGP_CONFIG/STEP_0/Impute_2023/Growth.R")
#'     source("SGP_CONFIG/STEP_0/Impute_2023/Status.R")
#'     Test_Data_LONG <-
#'         imputeScaleScore(
#'             long_data = data_to_impute,
#'             growth.config = growth_config_2023,
#'             status.config = status_config_2023,
#'             M = 1)
#'  }
#' @seealso 
#'  \code{\link[parallel]{detectCores}}
#'  \code{\link[data.table]{data.table-package}}, \code{\link[data.table]{setkey}}, \code{\link[data.table]{dcast.data.table}}, \code{\link[data.table]{setattr}}, \code{\link[data.table]{as.data.table}}, \code{\link[data.table]{melt.data.table}}, \code{\link[data.table]{setcolorder}}
#'  \code{\link[utils]{head}}
#'  \code{\link[arrow]{Table}}, \code{\link[arrow]{write_dataset}}, \code{\link[arrow]{open_dataset}}
#'  \code{\link[stats]{formula}}, \code{\link[stats]{setNames}}
#'  \code{\link[mice]{make.blocks}}, \code{\link[mice]{mice}}, \code{\link[mice]{densityplot.mids}}, \code{\link[mice]{complete.mids}}
#'  \code{\link[grDevices]{pdf}}, \code{\link[grDevices]{dev}}
#'  \code{\link[graphics]{plot.default}}
#'  \code{\link[dplyr]{filter-joins}}
#' @rdname imputeLongData
#' @export 
#' @importFrom parallel detectCores
#' @importFrom data.table data.table setkey setkeyv dcast setnames as.data.table melt setcolorder key
#' @importFrom utils tail head
#' @importFrom arrow arrow_table write_dataset open_dataset
#' @importFrom stats as.formula setNames
#' @importFrom mice make.blocks mice densityplot complete
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics plot
#' @importFrom dplyr semi_join
imputeLongData =
  function(
    long_data,
    return.all.data = TRUE,
    plot.dir = "Missing_Data",
    status.config = NULL,
    growth.config = NULL,
    focus.variable = "SCALE_SCORE",
    focus.years = NULL,
    student.factors = NULL,
    group = NULL,
    impute.method = NULL,
    formula.specs = list(simple = TRUE, random.slope = FALSE),
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
  diagn.dir <- file.path(plot.dir, "imputeLongData", "diagnostic_plots")
  if (!dir.exists(diagn.dir)) {
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
    unique(c(focus.variable, group, student.factors))

  ###  Cycle through configs to get results by cohort
  # res.list <- vector(mode = "list", length = length(configs))

  message("\n\tImputation with `imputeLongData` started: ", date())

  started.impute <- proc.time()

  for (K in seq(configs)) {
    cohort.iter <- configs[[K]]
    names(cohort.iter) <- gsub("^sgp[.]", "", names(cohort.iter))
    grade.length <- length(cohort.iter[["grade.sequences"]])

    cohort_lookup <-
      data.table::data.table(
        VALID_CASE = "VALID_CASE",
        CONTENT_AREA =
            utils::tail(cohort.iter[["content.areas"]], grade.length),
        YEAR = utils::tail(cohort.iter[["panel.years"]], grade.length),
        GRADE = cohort.iter[["grade.sequences"]]
      ) |>
        data.table::setkey(YEAR) |> # ensure lookup table is ordered by years.
          data.table::setkey(NULL)  # NULL out key so doesn't corrupt the join in dcast.

    if (!any(grepl("arrow", class(long_data), ignore.case = TRUE))) {
      arrow.tf <- FALSE
      data.table::setkeyv(long_data, getKey(long_data))
    } else {
      arrow.tf <- TRUE
      cohort_lookup_arw <- arrow::arrow_table(cohort_lookup)
      cohort_lookup_arw <-
        cohort_lookup_arw$cast(
          target_schema =
            long_data$schema[c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE")]
        )
    }

    prior.years <- utils::head(unique(cohort.iter[["panel.years"]]), -1)
    current.year <- utils::tail(unique(cohort.iter[["panel.years"]]), 1)
    current.grade <- utils::tail(unique(cohort.iter[["grade.sequences"]]), 1)
    prior.grades <- utils::head(unique(cohort.iter[["grade.sequences"]]), -1)

    sample.size <- cohort.iter[["sample.size"]]

    if (cohort.iter$analysis.type == "GROWTH") {
      ###  convert long to wide
      wide_data <-
        data.table::dcast(
          data =
            if (arrow.tf) {
              getLongArrow(
                db.con = long_data,
                config = list(cohort_lookup),
                cols = long.to.wide.vars
              )
            } else {
              long_data[cohort_lookup][,
                c(getKey(long_data), long.to.wide.vars),
                with = FALSE
              ]
            },
          formula = ID ~ GRADE + CONTENT_AREA,
          sep = ".", drop = FALSE,
          value.var = c("VALID_CASE", long.to.wide.vars)
        )

      ###  Exclude kids missing IN current AND most recent year/grade's data
      ###  (keeps them if in data with missing score, but set to "VALID_CASE")
      tmp.vcase <-
        paste0(
          paste0("VALID_CASE.", c(utils::tail(prior.grades, 1), current.grade)),
          ".", rep(unique(cohort.iter$content.areas), each = 2)
        )
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
    } else {  #  END "GROWTH"  --  Begin "STATUS"
      ##  Create `wide_data` and group level summary (`group_smry`)
      status_lookup <- cohort_lookup[YEAR == current.year]
      # cohort_lookup_arw |> filter(YEAR == current.year)
      if (arrow.tf) {
        status_lookup_arw <- arrow::arrow_table(status_lookup)
        status_lookup_arw <-
          status_lookup_arw$cast(target_schema = cohort_lookup_arw$schema)
      }

      wide.fmla <-
        stats::as.formula(
          paste("ID", ifelse(is.null(group), "", paste("+", group)),
                "~ GRADE + CONTENT_AREA"
          ))

      wide_data <-
        data.table::dcast(
          data =
            if (arrow.tf) {
              getLongArrow(
                db.con = long_data,
                config = list(status_lookup),
                cols = long.to.wide.vars
              )
            } else {
              long_data[status_lookup][,
                c(getKey(long_data), long.to.wide.vars),
                with = FALSE
              ]
            },
          formula = wide.fmla,
          sep = ".", drop = FALSE,
          value.var = c(focus.variable, student.factors)
        )

      if (!is.null(group)) {
        group_lookup <- if (arrow.tf) {
          getLongArrow(
              db.con = long_data,
              config = list(status_lookup),
              cols = group
            )[, c("ID", ..group)] |>
              unique() |> data.table::setkey()
        } else {
          long_data[status_lookup, c("ID", group), with = FALSE] |>
            unique() |> data.table::setkey()
        }
        wide_data <-
          wide_data[group_lookup[!duplicated(group_lookup, by = "ID")]]
          # !duplicated(...) above takes 1st entry for kids in multiple schools
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
        priors_lookup <- cohort_lookup[YEAR != current.year]

        ##  Prior grade(s) summaries to use
        smry.eval.expression <-
          paste0("PRIOR_IMV__", focus.variable,
                 " = ", "mean(", focus.variable, ", na.rm=TRUE)")
        smry.eval.expression <-
          stats::setNames(
            smry.eval.expression,
            sub("^(.*) = .*", "\\1", smry.eval.expression)
          )

        group_smry <-
          (if (arrow.tf) {
            getLongArrow(
              db.con = long_data,
              config = list(priors_lookup),
              cols = c("CONTENT_AREA", "GRADE", focus.variable, group)
            )
          } else {
            long_data[priors_lookup][,
              c("CONTENT_AREA", "GRADE", focus.variable, group),
              with = FALSE
            ]
          })[
            !is.na(get(group)),
              lapply(smry.eval.expression, \(f) eval(parse(text = f))),
            keyby = c("GRADE", "CONTENT_AREA", group)
          # ]
          # long_data[priors_lookup][,
          #     c("CONTENT_AREA", "GRADE", focus.variable, group),
          #     with = FALSE
          # ][
          #   !is.na(get(group)),
          #     lapply(smry.eval.expression, \(f) eval(parse(text = f))),
          #   keyby = c("GRADE", "CONTENT_AREA", group)
          ] |>
            data.table::dcast(
              formula = get(group) ~ GRADE + CONTENT_AREA,
              sep = ".",
              value.var = names(smry.eval.expression)
            )

        data.table::setnames(
          group_smry,
          c(group,
            paste0(names(smry.eval.expression), ".", names(group_smry)[-1])
          )
        )

        wide_data <- merge(wide_data, group_smry, by = group, all.x = TRUE)
        data.table::setnames(wide_data, group, group.vars)

        ##  Put in cross school mean for schools with no students in prior years
        for (prior.smry in grep("PRIOR_IMV__", names(wide_data), value = TRUE)) {
          tmp.grp.mean <- mean(wide_data[, get(prior.smry)], na.rm = TRUE)
          wide_data[is.na(get(prior.smry)), eval(prior.smry) := tmp.grp.mean]
        }
      }
    } ###  END "STATUS"

    ##  Subset out scale scores and student background (factors)
    ##  done above for "STATUS"
    if (cohort.iter$analysis.type == "GROWTH") {
      subset.vars <- focus.vars <-
        grep(paste0(focus.variable, "[.]", collapse = "|"),
             names(wide_data),
             value = TRUE
        )
      if (length(group)) {
        group.vars <- grep(group, names(wide_data), value = TRUE)
        subset.vars <- c(subset.vars, group.vars)
      }
      if (length(student.factors)) {
        student.vars <-
          grep(paste0(student.factors, "[.]", collapse = "|"),
               names(wide_data),
               value = TRUE
          )
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
        yr.len <- length(cohort.iter[["status.years"]])
        unq.ca <- unique(cohort.iter[["content.areas"]])

        priors_lookup <-
          data.table::data.table(
            VALID_CASE = "VALID_CASE",
            CONTENT_AREA = rep(unq.ca, yr.len),
            YEAR = rep(cohort.iter[["status.years"]], each = length(unq.ca)),
            GRADE = rep(cohort.iter[["grade.sequences"]], each = length(unq.ca))
          ) |>
            data.table::setkey(YEAR) |> # ensure lookup table is ordered by years.
              data.table::setkey(NULL)  # NULL out key so doesn't corrupt the join in dcast.

        if (!arrow.tf) {
          data.table::setkeyv(long_data, getKey(long_data))
        }

        ##  Prior grade(s) summaries to use
        smry.eval.expression <-
          paste0("PRIOR_IMV__", focus.variable,
                 " = ", "mean(", focus.variable, ", na.rm=TRUE)")
        smry.eval.expression <-
          stats::setNames(
            smry.eval.expression,
            sub("^(.*) = .*", "\\1", smry.eval.expression)
          )

        group_smry <-
          (if (arrow.tf) {
            getLongArrow(
              db.con = long_data,
              config = list(priors_lookup),
              cols = c("CONTENT_AREA", "GRADE", focus.variable, group)
            )
          } else {
            long_data[priors_lookup][,
              c("CONTENT_AREA", "GRADE", focus.variable, group),
              with = FALSE
            ]
          })[
          # long_data[priors_lookup][,
          #   c("CONTENT_AREA", "GRADE", focus.variable, group),
          #   with = FALSE
          # ][
            !is.na(get(group)),
              lapply(smry.eval.expression, \(f) eval(parse(text = f))),
            keyby = c("GRADE", "CONTENT_AREA", group)
          ] |>
            data.table::dcast(
              formula = get(group) ~ GRADE + CONTENT_AREA,
              sep = ".",
              value.var = names(smry.eval.expression)
            )

        ##  Merge in summaries for appropriate grade/subject
        for (prior.smry in names(group_smry)[-1]) {
          smry_names <-
            c(group, paste0("PRIOR_IMV__", focus.variable)) |>
              paste(prior.smry, sep = ".")
          if (!smry_names[1] %in% names(wide_data)) next
          sub_grp_smry <- group_smry[, c("group", prior.smry), with = FALSE]
          setnames(sub_grp_smry, smry_names)

          wide_data <-
            merge(wide_data, sub_grp_smry, by = smry_names[1], all.x = TRUE)

          # Use cross-school mean for schools with no students in prior years
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

      ##  Clean up all NA groups (if necessary)
      if (length(group) &
          toupper(impute.method) %in% c("PANIMPUTE", "JOMOIMPUTE")
      ) {
        for (gv in group.vars) {
          if (cohort.iter$analysis.type == "GROWTH") {
            tmp.imp.var <- paste0(focus.variable, ".", sub(".+?[.]", "", gv))
            tmp.fixd.ef <- focus.vars %w/o% tmp.imp.var

            ## Check for clusters with all NA values
            na_check <-
              wide_data[,
                list(ALL_NA = all(is.na(get(tmp.imp.var)))),
                keyby = gv
              ]

            if (length(all.na.ids <- na_check[ALL_NA == TRUE, get(gv)])) {
              wide_data[get(gv) %in% all.na.ids, eval(gv) := NA]
            }
          } else {
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
                wide_data[get(gv) %in% all.na.ids, eval(gv) := NA]
              }
            }
          }
        }
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
      }

      ##  arguments/objects for `mice`
      tmp.meth <-
        c(rep("sample", length(group.vars)),
          rep(impute.method, length(focus.vars))
         ) |>
          stats::setNames(c(group.vars, focus.vars))
      if (!is.null(focus.years)) {
        hifocus.vars <-
          paste(
            focus.variable,
            cohort_lookup[YEAR %in% focus.years][, GRADE],
            cohort_lookup[YEAR %in% focus.years][, CONTENT_AREA],
            sep = "."
          ) |> unique()
        tmp.meth[focus.vars %w/o% hifocus.vars] <- "pmm"
      } else {
        hifocus.vars <- focus.vars
      }

      tmp.pred <- names(wide_data) %w/o% "ID"

      ###   `mice`/... argument alignment
      if (!"blocks" %in% names(call)) {
        tmp.blok <-
          mice::make.blocks(
            data = c(group.vars, focus.vars),
            calltype = "formula"
          )
      } else {
        tmp.blok <- call[["blocks"]]
      }

      if (!"formulas" %in% names(call)) {
        tmp.fmla <- list()
        if (length(group)) {
          for (gv in group.vars) {
            if (!any(is.na(wide_data[[gv]]))) {
              tmp.meth[[gv]] <- ""
              next
            }
            tmp.fmla[[gv]] <-
              paste(gv, "~ 1") |>
                stats::as.formula()
          }
        }
        for (fv in focus.vars) {
          tmp.iv <- focus.vars %w/o% fv
          tmp.grd.subj <- sub(".+?[.]", "", fv)

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

            rhs <-
              ifelse(!formula.specs$simple,
                paste0(tmp.iv,
                       "*I(clusterMeans(", tmp.iv, ", ", tmp.gv, "))",
                       collapse = " + "
                ),
                paste0(tmp.iv, collapse = " + ")
              )

            if (length(tmp.pr <- grep("PRIOR_IMV__", tmp.pred))) {
              rhs <-
                c(rhs,
                  paste0(tmp.pred[tmp.pr], collapse = " + ")
                )
            }

            if (length(tmp.sf)) {
              rhs <- c(rhs, paste0(tmp.sf, collapse = " + "))
            }

            if (formula.specs$random.slope) {
              rslope <- paste(" +", tmp.iv, collapse = "")
            } else {
              rslope <- NULL
            }

            if (toupper(impute.method) %in% c("PANIMPUTE", "JOMOIMPUTE")) {
              rhs <- c(rhs, paste0("(1", rslope, "|", tmp.gv, ")"))
            }

            tmp.fmla[[fv]] <-
              paste(fv, "~", paste(rhs, collapse = " + ")) |>
                stats::as.formula()
          } else {
            tmp.fmla[[fv]] <-
              paste(fv, "~", paste(c(tmp.iv, tmp.sf), collapse = " + ")) |>
                stats::as.formula()
          }
        }
      } else {
        tmp.fmla <- call[["formulas"]]
      }

    if (parexecute == "SEQ") {
      imp_data <-
        suppressWarnings(
          mice::mice(
            data = wide_data,
            m = M, method = tmp.meth,
            visitSequence = c(group.vars, focus.vars),
            blocks = tmp.blok, formulas = tmp.fmla,
            maxit = maxit, seed = seed, print = verbose,
            ...
            # where = tmp.where, ignore = tmp.ign,
            # can't use `where` or `ignore` with panImpute
        ))
    } else {
      invisible(gc())
      imp_data <-
        suppressWarnings(
          parMICE(
            data = wide_data, m = M, maxit = maxit, meth = tmp.meth,
            bloks = tmp.blok, frmlas = tmp.fmla, visit = c(group.vars, focus.vars),
            seed = seed, par.config = parallel.config, sample.n = sample.size
        ))
    }

    ##  Save some diagnostic plots
    if (is.null(sample.size)) {
      if (cohort.iter$analysis.type == "GROWTH") {
        imp.type <- "GROWTH_"
      } else {
        imp.type <- "STATUS_"
      }

      grDevices::pdf(
        file =
          file.path(
            diagn.dir,
            paste0("Grade_", current.grade, "_", current.year, "_", imp.type,
                  gsub("[.]", "", impute.method),
                  "_M_", M, "__maxit_", maxit,
                  ifelse(is.null(group), "", paste0("_x_", group)),
                  "__converge.pdf"
            )
          )
      )
      print(graphics::plot(imp_data))
      invisible(grDevices::dev.off())

      grDevices::pdf(
        file =
          file.path(
            diagn.dir,
            paste0("Grade_", current.grade, "_", current.year, "_", imp.type,
                  gsub("[.]", "", impute.method),
                  "_M_", M, "__maxit_", maxit,
                  ifelse(is.null(group), "", paste0("_x_", group)),
                  "__density.pdf"
            )
          ),
        width = 11, height = 8
      )
      tryCatch(
        print(mice::densityplot(
          imp_data,
          eval(parse(text = paste0("~", paste(hifocus.vars, collapse = " + "))))
        )),
        error = function(e) TRUE
      ) -> err.tf
      if (is.logical(err.tf)) print(mice::densityplot(imp_data)) else rm(err.tf)
      invisible(grDevices::dev.off())
    }

    ###  Format and store results
    if (cohort.iter$analysis.type == "STATUS") {
        long.ids <- c("ID", ".imp", gv)
    } else {
        long.ids <- c("ID", ".imp")
    }

    if (is.null(sample.size) || parexecute == "SEQ") {
      imp_data <-
        data.table::as.data.table(
          mice::complete(imp_data, action = "long", include = TRUE)
        ) |>
          data.table::melt(
            id.vars = long.ids,
            value.name = focus.variable,
            variable.name = "GRADE",
            variable.factor = FALSE,
            measure.vars = focus.vars
          )
    } else {
        imp_data <-
          data.table::melt(
            data = imp_data,
            id.vars = long.ids,
            value.name = focus.variable,
            variable.name = "GRADE",
            variable.factor = FALSE,
            measure.vars = focus.vars
          )
    }

    if (cohort.iter$analysis.type == "STATUS") {
        imp_data[,
          CONTENT_AREA := gsub(paste0(".*", current.grade, "[.]"), "", GRADE)
        ][,
          GRADE := current.grade
        ][,
          YEAR := current.year
        ]
        if (!is.null(group)) {
          no_group <-
            imp_data[is.na(get(gv)), .(ID, GRADE, CONTENT_AREA)][,
              NA_GrpVar_IMP := TRUE
            ]
          data.table::setkey(no_group, ID, GRADE, CONTENT_AREA)
          data.table::setkey(imp_data, ID, GRADE, CONTENT_AREA)
          imp_data <- no_group[imp_data]
          imp_data[
            NA_GrpVar_IMP == TRUE, eval(gv) := NA
          ][,
            NA_GrpVar_IMP := NULL
          ]
        }
    } else {
        imp_data[, GRADE := gsub(paste0(focus.variable, "."), "", GRADE)]
        unq.ca <- unique(cohort.iter[["content.areas"]])

        imp_data[,
          CONTENT_AREA :=
            gsub(
              paste(
                paste0(".*", c(prior.grades, current.grade), "."),
                collapse = "|"
              ),
              "", GRADE
            )
        ][,
          GRADE :=
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

    imp_data <-
        data.table::dcast(
          data = imp_data,
          formula = as.formula(
            paste(
              paste(long.ids %w/o% ".imp", collapse = " + "),
                    "+ CONTENT_AREA + GRADE + YEAR ~ .imp"
            )),
          value.var = focus.variable
        )

    if (cohort.iter$analysis.type == "STATUS") {
      rekey <- key(imp_data) %w/o% gv
      imp_data[, eval(gv) := NULL]
      setkeyv(imp_data, rekey)

      if (arrow.tf) {
        long_final <-
          dplyr::semi_join(long_data, status_lookup_arw) |>
            data.table::as.data.table()
        long_final[,
            CONTENT_AREA := as.character(CONTENT_AREA)
        ][, YEAR := as.character(YEAR)
        ][, GRADE := as.character(GRADE)
        ]
        data.table::setcolorder(
          long_final,
          names(long_data)[names(long_data) %in% names(long_final)]
        )
      } else {
        long_final <-
          long_data[status_lookup][,
            unique(c(key(long_data), long.to.wide.vars)), with = FALSE]
        data.table::setcolorder(
          long_final,
          names(long_data)[names(long_data) %in% names(long_final)]
        )
      }
    } else {
      meas.list <- vector(mode = "list", length = length(long.to.wide.vars))
      meas.list <-
        lapply(
          long.to.wide.vars,
          \(f) {meas.list[[f]] <- grep(paste0("^", f, "[.]"), names(wide_data))}
        )
      names(meas.list) <- long.to.wide.vars

      long_final <-
        data.table::melt(
          data = wide_data,
          id.vars = "ID",
          variable.name = "GRADE",
          measure.vars = meas.list
        )[,
          VALID_CASE := "VALID_CASE"
      ][, CONTENT_AREA :=
            as.character(factor(GRADE, labels = cohort_lookup[["CONTENT_AREA"]]))
      ][, YEAR :=
            as.character(factor(GRADE, labels = cohort_lookup[["YEAR"]]))
      ][, GRADE :=
            as.character(factor(GRADE, labels = cohort_lookup[["GRADE"]]))
      ]
      data.table::setcolorder(
        long_final,
        names(long_data)[names(long_data) %in% names(long_final)]
      )
    }

    data.table::setnames(
      imp_data,
      as.character(0:M),
      c(focus.variable, paste0(focus.variable, ".IMP_", 1:M))
    )
    tmp.key <- c(data.table::key(imp_data), focus.variable)# "ID", ..., "YEAR"
    data.table::setkeyv(imp_data, tmp.key)
    data.table::setkeyv(long_final, tmp.key)
    if (nrow(imp_data) >= nrow(long_final)) {
      arrow::write_dataset(
          dataset =
            long_final[imp_data][,
                VALID_CASE := "VALID_CASE"
            ][, TEMP_Imp_Cohort_N := K],
          path = tempdir(),
          format = "parquet",
          existing_data_behavior = "overwrite",
          partitioning = "TEMP_Imp_Cohort_N"
      )
    } else {
      arrow::write_dataset(
          dataset = imp_data[long_final][, TEMP_Imp_Cohort_N := K],
          path = tempdir(),
          format = "parquet",
          existing_data_behavior = "overwrite",
          partitioning = "TEMP_Imp_Cohort_N"
      )
    }
    rm(list = c("long_final", "wide_data", "imp_data"))
    invisible(gc())
    message(
      "\n\t\tFinished with ", current.year, " Grade ", current.grade,
      ": ", date()
    )
  }  ###  END K

  ###   Compile and format final imputed data set
  final_imp_data <-
    arrow::open_dataset(
        sources = tempdir(),
        partitioning = "TEMP_Imp_Cohort_N"
    ) |> data.table::as.data.table()

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
    if (arrow.tf) {
      long_data <- data.table::as.data.table(long_data)
      long_data[,
          CONTENT_AREA := as.character(CONTENT_AREA)
      ][, YEAR := as.character(YEAR)
      ][, GRADE := as.character(GRADE)
      ]
    }

    fin.key <- c(getKey(final_imp_data), focus.variable)

    ##  Remove group, student.factors and other extraneous variables
    rm.nms <-
      c("TEMP_Imp_Cohort_N",
        names(final_imp_data)[
          names(final_imp_data) %in% names(long_data)
        ] %w/o% fin.key
      )

    final_imp_data[, (rm.nms) := NULL]

    data.table::setkeyv(long_data, fin.key)
    data.table::setkeyv(final_imp_data, fin.key)

    final_imp_data <- final_imp_data[long_data]
  }

  message(
    "\n\tImputation with `imputeLongData` completed in ",
    convertTime(timeTaken(started.impute))
  )
  return(final_imp_data)
}
