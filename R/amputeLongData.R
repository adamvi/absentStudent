#' @title amputeLongData
#' @description Function removes cases to produce desired patterns and proportions of missing data.
#' @param long_data The complete dataset in which to create missing values.
#' @param plot.dir Directory path for diagnostic plots to be placed. Default is
#'  the working directory. Additional subdirectories are created internally as needed.
#'  Default: 'Missing_Data'
#' @param growth.config An elongated `SGP` style config script with an entry for
#'  each grade/content_area/year cohort definition. Configs will be used to subset
#'  the `long_data` provided as required for cohort level data amputation.
#'  Default: NULL
#' @param status.config An elongated `SGP` style config script with an entry for
#'  each grade/content_area/year cohort definition. Configs will be used to subset
#'  the `long_data` provided as required for cohort level data amputation. Unlike
#'  a `growth.config` entry, `status.config` entries use data from the same grade,
#'  but from the prior year(s) (i.e. not individual variables). For example you
#'  might ampute 3rd grade ELA scores based on a previous year's 3rd grade
#'  school mean scale score, FRL status, etc. Default: NULL
#' @param focus.variable The variable to be amputed. Currently only a single
#'  variable is allowed. Default: 'SCALE_SCORE'
#' @param student.factors Demographic, (prior) academic achievement or other
#'  student level background information that will be used in the construction
#'  of the weighted scores that define the probability of being missing. Any
#'  i For example, if `c('SCALE_SCORE', 'FREE_REDUCED_LUNCH_STATUS')` is provided,
#'  student level scores and FRL status will be used to calculate weighted scores.
#'  The default (NULL) means that no factors are considered, creating a 'missing
#'  completely at random' (MCAR) missing data pattern. Default: NULL
#' @param group Grouping indicator (e.g. institution ID) used to construct
#'  group level means of any `student.factors` provided. Default: NULL
#' @param additional.data The function will return only data that is required for
#'  a SGP analysis (based on growth.config). This allows for the addition of more
#'  data (e.g. prior year(s) data not used for data amputation). Default: NULL
#' @param compact.results By default (`FALSE`), the function will return a list of
#'  longitudinal datasets with the current (amputed) and prior (unchanged) student
#'  records. This is helpful for diagnostics and ease of use, but also produces more
#'  redundant prior data than needed. Setting this argument to TRUE returns a single
#'  data.table object with a `TRUE/FALSE` indicator column added for each requested
#'  amputation. This flag can be used to make the `SCALE_SCORE` (and/or other variables)
#'  NA in subsequent use cases. Default: `FALSE`
#' @param ampute.var.weights Relative weights assigned to the `student.factors`. Default
#'  is NULL meaning the weighted sum scores will be calculated with equal weight (=1).
#'  A named list can be provided with the desired relative weights. For example,
#'  `list(SCALE_SCORE=3, FREE_REDUCED_LUNCH_STATUS=2, SCHOOL_NUMBER=1)` will weight
#'  a student's scale scores by a factor of 3 and FRL by 2, with all (any other `student.factors`)
#'  remaining at the default of 1. This includes school level aggregates in this
#'  example. Note that differential weights for institutions should be placed at
#'  the end of the list. If institution IDs (e.g., `SCHOOL_NUMBER`) are omitted from
#'  the list, the aggregates will be given the same weight as the associated student
#'  level variable. In the given example, `SCALE_SCORE` and school mean scale score
#'  would be given a relative weight of 3. This argument is ignored when
#'  `student.factors = NULL`. Default: NULL
#' @param reverse.weight The current default for `ampute.args$type` is `'RIGHT'`, which
#'  means that students with high weighted scores have the highest probability for
#'  amputation. This makes sense for high % FRL schools, but not for high achieving
#'  students and/or students in high achieving schools. This function inverses the
#'  variable(s) individual (and institutional mean) value(s) so that higher weight
#'  is given to lower scores/means. This argument is ignored when `student.factors = NULL`.
#'  Default: 'SCALE_SCORE'
#' @param ampute.args Variables to be used in the mice:ampute.continuous function.
#'  Currently only `prop` and `type` can be modified. See ?mice::ampute.continuous
#'  and ?mice::ampute for more information. The `prop` piece is inexact and has required
#'  some modification. Its still imprecise, particularly for values away
#'  from 0.5 (50% missing). Also, the max missingness is 85%, and for that you need
#'  to set prop=0.95 or greater. Note that the prop gives a total proportion missing,
#'  which accounts for missingness already included in the data due to students with
#'  incomplete longitudinal data histories. For example, if a cohort starts with 5%
#'  students missing due to incomplete histories, an additional 25% will be made
#'  missing to get to the 30% (default) missingness. This last point has been dealt
#'  with in some regards with the next argument, which removes these cases first.
#'  Default: list(prop = 0.3, type = 'RIGHT')
#' @param complete.cases.only Should cases without the most recent prior and current
#'  score be removed? This removes students with partial longitudinal histories
#'  from the most recent prior (e.g., 2019) to the current year (e.g., 2021),
#'  producing a 'complete' dataset that is easier to interpret. Default: TRUE
#' @param partial.fill Should an attempt be made to fill in some of the demographic
#'  and institution ID information based on students previous values? Part of the
#'  process of the amputeScaleScore function is to take the longitudinal data (`long_data`)
#'  and then first widen and then re-lengthen the data, which creates holes for
#'  students with incomplete longitudinal records. This part of the function fills
#'  in these holes for students with existing missing data. Default: TRUE
#' @param invalidate.repeater.dups Students who repeat a grade will get missing data
#'  rows inserted for the grade that they 'should' be in for any given year. This leads
#'  to duplicated cases that can lead to problems in the SGP analyses. This argument
#'  returns those cases with the `VALID_CASE` variable set to `'INVALID_CASE'`. Default: TRUE
#' @param seed A random seed set for the amputation process to allow for replication
#'  of results, or for alternative results using the same code. Default: 4224
#' @param R The number of amputed datasets to return. Default: 10
#' @return Function returns either a list (default) of amputed longitudinal data
#'  sets or a single data set containing additional columns of indicators for
#'  records to be removed.
#' @details From the amputation process specified, the function returns either a list
#'  (default) of `R` amputed datasets or a single data set with `R` columns of missing
#'  record indicators. The datasets will exclude data for students not used in any of
#'  the specified growth.config or status.config cohorts, unless the additional.data 
#'  argument has been included.
#' @examples 
#'  \dontrun{
#'    data_to_ampute <- SGPdata::sgpData_LONG_COVID
#'
#'    ###   Read in STEP 0 SGP configuration scripts
#'    source("SGP_CONFIG/STEP_0/Ampute_2021/Growth.R")
#'    source("SGP_CONFIG/STEP_0/Ampute_2021/Status.R")
#'
#'    Test_Data_LONG <-
#'        amputeScaleScore(
#'            long_data = data_to_ampute,
#'            growth.config = growth_config_2021,
#'            status.config = status_config_2021,
#'            R = 1)
#'  }
#' @seealso 
#'  \code{\link[mice]{ampute}}, \code{\link[mice]{ampute.continuous}}
#' @rdname amputeLongData
#' @export 
#' @author Adam R. VanIwaarden \email{avaniwaarden@nciea.org}
#' @keywords misc
#' @keywords models
#' @importFrom data.table SJ setnames setkey setkeyv dcast
#' @importFrom mice ampute ampute.continuous

amputeLongData =
  function(
    long_data,
    plot.dir = "Missing_Data",
    growth.config = NULL,
    status.config = NULL,
    focus.variable = "SCALE_SCORE",
    student.factors = NULL,
    group = NULL,
    additional.data = NULL,
    compact.results = FALSE,
    ampute.var.weights = NULL,
    reverse.weight = "SCALE_SCORE",
    ampute.args = list(prop = 0.3, type = "RIGHT"),
    complete.cases.only = TRUE,
    partial.fill = TRUE,
    invalidate.repeater.dups = TRUE,
    seed = 4224L,
    R = 10
  ) {

   ###   Avoid "Undefined global functions or variables:" from R CMD check
    GRADE <- ID <- VALID_CASE <- YEAR <- SCALE_SCORE <-
    ACHIEVEMENT_LEVEL <- TMP_MCAR_PROB <- NULL

    ###  Initial checks
    if (!dir.exists(
      diagn.dir <- file.path(plot.dir, "amputeLongData", "diagnostic_plots")
    )) {
      dir.create(diagn.dir, recursive = TRUE)
    }

    ###   Combine and augment config lists
    configs <- growth.config

    if (!is.null(growth.config)) {
      for (f in seq(configs)) configs[[f]][["analysis.type"]] <- "GROWTH"
    }
    if (!is.null(status.config)) {
      for (f in seq(status.config)) status.config[[f]][["analysis.type"]] <- "STATUS"
      configs <- c(configs, status.config)
    }

    if (is.null(configs)) stop("Either a 'growth.config' or 'status.config' must be supplied")

    long.to.wide.vars <- c(focus.variable, group, student.factors) # c("CONTENT_AREA", "GRADE",


    ###  Cycle through configs to get results by cohort
    res.list <- vector(mode = "list", length = R)

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
          c(getKey(long_data), long.to.wide.vars),
          with = FALSE
        ]

      # prior.years <- head(unique(cohort.iter$panel.years), -1)
      prior.year <- utils::tail(unique(cohort.iter$panel.years), 2)[1]
      current.year <- utils::tail(unique(cohort.iter$panel.years), 1)
      current.grade <- utils::tail(unique(cohort.iter$grade.sequences), 1)

      if (cohort.iter$analysis.type == "GROWTH") {
        ###   convert long to wide
        wide_data <-
          data.table::dcast(
            data = tmp_long,
            formula = ID ~ YEAR + CONTENT_AREA,
            sep = ".", drop = FALSE,
            value.var = c("GRADE", long.to.wide.vars)
          )
      } else {  #  END "GROWTH"  --  Begin "STATUS"
        ###   Create institution level summaries
        current.year <- utils::tail(cohort.iter$panel.years, 1)

        wide_data <-
          tmp_long[YEAR %in% current.year, c("ID", student.factors), with = FALSE] # assuming not using any other current year data.

        if (!is.null(reverse.weight)) {
          for (rev.var in reverse.weight) {
            # tmp_long.priors[, eval(rev.var) := -1*get(rev.var)]
            wide_data[, eval(rev.var) := -1 * get(rev.var)]
          }
        }

        if (length(student.factors)) {
          for (demog in student.factors) {
            wide_data[, eval(demog) := as.integer(factor(get(demog)))-1L]
          }
        }

        if (length(group)) {
          for (wav in student.factors) {
            wide_data[,
              paste0("TMP_IMV__", wav, "_", group) := mean(get(wav), na.rm = TRUE),
              by = list(get(group))
            ]
          }

          ##    Put in cross group mean for group with no students in prior years
          for (tmp.inst.smry in grep("TMP_IMV__", names(wide_data), value = TRUE)) {
            tmp.mean <- mean(wide_data[, get(tmp.inst.smry)], na.rm = TRUE)
            wide_data[is.na(get(tmp.inst.smry)), eval(tmp.inst.smry) := tmp.mean]
          }
          wide_data[, eval(group) := NULL]
        }

        #  remove columns that are all NA (e.g., SGP for 3rd grade priors)
        wide_data <- wide_data[,
          names(wide_data)[!unlist(lapply(names(wide_data), \(f) all(is.na(wide_data[,get(f)]))))], with = FALSE]
        wide_data <- stats::na.omit(wide_data)

        ###   Create long_final with only the "current" year (last elements of the config)
        ###   More thorough to do it with config than just long_final <- tmp_long[YEAR %in% current.year]
        cohort.lookup <- SJ("VALID_CASE", utils::tail(cohort.iter[["content.areas"]], 1),
          utils::tail(cohort.iter[["panel.years"]], 1), utils::tail(cohort.iter[["grade.sequences"]], 1))
        setkeyv(long_data, getKey(long_data))

        long_final <- long_data[cohort.lookup][, unique(c("VALID_CASE", "ID", "YEAR", long.to.wide.vars)), with = FALSE]
      } ###  END "STATUS"

      ##    Subset out scale scores and student.factors
      if (cohort.iter$analysis.type == "GROWTH") {  #  done above for "STATUS"
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
        }

        #  create subset of wide data with only variables to be used in imputation
        wide_data <- wide_data[,
                            grep(paste0("^ID$|", paste(subset.vars, collapse = "|")),
                                names(wide_data)), with = FALSE]

        ##    Create institutional level averages of achievement and student.factors
        ##    TMP_IMV__  -  temp institutional mean variable

# Just for most recent prior year (of current group)?

        if (length(group)) {
          tmp.inst.var <- paste0(group, ".", prior.year)
          for (wav in subset.vars) {
            wide_data[, paste0("TMP_IMV__", wav, "_", group) :=
                              mean(get(wav), na.rm = TRUE),
                            by = list(get(tmp.inst.var))]
          }
        }

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
      }  ###  END "GROWTH"

      ###   AMPUTE
      ###   https://rianneschouten.github.io/mice_ampute/vignette/ampute.html

      test.result <- mice::ampute(wide_data)

      ##    Reduce amputation analysis data to non-NA data and standardize
      wide_data.std <- scale(wide_data[,-1,])

      tmp.weights <- matrix(rep(1, ncol(wide_data.std)), nrow = 1)
      if (!is.null(ampute.var.weights)) {
        for (n in names(ampute.var.weights)) {
          tmp.weights[grep(n, names(wide_data)[-1])] <- ampute.var.weights[[n]]
        }
      }
      #  moderate the STATUS weight for (now current) achievement (0.65 is approx low correlation between current & prior)
      if (cohort.iter$analysis.type == "STATUS") {
        if (length(ss.indx <- grep("^SCALE_SCORE$", names(wide_data)[-1])) > 0) {
          tmp.weights[ss.indx] <- tmp.weights[ss.indx]*0.65
        }
      }

      tmp.scores <- apply(wide_data.std, 1, function (x) tmp.weights %*% x)
      if (is.null(ampute.args$type)) ampute.args$type <- "RIGHT"

      ##    mice::ampute.continuous doesn't do a good job at getting the right `prop`
      ##    value down  :(  Do a "burn in" to get it closer

      if (!too.low.tf) {
        adj.prop <- target.prop
        fin.props <- c()
        res.props <- c()
        shrink.tol <- 1L
        for (j in 1:25) {
          adj_mask <- mice::ampute.continuous(
                                P = rep(2, nrow(wide_data)), prop = adj.prop,
                                scores = list(tmp.scores), type = ampute.args$type)
          res.prop <- (sum(adj_mask[[1]] == 0)/nrow(wide_data))
          if (!inrange(res.prop, ltol[shrink.tol], utol[shrink.tol])) {
            # constrain adjustment to keep from getting too big/small too fast
            tmp.ratio <- target.prop/res.prop
            if (tmp.ratio > 1.15) tmp.ratio <- 1.15
            if (tmp.ratio < 0.85) tmp.ratio <- 0.85
            adj.prop <- adj.prop * tmp.ratio
            if (adj.prop > 0.999) adj.prop <- 0.999
            if (adj.prop < 0.001) adj.prop <- 0.001
          } else {
            fin.props <- c(fin.props, adj.prop)
            res.props <- c(res.props, res.prop)
            shrink.tol <- shrink.tol + 1
          }
        }
        if (!is.null(fin.props)) {
          fin.prop <- stats::weighted.mean(x = fin.props, w = (1/abs(1-target.prop/res.props)))
          if (is.na(fin.prop)) fin.prop <- mean(fin.props)
        } else fin.prop <- target.prop

        if (compact.results) {
          amp.tf <- data.table(VALID_CASE = "VALID_CASE", ID = wide_data$ID, CONTENT_AREA = utils::tail(cohort.iter[["content.areas"]], 1),
                    GRADE = utils::tail(cohort.iter[["grade.sequences"]], 1), YEAR = utils::tail(cohort.iter[["panel.years"]], 1))
          setkey(amp.tf)
        }

        ##    Get R sets of amputation candidates
        for (amp.m in seq(R)) {
          set.seed(seed*amp.m)
          mask_var <- mice::ampute.continuous(
                                  P = rep(2, nrow(wide_data)), prop = fin.prop,
                                  scores = list(tmp.scores), type = ampute.args$type)

          if (compact.results) {
            amp.tf[, eval(paste0("AMP_", amp.m)) := mask_var[[1]] == 0L]
          } else {
            amp.ids <- wide_data[which(mask_var[[1]] == 0L), ID]

            res.list[[amp.m]][[K]] <- copy(long_final)
            res.list[[amp.m]][[K]][YEAR == current.year & ID %in% amp.ids, SCALE_SCORE := NA]
            if ("ACHIEVEMENT_LEVEL" %in% names(res.list[[amp.m]][[K]])) {
              res.list[[amp.m]][[K]][YEAR == current.year & ID %in% amp.ids, ACHIEVEMENT_LEVEL := NA]
            }
          }
        }
        if (compact.results) {
          setkeyv(long_final, key(amp.tf))
          res.list[[K]] <- amp.tf[long_final]
        }
      } else {
        max.scores <- utils::head(rev(sort(tmp.scores)), pick.miss*R)
        id.list <- wide_data[which(tmp.scores %in% max.scores), ID]
        score.prob <- tmp.scores[which(tmp.scores %in% max.scores)]
        for (amp.m in seq(R)) {
          set.seed(seed*amp.m)
          amp.ids <- sample(x = id.list, size = pick.miss, prob = score.prob)
          res.list[[amp.m]][[K]] <- copy(long_final)
          res.list[[amp.m]][[K]][YEAR == current.year & ID %in% amp.ids, SCALE_SCORE := NA]
          if ("ACHIEVEMENT_LEVEL" %in% names(res.list[[amp.m]][[K]])) {
            res.list[[amp.m]][[K]][YEAR == current.year & ID %in% amp.ids, ACHIEVEMENT_LEVEL := NA]
          }
        }
      }
    }

    if (compact.results) {
      final.amp.list <- rbindlist(res.list, use.names = TRUE)
      if (remove.tmp.amp.var) invisible(final.amp.list[, TMP_MCAR_PROB := NULL])
      if (!is.null(additional.data)) {
        final.amp.list <- rbindlist(
          list(final.amp.list, additional.data[, names(final.amp.list) %w/o% paste0("AMP_", seq(R)), with = FALSE]), fill = TRUE)
      }
      if (invalidate.repeater.dups) {
        setkeyv(final.amp.list, getKey(final.amp.list))
        setkeyv(final.amp.list, key(final.amp.list) %w/o% "GRADE")
        dup.ids <- final.amp.list[which(duplicated(final.amp.list, by = key(final.amp.list))), ID]
        final.amp.list[ID %in% dup.ids & is.na(SCALE_SCORE), VALID_CASE := "INVALID_CASE"]
      }
    } else {
      final.amp.list <- vector(mode = "list", length = R)

      for (L in seq(R)) {
        final.amp.list[[L]] <- rbindlist(res.list[[L]], fill = TRUE)

        if (remove.tmp.amp.var) invisible(final.amp.list[[L]][, TMP_MCAR_PROB := NULL])

        if (!is.null(additional.data)) {
          final.amp.list[[L]] <- rbindlist(
            list(final.amp.list[[L]], additional.data[, names(final.amp.list[[L]]), with = FALSE]), fill = TRUE)
        }
        if (invalidate.repeater.dups) {
          setkeyv(final.amp.list[[L]], getKey(final.amp.list[[L]]))
          setkeyv(final.amp.list[[L]], key(final.amp.list[[L]]) %w/o% "GRADE")
          dup.ids <- final.amp.list[[L]][which(duplicated(final.amp.list[[L]], by = key(final.amp.list[[L]]))), ID]
          final.amp.list[[L]][ID %in% dup.ids & is.na(SCALE_SCORE), VALID_CASE := "INVALID_CASE"]
        }
      }
    }
    return(final.amp.list)
}
