#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param long_data PARAM_DESCRIPTION
#' @param shift.group PARAM_DESCRIPTION, Default: c("ID", "CONTENT_AREA")
#' @param shift.period PARAM_DESCRIPTION, Default: 'YEAR'
#' @param shift.max PARAM_DESCRIPTION, Default: 2
#' @param invariant.vars PARAM_DESCRIPTION, Default: NULL
#' @param conditional.vars PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[data.table]{J}}, \code{\link[data.table]{setkey}}, \code{\link[data.table]{shift}}, \code{\link[data.table]{rbindlist}}
#' @rdname completeData
#' @export 
#' @author AUTHOR [AUTHOR_2]
#' @keywords KEYWORD_TERM
#' @author AUTHOR [AUTHOR_2]
#' @keywords KEYWORD_TERM
#' @importFrom data.table CJ setkeyv shift rbindlist
completeData =
  function(
    long_data,
    shift.group = c("ID", "CONTENT_AREA"),
    shift.period = "YEAR",
    shift.max = 2,
    invariant.vars = NULL,
    conditional.vars = NULL
  ) {
    # GRADE_LAG_1 <- GRADE_LAG_2 <- GRADE_LEAD_1 <- GRADE_LEAD_2 <-
    VALID_CASE <- GRADE <- TMP_GRADE <- NULL
    
    # Set data key as necessary
    shift.key <- c("ID", "CONTENT_AREA", "YEAR", "GRADE", "VALID_CASE")
    setkeyv(long_data, shift.key[shift.key %in% names(long_data)])

    ### Utility functions
    `completeDT` =
      function(DT, cols, defs = NULL) {
        mDT = do.call(data.table::CJ, c(DT[, ..cols], list(unique = TRUE)))
        res = DT[mDT, on = names(mDT)]
        if (length(defs)) {
          res[,
            names(defs) := Map(replace, .SD, lapply(.SD, is.na), defs),
            .SDcols = names(defs)
          ]
        }
        res[]
      }

    if ("VALID_CASE" %in% names(long_data)) {
      valid.case <- "VALID_CASE"
    } else {
      valid.case <- NULL
    }
    grade.levels <- unique(long_data[["GRADE"]]) %w/o% NA

    long_data_complete <- completeDT(long_data, cols = c(shift.period, shift.group))
    if (!is.null(valid.case)) {
      long_data_complete[, COMPLETE_CASE_ADDED := is.na(VALID_CASE)]
      long_data_complete[is.na(VALID_CASE), VALID_CASE := "VALID_CASE"]
    }

    ###   Fill in within-year (shift.period) missing data
    data.table::setkeyv(long_data_complete, c("ID", shift.period, valid.case))

    shift.vars <- c("GRADE", invariant.vars)
    shift.var.names <-
      paste0(c(shift.vars, names(conditional.vars)), "_LAG_", 0)
    # takes care of 2 content area case. might not for 3...
    # shift.var.names <- paste(rep(shift.vars, each = 2),
    #                          c("LAG", "LEAD"), 0, sep = "_")
    long_data_complete[,
      (shift.var.names) := data.table::shift(.SD, type = "cyclic"), # n = c(1, -1),
      .SDcols = c(shift.vars, names(conditional.vars)),
      by = c("ID", shift.period)
    ]
    for (svr in c(shift.vars, names(conditional.vars))) {
      long_data_complete[is.na(get(svr)),
        eval(svr) := get(paste0(svr, "_LAG_0"))
      ]
    }
    long_data_complete[, (shift.var.names) := NULL]

    ###   Fill in with proximate year(s)
    ##    GRADE values
    data.table::setkeyv(long_data_complete, c(shift.group, shift.period, valid.case))
    shift.var.names <-
      paste("GRADE",
            rep(c("LAG", "LEAD"), each = shift.max),
            rep(1:shift.max, 2),
        sep = "_"
      )
    long_data_complete[,
      (shift.var.names) := data.table::shift(.SD, n = c(1:shift.max, -(1:shift.max))),
      .SDcols = "GRADE",
      by = shift.group
    ]

    long_data_complete[, TMP_GRADE := as.numeric(GRADE)]
    for (smx in 1:shift.max) {
      long_data_complete[is.na(TMP_GRADE),
        TMP_GRADE := as.numeric(get(paste0("GRADE_LAG_", smx))) + smx
      ]
      long_data_complete[is.na(TMP_GRADE),
        TMP_GRADE := as.numeric(get(paste0("GRADE_LEAD_", smx))) - smx
      ]
    }
    mode(long_data_complete$TMP_GRADE) <- mode(long_data_complete$GRADE)
    long_data_complete[is.na(GRADE), GRADE := TMP_GRADE]
    long_data_complete[, TMP_GRADE := NULL]
    long_data_complete[, (shift.var.names) := NULL]
    long_data_complete <- long_data_complete[GRADE %in% grade.levels, ]

    ##    Time-invariant variable values
    if (length(invariant.vars)) {
      data.table::setkeyv(
        long_data_complete,
        c(shift.group, shift.period, valid.case)
      )
      shift.var.names <-
        paste(rep(invariant.vars, each = 2*shift.max),
              rep(c("LAG", "LEAD"), each = shift.max),
              rep(1:shift.max, 2), sep = "_"
        )
      long_data_complete[,
        (shift.var.names) :=
          data.table::shift(.SD, n = c(1:shift.max, -(1:shift.max)), type = "cyclic"),
        .SDcols = invariant.vars, by = shift.group
      ]

      for (invar in invariant.vars) {
        for (smx in 1:shift.max) {
          long_data_complete[is.na(get(invar)),
            eval(invar) := get(paste0(invar, "_LAG_", smx))
          ]
          long_data_complete[is.na(get(invar)),
            eval(invar) := get(paste0(invar, "_LEAD_", smx))
          ]
        }
      }
    }
    long_data_complete[, (shift.var.names) := NULL]

    ##    Variable values with conditional statements
    if (!is.null(conditional.vars)) {
      for (cndvar in names(conditional.vars)) {
        shift.var.names <-
          paste(cndvar,
                rep(c("LAG", "LEAD"), each = shift.max),
                rep(1:shift.max, 2), sep = "_"
          )
        excl.expr <-
          paste0("!(", paste(conditional.vars[[cndvar]], collapse = " | "), ")")
        tmp_data_all <- long_data_complete[eval(parse(text = excl.expr)), ]

        for (k in seq_along(conditional.vars[[cndvar]])) {
          cnd <- conditional.vars[[cndvar]][k]
          tmp_data <- long_data_complete[eval(parse(text = cnd)), ]
          tmp_data[, (shift.var.names) := data.table::shift(.SD, n = c(1:shift.max, -(1:shift.max))),
                     .SDcols = cndvar, by = shift.group]

          for (smx in shift.max) {
            tmp_data[is.na(get(cndvar)),
              eval(cndvar) := get(paste0(cndvar, "_LAG_", smx))
            ]
            tmp_data[is.na(get(cndvar)),
              eval(cndvar) := get(paste0(cndvar, "_LEAD_", smx))
            ]
          }

          tmp_data[, (shift.var.names) := NULL]
          tmp_data_all <- data.table::rbindlist(list(tmp_data_all, tmp_data))
        }

        long_data_complete <- tmp_data_all
        rm(tmp_data_all); gc()
      }
    }
    long_data_complete[]
}
