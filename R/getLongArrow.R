
#' @importFrom data.table data.table as.data.table setcolorder setkeyv rbindlist
#' @importFrom arrow arrow_table
#' @importFrom dplyr semi_join
getLongArrow =
    function(
        db.con,
        config,
        cols = NULL
    ){
        if (length(config) > 1) return_data <- data.table()

        for (K in seq(config)) {
            config.iter <- config[[K]]
            names(config.iter) <- gsub("^sgp[.]", "", names(config.iter))

            if ("list" %in% class(config.iter[["grade.sequences"]])) {
                for (G in seq(config.iter[["grade.sequences"]])) {
                    cohort.iter <- config.iter
                    cohort.iter[["grade.sequences"]] <-
                        cohort.iter[["grade.sequences"]][[G]]
                    grade.length <- length(cohort.iter[["grade.sequences"]])
                    cohort_lookup <-
                        data.table::data.table(
                            VALID_CASE = "VALID_CASE",
                            CONTENT_AREA =
                                tail(cohort.iter[["content.areas"]], grade.length),
                            YEAR =
                                tail(cohort.iter[["panel.years"]], grade.length),
                            GRADE =
                                cohort.iter[["grade.sequences"]]
                        ) |>
                            arrow::arrow_table()

                    cohort_lookup <-
                        cohort_lookup$cast(
                            target_schema =
                                db.con$schema[c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE")]
                        )

                    cohort_data <- db.con |>
                        dplyr::semi_join(cohort_lookup) |>
                            data.table::as.data.table(key = getKey(db.con)) |>
                                data.table::setcolorder(getKey(db.con))

                    if (!is.null(cols)) {
                        cohort_data <-
                            cohort_data[,
                            (unique(c(getKey(db.con), cols))),
                            with = FALSE
                            ]
                    }
                    if (length(config) > 1) {
                        data.table::setkeyv(cohort_data, getKey(cohort_data))
                        return_data <-
                            data.table::rbindlist(list(return_data, cohort_data))
                    }
                }
            } else {
                if ("list" %in% class(config.iter)) {
                    grade.length <- length(config.iter[["grade.sequences"]])
                    cohort_lookup <-
                        data.table::data.table(
                            VALID_CASE = "VALID_CASE",
                            CONTENT_AREA =
                                tail(config.iter[["content.areas"]], grade.length),
                            YEAR =
                                tail(config.iter[["panel.years"]], grade.length),
                            GRADE =
                                config.iter[["grade.sequences"]]
                        ) |>
                            arrow::arrow_table()
                } else {
                    cohort_lookup <- arrow::arrow_table(config.iter)
                }
                cohort_lookup <-
                    cohort_lookup$cast(
                        target_schema =
                            db.con$schema[c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE")]
                    )

                cohort_data <- db.con |>
                    dplyr::semi_join(cohort_lookup) |>
                        data.table::as.data.table(key = getKey(db.con)) |>
                            data.table::setcolorder(getKey(db.con))

                if (!is.null(cols)) {
                    cohort_data <-
                        cohort_data[,
                          (unique(c(getKey(db.con), cols))),
                          with = FALSE
                        ]
                }
                if (length(config) > 1) {
                    # data.table::setkeyv(cohort_data, getKey(cohort_data))
                    return_data <-
                        data.table::rbindlist(list(return_data, cohort_data))
                }
            }
        }
        if (length(config) > 1) return_data[] else cohort_data[]
    }
