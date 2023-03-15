#' @importFrom data.table rbindlist as.data.table
#' @importFrom parallel makeCluster clusterSetRNGStream clusterExport clusterEvalQ parLapply stopCluster
#' @importFrom mice mice ibind complete
parMICE =
  function(
    data, m, maxit,
    meth, bloks, frmlas,
    visit = NULL,
    seed = NA,
    par.config,
    sample.n
  ) {

    packages <- par.config$packages
    if (is.null(packages)) packages <- "mice"
    if (is.null(par.config$cluster.type)) par.config$cluster.type <- "PSOCK"

    sampleData <-
        function(data, samp.n) {
            cdata <- na.omit(data)
            if (is.character(samp.n)) {
                samp.n <-
                ceiling(
                    nrow(cdata)*as.numeric(gsub("%", "", samp.n)) / 100
                )
            }
            if (samp.n >= nrow(cdata)) return(data)
            data.table::rbindlist(
                list(
                cdata[sample(1:nrow(cdata), samp.n), ],
                na.omit(data, invert = TRUE)
                )
            )
        }

    tmp.cl <-
        parallel::makeCluster(
            spec = par.config$cores,
            type = par.config$cluster.type
        )
    if (!is.na(seed)) {
        parallel::clusterSetRNGStream(cl = tmp.cl, iseed = seed)
    }

    if (maxit == 0) {
        stop("The argument maxit = 0 is not relevant for parallel calculation.")
    }

    parallel::clusterExport(
        cl = tmp.cl,
        varlist =
            list(
                "maxit",
                "meth",
                "bloks",
                "frmlas",
                "visit",
                "packages",
                "sampleData"
            ),
        envir = environment()
    )

    if (par.config$cluster.type != "FORK") {
        parallel::clusterExport(
            cl = tmp.cl,
            varlist = list("data"),
            envir = environment()
        )
    }

    parallel::clusterEvalQ(
        cl = tmp.cl,
        expr = eval(parse(text =
            paste0("require(", packages, ")", collapse = ";")
        ))
    )

    res <-
        parallel::parLapply(
            cl = tmp.cl,
            X = 1:m,
            fun =
                function(f) {
                    res.mice <-
                        mice::mice(
                            data =
                                if (is.null(sample.n)) {
                                    data
                                } else {
                                    sampleData(data, sample.n)
                                },
                            m = 1,
                            maxit = maxit,
                            method = meth,
                            visitSequence = visit,
                            blocks = bloks,
                            formulas = frmlas
                        )
                }
        )

    parallel::stopCluster(cl = tmp.cl)

    if (is.null(sample.n)) {
        res_out <- res[[1]]
        for (i in 2:m) res_out <- mice::ibind(res_out, res[[i]])
    } else {
        res_out <- data[, .imp := 0L]
        for (i in 1:m) {
        res_out <-
            data.table::rbindlist(
                list(
                    res_out,
                    data.table::as.data.table(
                        mice::complete(res[[i]], action = "long"),
                        key = "ID"
                    )[, .id := NULL][, .imp := i]
                ),
                use.names = TRUE
            )
        }
    }
    return(res_out)
}
