#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data PARAM_DESCRIPTION
#' @param m PARAM_DESCRIPTION, Default: 5
#' @param maxit PARAM_DESCRIPTION, Default: 5
#' @param meth PARAM_DESCRIPTION
#' @param bloks PARAM_DESCRIPTION
#' @param frmlas PARAM_DESCRIPTION
#' @param visit PARAM_DESCRIPTION, Default: NULL
#' @param seed PARAM_DESCRIPTION, Default: NA
#' @param nnodes PARAM_DESCRIPTION, Default: 5
#' @param cluster.type PARAM_DESCRIPTION, Default: NULL
#' @param packages PARAM_DESCRIPTION, Default: NULL
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[parallel]{makeCluster}}, \code{\link[parallel]{RNGstreams}}, \code{\link[parallel]{clusterApply}}
#'  \code{\link[mice]{mice}}
#'  \code{\link[abind]{abind}}
#' @rdname parMICE
#' @export 
#' @author AUTHOR [AUTHOR_2]
#' @keywords KEYWORD_TERM
#' @author AUTHOR [AUTHOR_2]
#' @keywords KEYWORD_TERM
#' @importFrom parallel makeCluster clusterSetRNGStream clusterExport clusterEvalQ parSapply stopCluster
#' @importFrom mice mice
#' @importFrom abind abind
parMICE =
  function(
    data,
    m = 5,
    maxit = 5,
    meth,
    bloks,
    frmlas,
    visit = NULL,
    seed = NA,
    nnodes = 5,
    cluster.type = NULL,
    packages = NULL
  ) {

    if (missing(bloks) || missing(frmlas)) {
        stop("The `bloks` OR `frmlas` arguments are not defined. Currently, this case is not handled by the parMICE function")
    }

    if (is.null(cluster.type)) cluster.type <- "PSOCK"

    tmp.cl <- parallel::makeCluster(spec = nnodes, type = cluster.type)
    if (!is.na(seed)) {
        parallel::clusterSetRNGStream(cl = tmp.cl, iseed = seed)
    }

    if (maxit == 0) {
        stop("The argument maxit = 0 is not relevant for parallel calculation, use the mice function from the mice package")
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
                "packages"
            ),
        envir = environment()
    )

    if (cluster.type != "FORK") {
        parallel::clusterExport(
            cl = tmp.cl,
            varlist = list("data"),
            envir = environment()
        )
    }

    parallel::clusterEvalQ(
        cl = tmp.cl,
        expr = eval(parse(text = paste0("require(", packages, ")", collapse = ";")))
    )

    res <-
        parallel::parSapply(
            cl = tmp.cl,
            X = 1:m,
            FUN =
                \(f) {
                res.mice <-
                    mice::mice(
                        data = data,
                        m = 1,
                        maxit = maxit,
                        method = meth,
                        visitSequence = visit,
                        blocks = bloks,
                        formulas = frmlas
                    )
                },
            simplify = FALSE
        )

    parallel::stopCluster(cl = tmp.cl)

    res.out <- res[[1]]
    res.out$call <- match.call(expand.dots = FALSE)
    res.out$m <- m
    res.out$imp <-
        mapply(
          as.list(colnames(data)),
          FUN =
              \(xx, res) {
                do.call(cbind, lapply(lapply(res, "[[", "imp"), "[[", xx))
              },
          MoreArgs = list(res = res),
          SIMPLIFY = FALSE
        )
    names(res.out$imp) <- colnames(data)
    res.out$imp <-
        lapply(
            res.out$imp,
            \(xx) {
                if (!is.null(xx)) {
                    yy <- xx
                    colnames(yy) <- as.character(seq(ncol(xx)))
                    return(yy)
                } else {
                    return(xx)
                }
            }
        )
    res.out$seed <- seed
    res.out$lastSeedValue <- lapply(res, "[[", "lastSeedValue")
    res.out$chainMean <-
        do.call(abind::abind, lapply(res, "[[", "chainMean"), 3)
    dimnames(res.out$chainMean)[[3]] <- paste("Chain", seq(m))
    res.out$chainVar <-
        do.call(abind::abind, lapply(res, "[[", "chainVar"), 3)
    dimnames(res.out$chainVar)[[3]] <- paste("Chain", seq(m))
    res.out$loggedEvents <- lapply(res, "[[", "loggedEvents")
    class(res.out) <- "mids"
    return(res.out)
}
