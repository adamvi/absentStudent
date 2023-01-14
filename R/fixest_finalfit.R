#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param .data PARAM_DESCRIPTION
#' @param condense PARAM_DESCRIPTION, Default: TRUE
#' @param metrics PARAM_DESCRIPTION, Default: FALSE
#' @param remove_intercept PARAM_DESCRIPTION, Default: TRUE
#' @param explanatory_name PARAM_DESCRIPTION, Default: 'explanatory'
#' @param estimate_name PARAM_DESCRIPTION, Default: 'OR'
#' @param estimate_suffix PARAM_DESCRIPTION, Default: ''
#' @param p_name PARAM_DESCRIPTION, Default: 'p'
#' @param digits PARAM_DESCRIPTION, Default: c(2, 2, 3)
#' @param exp PARAM_DESCRIPTION, Default: TRUE
#' @param confint_type PARAM_DESCRIPTION, Default: 'profile'
#' @param confint_level PARAM_DESCRIPTION, Default: 0.95
#' @param confint_sep PARAM_DESCRIPTION, Default: '-'
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
#'  \code{\link[finalfit]{condense_fit}}, \code{\link[finalfit]{remove_intercept}}
#' @rdname fit2df.fixest
#' @author Adam R. VanIwaarden \email{avaniwaarden@nciea.org}
#' @importFrom finalfit condense_fit remove_intercept
fit2df.fixest =
  function(
    .data,
    condense = TRUE,
    metrics = FALSE,
    remove_intercept = TRUE,
    explanatory_name = "explanatory",
    estimate_name = "OR",
    estimate_suffix = "",
    p_name = "p",
    digits = c(2, 2, 3),
    exp = TRUE,
    confint_type = "profile",
    confint_level = 0.95,
    confint_sep = "-",
    ...
  ) {
    df.out <- extract_fit.fixest(
        .data = .data, explanatory_name = explanatory_name,
        estimate_name = estimate_name, estimate_suffix = estimate_suffix,
        exp = exp,
        confint_type = confint_type,
        confint_level = confint_level,
        p_name = p_name
    )

    if (condense == TRUE) {
        df.out <- finalfit::condense_fit(df.out,
            explanatory_name = explanatory_name,
            estimate_name = estimate_name, estimate_suffix = estimate_suffix,
            p_name = p_name, digits = digits, confint_sep = confint_sep
        )
    }

    if (remove_intercept == TRUE) {
        df.out <- finalfit::remove_intercept(df.out)
    }

    # Extract model metrics
    if (metrics == TRUE) {
        metrics.out <- ff_metrics.fixest(.data)
        return(list(df.out, metrics.out))
    } else {
        return(df.out)
    }
  }

#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param .data PARAM_DESCRIPTION
#' @param explanatory_name PARAM_DESCRIPTION, Default: 'explanatory'
#' @param estimate_name PARAM_DESCRIPTION, Default: 'OR'
#' @param estimate_suffix PARAM_DESCRIPTION, Default: ''
#' @param p_name PARAM_DESCRIPTION, Default: 'p'
#' @param exp PARAM_DESCRIPTION, Default: TRUE
#' @param confint_type PARAM_DESCRIPTION, Default: 'profile'
#' @param confint_level PARAM_DESCRIPTION, Default: 0.95
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
#'  \code{\link[dplyr]{reexports}}, \code{\link[dplyr]{select}}
#' @rdname extract_fit.fixest
#' @author Adam R. VanIwaarden \email{avaniwaarden@nciea.org}
#' @importFrom dplyr tibble select
extract_fit.fixest =
  function(
    .data,
    explanatory_name = "explanatory",
    estimate_name = "OR",
    estimate_suffix = "",
    p_name = "p",
    exp = TRUE,
    confint_type = "profile",
    confint_level = 0.95,
    ...
  ) {
    x <- .data
    explanatory <- names(coef(x))
    estimate <- coef(x)
    if (confint_type == "profile") {
        confint <- stats::confint(x, level = confint_level)
    } else if (confint_type == "default") {
        confint <- stats::confint.default(x, level = confint_level)
    }
    p_col <- dimnames(x$coeftable)[[2]] %in% c("Pr(>|t|)", "Pr(>|z|)")
    p <- x$coeftable[, p_col]
    L_confint_name <- paste0("L", confint_level * 100)
    U_confint_name <- paste0("U", confint_level * 100)

    df.out <- dplyr::tibble(explanatory, estimate, confint[, 1], confint[, 2], p)
    colnames(df.out) <- c(
        explanatory_name, paste0(estimate_name, estimate_suffix),
        L_confint_name, U_confint_name, p_name
    )
    if (exp) {
        df.out[, 2:4] <- exp(df.out[, 2:4])
    }
    if (confint_level != 0.95) {
        df.out <- df.out %>% dplyr::select(-p_name)
    }
    df.out <- data.frame(df.out)
    return(df.out)
  }


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param .data PARAM_DESCRIPTION
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso 
#'  \code{\link[stats]{AIC}}
#'  \code{\link[pROC]{are.paired}}, \code{\link[pROC]{auc}}, \code{\link[pROC]{ci}}, \code{\link[pROC]{ci.auc}}, \code{\link[pROC]{ci.coords}}, \code{\link[pROC]{ci.se}}, \code{\link[pROC]{ci.sp}}, \code{\link[pROC]{ci.thresholds}}, \code{\link[pROC]{coords}}, \code{\link[pROC]{cov.roc}}, \code{\link[pROC]{has.partial.auc}}, \code{\link[pROC]{lines.roc}}, \code{\link[pROC]{multiclass.roc}}, \code{\link[pROC]{pROC-package}}, \code{\link[pROC]{plot.ci}}, \code{\link[pROC]{plot.roc}}, \code{\link[pROC]{power.roc.test}}, \code{\link[pROC]{print}}, \code{\link[pROC]{roc}}, \code{\link[pROC]{roc.test}}, \code{\link[pROC]{smooth}}, \code{\link[pROC]{var.roc}}
#'  \code{\link[finalfit]{metrics_hoslem}}
#' @rdname ff_metrics.fixest
#' @author Adam R. VanIwaarden \email{avaniwaarden@nciea.org}
#' @importFrom stats BIC
#' @importFrom pROC roc
#' @importFrom finalfit metrics_hoslem
ff_metrics.fixest =
  function(.data){
    x = .data
    n_data = x$nobs_origin
    n_model = x$nobs
    n_groups = x$fixef_sizes
    log_like = round(x$loglik, 1)
    pseudo_r2 = round(x$pseudo_r2, 3)
    sq_cor = round(x$sq.cor, 3)
    bic = round(stats::BIC(x), 1)
    if (!"lean" %in% names(x)) {
        auc = round(pROC::roc(x$y, x$fitted.values)$auc[1], 3)
        h_l = finalfit::metrics_hoslem(x$y, x$fitted.values)
    } else {
        h_l <- auc <- "NA - 'lean' model object provided"
    }
    metrics.out = paste0(
        "Number in dataframe = ", n_data,
        ", Number in model = ", n_model,
        ", Missing = ", n_data-n_model,
        ", BIC = ", bic,
        ", Log-Likelihood = ", log_like,
        ", Adj. Pseudo R2 = ", pseudo_r2,
        ", Squared Cor. = ", sq_cor,
        ", C-statistic = ", auc,
        ", H&L = ", h_l) |> 
        as.data.frame(stringsAsFactors = FALSE) |> 
        unname()
    class(metrics.out) = c("data.frame.ff", class(metrics.out))
    return(metrics.out)
  }


#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param .data PARAM_DESCRIPTION
#' @param dependent PARAM_DESCRIPTION
#' @param explanatory PARAM_DESCRIPTION
#' @param fixed_effect PARAM_DESCRIPTION, Default: NULL
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
#'  \code{\link[stats]{formula}}
#'  \code{\link[fixest]{feglm}}
#' @rdname glmfixed
#' @author Adam R. VanIwaarden \email{avaniwaarden@nciea.org}
#' @importFrom stats as.formula
#' @importFrom fixest feglm
glmfixed =
  function(
    .data,
    dependent,
    explanatory,
    fixed_effect = NULL,
    ...
  ) {
    if (length(fixed_effect)) {
        if (!grepl("\\|", fixed_effect)) {
            fixed_effect <- paste0(" | ", fixed_effect)
        }
    }
    ## fixest::feglm requires DV be 0 or 1
    if (is.factor(.data[[dependent]])) {
        .data[[dependent]] <- as.integer(.data[[dependent]]) - 1L
    }
    mod.fmla <- stats::as.formula(
        paste0(dependent, " ~ ",
               paste(explanatory, collapse = " + "),
               fixed_effect
        )
    )
    fixest::feglm(fml = mod.fmla, data = .data, family = "binomial", ...)
  }
