###  Utility functions from SGP package
`%w/o%` = function(x, y) x[!x %in% y]

`getKey` = function(data) {
  if ("YEAR_WITHIN" %in% names(data)) {
    return(c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE", "ID", "YEAR_WITHIN"))
  } else {
    return(c("VALID_CASE", "CONTENT_AREA", "YEAR", "GRADE", "ID"))
  }
}

`Z` = function(data_table, var.to.standardize, reference.year = NULL, rm.na = TRUE) {
  YEAR <- NULL
  x <- data_table[, get(var.to.standardize)]
  if (!is.null(reference.year)) {
    y <- data_table[YEAR == reference.year, get(var.to.standardize)]
  } else {
    y <- x
  }
  (x - mean(y, na.rm = rm.na)) / stats::sd(y, na.rm = rm.na)
}

##  parMICE Internal function
`miceInt` = function(data, maxit, meth, visit, bloks, frmlas) {
    mice::mice(
        data = data,
        m = 1,
        maxit = maxit,
        method = meth,
        visitSequence = visit,
        blocks = bloks,
        formulas = frmlas
    )
}

