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

##  timing functions
`timeTaken` =
  function(started.at) {
    format = function(secs) {
      secs.integer = as.integer(secs)
      sprintf("%02d:%02d:%02d:%.3f",
        secs.integer%/%86400L,
        (secs.integer%/%3600L)%%24L,
        (secs.integer%/%60L)%%60L,
        secs%%60L)
    }
    tt = proc.time() - started.at
    format(tt[3L])
  } ###   From SGP package: timetakenSGP function

`convertTime` =
function(tmp.time) {
  tmp <- tail(c(0, 0, 0, as.numeric(unlist(strsplit(tmp.time, ":")))), 4)
  tmp.label <- c("Day", "Hour", "Minute", "Second")
  tmp.label[which(tmp!=1)] <- paste0(tmp.label, "s")[which(tmp!=1)]
  return(paste(paste(tmp[tmp!=0], tmp.label[tmp!=0]), collapse=", "))
} ###   From SGP package: convertTime
