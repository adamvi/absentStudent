#' @importFrom crayon magenta
#' @importFrom utils packageVersion

`.onAttach` =
  function(libname, pkgname) {
    if (interactive()) {
      packageStartupMessage(crayon::magenta$bold("absentStudent",
        paste(paste0(unlist(strsplit(as.character(utils::packageVersion("absentStudent")), "[.]")),
        c(".", "-", ".", "")), collapse = ""),
        " (5-30-2023). For help visit https://adamvi.github.io/absentStudent"))
  }
}
