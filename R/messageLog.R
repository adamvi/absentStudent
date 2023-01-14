#' @title A function for creating customized log reports and saving message info.
#' @description Function produces a text file (custom message log).
#' @param log.message Text string to be exported to a text file log.
#' @param log.directory File path where log file should be saved.
#'  Default: 'Logs'
#' @param logfile Name of file to be saved.  If NULL, then file will be named
#'  as \code{messageLog_*DATE*.txt}, where \code{*DATE*} is the current date.
#'  Default: NULL
#' @param add.date Add the date to a user provided (non-NULL) \code{logfile} name?
#'  Default: TRUE
#' @param appendLF Add log message to a file (if it already exists)?
#'  Default: TRUE
#' @return A text file with name provided in "logfile", saved in the "log.directory".
#' @details Creates (or adds to) a log file with the text given in the "log.message" argument.
#' @rdname messageLog
#' @export 
#' @author Adam R. VanIwaarden \email{avaniwaarden@nciea.org}
#' @keywords misc
#' @keywords models

messageLog =
  function(
	log.message,
	log.directory = "Logs",
	logfile = NULL,
	add.date = TRUE,
	appendLF = TRUE
  ) {

	PrintLogMessage <- function() {
		# print log message to file
		if (!dir.exists(log.directory)) dir.create(log.directory, recursive = TRUE)
		if (is.null(logfile)) {
			logfile <- file.path(log.directory, paste0("messageLog_", gsub("-", "_", Sys.Date()), ".txt"))
		} else {
			if (add.date) {
				logfile <- paste0(logfile, "_", paste(strsplit(as.character(Sys.Date()), "-")[[1]][c(2,3,1)], collapse="_"), ".txt")
			} else logfile <- paste0(logfile, ".txt")
			logfile <- file.path(log.directory, logfile)
		}

		if (is.call(log.message)) {
			log.message2 <- c(paste0("\n\n\t", as.character(log.message)[1L], "(\n\t\t"), paste(names(log.message)[-1L], as.character(log.message)[-1L], sep=" = ", collapse="\n\t\t"), ")\n\n")
			cat(log.message2, file = logfile, append=appendLF)
		} else cat(log.message, "\n", file=logfile, sep="", append=appendLF)
	}

	if (!is.call(log.message)) {
		message(log.message)
	}
	PrintLogMessage()
	invisible()
  }
