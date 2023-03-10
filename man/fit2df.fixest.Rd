% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fixest_finalfit.R
\name{fit2df.fixest}
\alias{fit2df.fixest}
\title{FUNCTION_TITLE}
\usage{
fit2df.fixest(
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
)
}
\arguments{
\item{.data}{PARAM_DESCRIPTION}

\item{condense}{PARAM_DESCRIPTION, Default: TRUE}

\item{metrics}{PARAM_DESCRIPTION, Default: FALSE}

\item{remove_intercept}{PARAM_DESCRIPTION, Default: TRUE}

\item{explanatory_name}{PARAM_DESCRIPTION, Default: 'explanatory'}

\item{estimate_name}{PARAM_DESCRIPTION, Default: 'OR'}

\item{estimate_suffix}{PARAM_DESCRIPTION, Default: ''}

\item{p_name}{PARAM_DESCRIPTION, Default: 'p'}

\item{digits}{PARAM_DESCRIPTION, Default: c(2, 2, 3)}

\item{exp}{PARAM_DESCRIPTION, Default: TRUE}

\item{confint_type}{PARAM_DESCRIPTION, Default: 'profile'}

\item{confint_level}{PARAM_DESCRIPTION, Default: 0.95}

\item{confint_sep}{PARAM_DESCRIPTION, Default: '-'}

\item{...}{PARAM_DESCRIPTION}
}
\value{
OUTPUT_DESCRIPTION
}
\description{
FUNCTION_DESCRIPTION
}
\details{
DETAILS
}
\examples{
\dontrun{
if(interactive()){
 #EXAMPLE1
 }
}
}
\seealso{
\code{\link[finalfit]{condense_fit}}, \code{\link[finalfit]{remove_intercept}}
}
\author{
Adam R. VanIwaarden \email{avaniwaarden@nciea.org}
}
