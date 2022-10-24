\name{combine_log}
\alias{combine_log}

\title{
Combine and filter (.p) log files from Mr.Bayes
}
\description{
Imports parameter (.p) log files from Mr. Bayes and combines them into a single data frame. Samples can be dropped from the start of each log file (i.e., discarded as burn-in) and/or downsampled to reduce the size of the output object.
}
\usage{
combine_log(path = ".", burnin = 0.25, downsample = 10000)
}
\arguments{
  \item{path}{
The path to a folder containing (.p) log files or a character vector of log files to be read.
}
  \item{burnin}{
Either the number or a proportion of generations to drop from the beginning of each log file.
}
  \item{downsample}{
Either the number or the proportion of generations the user wants to keep after downsampling for the final (combined) log file. Generations will be dropped in approximately equally-spaced intervals.
}
}
\details{
\code{combine_log()} imports log files produced by Mr.Bayes, ignoring the first row of the file (which contains an ID number). The files are appended together, optionally after removing burn-in generations from the beginning and/or by further filtering throughout the rest of each file. When \code{burnin} is greater than 0, the number or propotion of generations corresponding to the supplied value will be dropped from the beginning of each file as it is read in. For example, setting \code{burnin = .25} (the default) will drop the first 25\% of generations from each file. When \code{downsample} is greater than 0, the file will be downsampled until the number or proportion of generations corresponding to the supplied value is reached. For example, if \code{downsample = 10000} generations (the default) for log files from 4 independent runs (i.e., 4 (.p) files), each log file will be downsampled to 2500 generations, and the final combined data frame will contain 10000 samples, selected in approximately equally spaced intervals from the original data.

The output can be supplied to \code{\link{get_pwt_rates_MrBayes}} and to \code{\link{FBD_reshape}}. The latter will convert the log data frame from my wide to long format, which is necessary to be used as input for downstream analyses using \code{\link{FBD_summary}}, \code{\link{FBD_dens_plot}}, \code{\link{FBD_normality_plot}}, \code{\link{FBD_tests1}}, or \code{\link{FBD_tests2}}.

}
\value{
A data frame with columns corresponding to the columns in the supplied log files and rows containing the sampled parameter values. Examples of the kind of output produced can be accessed using \code{\link[=posterior1p]{data("posterior1p")}} and \code{\link[=posterior3p]{data("posterior3p")}}.
}
\seealso{
\code{vignette("fbd-params")} for the use of this function as part of an analysis pipeline.

\code{\link{FBD_reshape}}, which reshapes a combined parameter log file for use in some other package functions.
}
\examples{
# See vignette("fbd-params") for how to use this
# function as part of an analysis pipeline
\dontrun{
posterior <- combine_log("path/to/folder", burnin = .25,
                         downsample = 10000)
}
}

