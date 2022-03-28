\name{import_log}
\alias{import_log}

\title{
Import MrBayes FBD parameter log files
}
\description{
Imports fossilized birth–death process (FBD) parameter log (.p) files and creates a data frame containing the MCMC samples. Samples can be dropped from the start of each log file (i.e., as burn-in) or downsampled to reduce the size of the output object.
}
\usage{
import_log(path = ".", burnin = 0.25, downsample = 10000)
}
\arguments{
  \item{path}{
the path to a folder containing log (.p) files or a character vector of log files to be read.
}
  \item{burnin}{
either the number of runs or a proportion of runs to drop from the beginning of each log file.
}
  \item{downsample}{
either the number of runs remaining or the proportion of runs remaining after downsampling (i.e., dropping runs from the data). Runs will be dropped in approximately equally-spaced intervals.
}
}
\details{
\code{import_log()} imports FBD phylogenetic analysis log files produced by MrBayes, ignoring the first row of the file (which contains an ID number). The files are appended together, optionally after removing runs at the beginning or throughout the file. When \code{burnin} is greater than 0, the number of or propotion of runs corresponding to the supplied value will be dropped from the beginning of each file as it is read. For example, setting \code{burnin = .25} (the default) will drop the first 25\% of runs from each file. When \code{downsample} is greater than 0, the file will be downsampled (i.e., runs will be dropped) until the number or proportion of runs corresponding to the supplied value is reached. For example, if \code{downsample = 10000} (the default), the returned data frame will contain 10000 samples, selected in approximately equally spaced intervals from the original data.

The log files must contain columns for the FBD parameters \code{"net_speciation"}, \code{"relative_extinction"}, and \code{"relative_fossilization"}.

After reading in the parameter samples the data are coverted to "long" format with one row per generation per time bin and a column for each parameter.

The output can be supplied to \code{\link{FBD_summary}}, \code{\link{FBD_dens_plot}}, \code{\link{FBD_normality_plot}}, \code{\link{FBD_tests1}}, or \code{\link{FBD_tests2}}.
}
\value{
A data frame with columns corresponding to the columns in the supplied log files and rows containing the sampled parameter values, reshaped into longer format with one row per generation per time bin and a new column \code{"Time_bin"} containing the time bins as a factor. An example of the kind of output produced can be accessed using \code{\link[=posterior1p]{data("posterior1p")}}.
}
\examples{
\dontrun{
posterior <- import_log("path/to/folder", burnin = .25,
                        downsample = 1e4)
}
}

