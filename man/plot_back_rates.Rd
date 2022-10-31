\name{plot_back_rates}
\alias{plot_back_rates}

\title{
Plots distribution of background rates extracted from posterior log files.
}
\description{
Plots The distribution of background rates extracted from the posterior log files from Mr. Bayes or BEAST2, as well as the distribution of background rates if log transformed to test for normality of data distribution.
}
\usage{
plot_back_rates(type = c("MrBayes", "BEAST2"),
                posterior,
                clock = 1,
                trans = c("none", "log", "log10"),
                size = 12, quantile = 0.95)
}
\arguments{
    \item{type}{
Whether to use data output from "Mr.Bayes" or "BEAST2".
}
  \item{posterior}{
A data frame of posterior parameter estimates (log file). From Mr.Bayes, it includes a "clockrate" column indicating the mean (background) clock rate estimate for each generation that will be used for pairwise t-tests. Such data frame can be imported using \code{\link{combine_log}} (no need to reshape from wide to long). See the \code{\link{posterior1p}} or \code{\link{posterior3p}} datasets for an examples of how the input file should look.
     From BEAST2, it will include at least one "rate<filename>.mean" column indicating the mean (background) clock rate estimate for each generation. If there are "P" unlinked clock partitions in BEAST2, there will be P x "rate<filename>.mean" columns (one for each partition) in the posterior log file.
}
  \item{clock}{
The clock partition number to calculate selection mode. Ignored if only one clock is available.
}
  \item{trans}{
Type of data transformation to perform on background rates extracted from the posterior log file from Mr. Bayes or BEAST2. Options include "none" (if rates are normally distributed), natural log transformation "log", and log of base 10 transformation "log10". The necessity of using data transformation can be tested using the function \code{\link{plot_back_rates}}.
  }
  \item{size}{
Font size for title of plot
}
  \item{quantile}{
Upper limit for X axis (passed on to 'xlim') to remove outliers from histogram. The quantile can be any value between "0" and "1", but values equal or above "0.95" provide good results in most cases in which the data distribution is right skewed.
}
}
\details{
Plots The distribution of background rates extracted from the posterior log files from Mr. Bayes or BEAST2, as well as the distribution of background rates if log transformed.  Background rates should be normally distributed for meeting the assumptions of t-tests and other tests passed on by downstream functions, including \code{\link{get_pwt_rates_MrBayes}}, \code{\link{get_pwt_rates_BEAST2}}, and \code{\link{plot_treerates_sgn}}.
}

\value{
It produces a \code{ggplot} object that can be manipulated using \pkg{ggplot2} syntax (e.g., to change the \code{theme} or labels).
}

\seealso{
\code{vignette("rates-selection")} for the use of this function as part of an analysis pipeline.
}
\examples{
# See vignette("rates-selection") for how to use this
# function as part of an analysis pipeline

## MrBayes example
# Load example tree and posterior

data("posterior3p")

P <- plot_back_rates (type = "MrBayes", posterior, clock = 1,
                      trans = "log10", size = 10, quantile = 0.95)
P
}

