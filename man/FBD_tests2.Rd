\name{FBD_tests2}
\alias{FBD_tests2}

\title{
Tests for differences in FBD parameter values
}
\description{
\code{FBD_tests2()} performs t-tests and Mann-Whitney U-tests to compare the average value of fossilized birth–death process (FBD) parameters between time bins.
}
\usage{
FBD_tests2(posterior, p.adjust.method = "fdr")
}
\arguments{
  \item{posterior}{
a dataset of FBD parameter MCMC posterior samples imported using \code{\link{import_log}}. See the \code{\link{posterior1p}} dataset for an example of how the resulting file should look.
}
  \item{p.adjust.method}{
the method use to adjust the p-values for multiple testing. See \code{\link{p.adjust}} for details and options. Default if \code{"fdr"} for the Benjamini-Hochberg false discovery rate correction.
}
}
\details{
\code{\link{pairwise.t.test}} and \code{\link{pairwise.wilcox.test}} are used to compute the t-test and Mann-Whitney U-tests statistics and p-values. Because the power of these tests depends on the number of MCMC samples, it can be helpful to examine the distributions of FBD parameter posteriors using \code{\link{FBD_dens_plot}} instead of relying heavily on the tests.
}
\value{
A list with an element for each test, each of which contains a list of test results for each parameter. The results are in the form of a data frame containing the sample sizes and unadjusted and adjusted p-values for each comparison.
}

\seealso{
\code{\link{import_log}} for producing a data set of FBD parameter posterior samples from parameter log files.

\code{\link{FBD_dens_plot}}, \code{\link{FBD_normality_plot}}, \code{\link{FBD_tests1}}, and \code{\link{FBD_tests2}} for other functions used to summarize and display the distributions of the parameter posteriors.

\code{\link{pairwise.t.test}} and \code{\link{pairwise.wilcox.test}} for the tests used.
}
\examples{
data("posterior1p")

FBD_tests2(posterior1p)
}
