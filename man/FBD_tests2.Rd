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
A data frame of posterior parameter estimates containing a single "Time_bin" column and one column for each FBD parameter value (e.g.,  "net_speciation", "relative extinction", and "relative fossilization"). Such data frame can be imported using \code{\link{combine_log}} followed by \code{\link{FBD_reshape}} . See the \code{\link{posterior1p_long}} or \code{\link{posterior3p_long}} datasets for an examples of how the input file should look.
}
  \item{p.adjust.method}{
The method use to adjust the p-values for multiple testing. See \code{\link{p.adjust}} for details and options. Default if \code{"fdr"} for the Benjamini-Hochberg false discovery rate correction.
}
}
\details{
\code{\link{pairwise.t.test}} and \code{\link{pairwise.wilcox.test}} are used to calculate, respectively, the t-test and Mann-Whitney U-tests statistics and p-values. Because the power of these tests depends on the number of posterior samples, it can be helpful to examine the distributions of FBD parameter posteriors using \code{\link{FBD_dens_plot}} instead of relying heavily on the tests.
}
\value{
A list with an element for each test, each of which contains a list of test results for each parameter. The results are in the form of a data frame containing the sample sizes and unadjusted and adjusted p-values for each comparison.
}

\seealso{
\code{\link{combine_log}} for producing a single data set of parameter posterior samples from individual parameter log files.

\code{\link{FBD_reshape}} for converting posterior parameter table from wide to long format.

\code{\link{FBD_dens_plot}}, \code{\link{FBD_normality_plot}}, \code{\link{FBD_tests1}}, and \code{\link{FBD_tests2}} for other functions used to summarize and display the distributions of the parameter posteriors.

\code{\link{pairwise.t.test}} and \code{\link{pairwise.wilcox.test}} for the tests used.
}
\examples{
data("posterior3p_long")

FBD_tests2(posterior3p_long)
}
