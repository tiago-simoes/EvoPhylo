\name{FBD_tests2}
\alias{FBD_tests2}

\title{
Tests for differences in parameter values
}
\description{
\code{FBD_tests2()} performs t-tests and Mann-Whitney U-tests to compare the average parameter value between time bins.
}
\usage{
FBD_tests2(AllRunsrMelted_MC)
}
\arguments{
  \item{AllRunsrMelted_MC}{
an AllRunsMelted_MC (long) dataset with a Time_bin column and columns for each parameter (net_speciation, relative_exinction, and relative_fossilization).
}
}
\details{
\code{\link{pairwise.t.test}} and \code{\link{pairwise.wilcox.test}} are used to compute the t-test and Mann-Whitney U-tests statistics and p-values. The p-values within each parameter and test are adjusted using the Benjamini-Hochberg false discovery rate correction using \code{\link{p.adjust}} with \code{method = "fdr"}.
}
\value{
A list with an element for each test, each of which contains a list of test results for each parameter. The results are in the form of a data frame containing the sample sizes and unadjusted and adjusted p-values for each comparison.
}

\seealso{
\code{\link{FBD_tests1}} for tests of assumptions on the distributions of the parameters.

\code{\link{pairwise.t.test}} and \code{\link{pairwise.wilcox.test}} for the tests used.
}
\examples{

}
