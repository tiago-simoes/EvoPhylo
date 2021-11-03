\name{FBD_tests1}
\alias{FBD_tests1}

\title{
Test assumptions of normality and homoscedasticity for parameters
}
\description{
Produces tests of normality (within time bin, ignoring time bin, and pooling within-time bin values) and homoscedasticity (equality of variances) for each parameter.
}
\usage{
FBD_tests1(AllRunsrMelted_MC, downsample = TRUE)
}

\arguments{
  \item{AllRunsrMelted_MC}{
an AllRunsMelted_MC (long) dataset with a Time_bin column and columns for each parameter (net_speciation, relative_exinction, and relative_fossilization).
}
  \item{downsample}{
whether to downsample the observations to ensure Shapiro-Wilk normality tests can be run. If \code{TRUE}, observations will be dropped so that no more than 5000 observations are used for the tests on the full dataset, as required by \code{\link{shapiro.test}}. They will be dropped in evenly spaced intervals. If \code{FALSE} and there are more than 5000 observations for any test, that test will not be run.
}
}
\details{
\code{FBD_tests1()} performs several tests on the distributions of parameter values within and across time bins. It produces the Shapiro-Wilk test for normality using \code{\link{shapiro.test}} and the Bartlett and Fligner tests for homogeneity of variance using \code{\link{bartlett.test}} and \code{\link{fligner.test}}, respectively. Note that these tests are likely to be significant even if the observations are approximately normally distributed or have approximately equal variance; therefore, they should be supplemented with visual inspection using \code{\link{FBD_normality_plot}}.
}
\value{
A list containing the results of the three tests with the following elements:
\item{shapiro}{A list with an element for each parameter. Each element is a data frame with a row for each time bin and the test statistic and p-value for the Shipro-Wilk test for normality. In addition, there will be a row for an overall test, combining all observations ignoring time bin, and a test of the residuals, which combines the group-mean-centered observations (equivalent to the residuals in a regression of the parameter on time bin).}
\item{bartlett}{A data frame of the Bartlett test for homogeneity of variance across time bins with a row for each parameter and the test statistic and p-value for the test.}
\item{fligner}{A data frame of the Fligeer test for homogeneity of variance across time bins with a row for each parameter and the test statistic and p-value for the test.}
}

\seealso{
\code{\link{FBD_normality_plot}} for visual tests.

\code{\link{FBD_tests2}} for tests of differences in parameter means.

\code{\link{shapiro.test}}, \code{\link{bartlett.test}}, and \code{\link{fligner.test}} for the statistical tests used.
}
\examples{

}
