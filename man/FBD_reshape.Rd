\name{FBD_reshape}
\alias{FBD_reshape}

\title{
Reshape Allruns data from wide to long
}
\description{
Reshapes Allruns data from wide, where each time bin has its own column for each parameter, to long, where time bin is stored in its own column and each row corresponds to a parameter value at a given time bin. Essentially a wrapper for \code{\link{reshape}}.
}
\usage{
FBD_reshape(AllRuns)
}
\arguments{
  \item{AllRuns}{
An Allrund dataset with a column for each parameter (net_speciation, relative_exinction, and relative_fossilization) at each time bin.
}
}
\details{

}
\value{
A long dataset with a new "Time_bin" column and one column for each parameter. The \code{"reshapeLong"} attribute contains information for automatically creating a wide dataset using \code{reshape(out, direction = "wide")}.
}

\author{

}
\note{

}

\seealso{
\code{\link{reshape}}
}
\examples{

}
