\name{get_clockrate_table}
\alias{get_clockrate_table}

\title{
Extract clock rates from a Tree
}
\description{
Extracts the clock rates for each node from a tree and stores them in a table.
}
\usage{
get_clockrate_table(tree, summary = "mean",
                    drop.dummyextent = TRUE, file = NULL)
}
\arguments{
  \item{tree}{
a \code{tidytree} object; the output of a call to \code{\link[treeio:read.mrbayes]{treeio::read.mrbayes}}.
}
  \item{summary}{
the name of the rate summary. Should be one of \code{"mean"} or \code{"median"}.
}
  \item{drop.dummyextent}{
\code{logical}; whether to drop the "DummyExtent" tip of the tree prior to extracting the clock rates. Default is \code{TRUE}.
}
  \item{file}{
an optional file path where the resulting table will be stored using \code{\link{write.csv}}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A data frame containing a column containing the node identifer (\code{node}) and one column for each clock in the tree object containing the rates.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
