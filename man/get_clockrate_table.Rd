\name{get_clockrate_table}
\alias{get_clockrate_table}

\title{
Extract clock rates from a Tree
}
\description{
Extracts the clock rates for each node from a tree and stores them in a table.
}
\usage{
get_clockrate_table(tree, summary = "mean", file = NULL)
}
\arguments{
  \item{tree}{
a \code{tidytree} object; the output of a call to \code{\link[treeio:read.mrbayes]{treeio::read.mrbayes}}.
}
  \item{summary}{
the name of the rate summary. Should be one of \code{"mean"} or \code{"median"}.
}
  \item{file}{
an optional file path where the resulting table will be stored using \code{\link{write.csv}}.
}
}
\details{
This function converts metadata associated with each tree node (e.g., rate and node age values) into numerical values that can be used for subsequent analyses. It subsequently extracts all clock rate values and outputs them to a table with associated node numbers. 
Note that for evolutionary trees inferred by Mr. Bayes including only fossil taxa, it might be useful to include a dummy extant taxon for correctly outputting node age values (Simões and Pierce, 2021). In such cases, the user may label the latter as “Dummyextant”, which will be automatically pruned by default in this function.  Note that the introduction of a dummy extant taxon is not necessary with BEAST2 if using the treewoffset operator for fossil trees (Barido-Sottani et al., 2019).

}
\value{
A data frame containing a column containing the node identifer (\code{node}) and one column for each clock rate in the tree object containing the rates.
}
\references{
Barido-Sottani, J., Aguirre-Fernández, G., Hopkins, M. J., Stadler, T., Warnock, R. 2019. Ignoring stratigraphic age uncertainty leads to erroneous estimates of species divergence times under the fossilized birth-death process. Proc. R. Soc. Lond., Ser. B: Biol. Sci. 286, 20190685.
Simões, T. R., Pierce, S. E. 2021. Sustained High Rates of Morphological Evolution During the Rise of Tetrapods. Nat. Ecol. Evol. 5, 1403–1414.
}
\author{

Tiago R. Simões & Noah Greifer
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
