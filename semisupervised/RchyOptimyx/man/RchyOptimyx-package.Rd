\name{RchyOptimyx-package}
\alias{RchyOptimyx-package}
\docType{package}
\title{
Gating Hierarchy Optimization for Flow Cytometry
}
\description{
RchyOptimyx is a computational tool which determines the minimal sets of
markers that can identify a target population within a given purity
standard. This can enable panel redesign using smaller subsets of
markers, with consequent savings on reagents, potential for adaptation to older instruments, or the addition of other markers of interest to the panel design.
}
\details{
\tabular{ll}{
Package: \tab RchyOptimyx\cr
Type: \tab Package\cr
Version: \tab 0.99.1\cr
Date: \tab 2011-10-12\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
}
See the \code{RchyOptimyx} function.
}

\keyword{ package }
\seealso{
  RchyOptimyx
}

\references{
D. Eppstein, Finding the k shortest paths, SIAM J. Comput., 1998 - Citeseer.
}



\examples{
data(HIVData)

res<-RchyOptimyx(PhenoCodes, -log10(LogRankPvals), '2000102001', 5,FALSE) 
plot(res,phenotypeCodes=PhenoCodes,phenotypeScores=-log10(LogRankPvals),marker.names=marker.names)
}

\author{
  Nima Aghaeepour \email{<naghaeep@gmail.com>} and Adrin Jalali \email{<ajalali@bccrc.ca>}
}

\keyword{ FlowCytData }
