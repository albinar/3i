\name{RchyOptimyx}
\docType{methods}
\alias{RchyOptimyx}
\alias{RchyOptimyx-methods}
\alias{RchyOptimyx,character,numeric,character,numeric,logical-method}
\title{
  Gating Hierarchy Optimization for Flow Cytometry
}
\description{
  RchyOptimyx is a computational tool which determines the minimal sets of
  markers that can identify a target population within a given purity
  standard. This can enable panel redesign using smaller subsets of
  markers, with consequent savings on reagents, potential for adaptation to older instruments, or the addition of other markers of interest to the panel design.
}
\usage{
\S4method{RchyOptimyx}{character,numeric,character,numeric,logical}(pheno.codes, phenotypeScores, startPhenotype, pathCount=1, trimPaths=FALSE, trim.tolerance=0, trim.level=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pheno.codes}{
    A vector of strings of length N (the number of markers) for each phenotype measured. For every phenotype, the character corresponding to a given marker can be 0, 1, 2, etc for
    neutral, negative, positive, bright, etc. See
    the provided vigenette (inst/doc/RchyOptimyx.pdf) for more details
    and examples.
  }
  \item{phenotypeScores}{
    A vector of the scores assigned to every phenotype. The optimization
    algorithm will try to maximize this score. See the provided vigenette (inst/doc/RchyOptimyx.pdf) for more details
    and examples.
  }
  \item{startPhenotype}{
    The terminal phenotype which includes every marker that must be
    considered in the analysis. This variable is a vector of the length
    of the number of markers. Every element of the vector can be 0, 1,
    or 2 for negative, neutral, or positive respectively (see details).
  }
  \item{pathCount}{
    The number of paths that must be generated.
  }
  \item{trimPaths}{
    A boolean indicating wether the paths should be trimmed or not. If
    set to TRUE, the paths will be terminated as soon as adding
    \code{trim.tolerance} extra number of edges to
    the path doesn't result in an improvement in the score.
  }
  \item{trim.tolerance}{
    An integer indicating the number of levels before the path is trimed if the score is decreasing by adding extra edges.
  }
  \item{trim.level}{
    An integer indicating from which level the tree is going to be
    trimed. This parameter has no effect if set to zero.
  }
}
\details{
  Every marker can be positive, negative, or neutral. Neutral markers
  are not included in the measurement of the respective phenotype.
  
  A path is trimed in one of the following conditions:
  1. If the \code{trimPaths} parameter is \code{TRUE}, the path is trimed as soon as it decreases for equal or more levels that \code{trim.trolerance} parameter. If the path is decreasing while reaching the \code{trim.level}, the path is trimed from the point that it started to decrease.
  2. If the path reaches the \code{trim.level}.
}
\value{
  \item{OptH}{An OptimizedHierarchy object}
}
\author{
  Nima Aghaeepour \email{<naghaeep@gmail.com>} and Adrin Jalali \email{<ajalali@bccrc.ca>}
}

\references{
  D. Eppstein, Finding the k shortest paths, SIAM J. Comput., 1998 - Citeseer.

  Nima Aghaeepour, Pratip K. Chattopadhyay, Anuradha Ganesan, Kieran
  O'Neill, Habil Zare, Holger H. Hoos, Mario Roederer, and Ryan
  R. Brinkman. Early Immunologic Correlates of HIV Protection can be
  Identified from Computational Analysis of Complex Multivariate T-cell Flow Cytometry Assays. Bioinformatics, Feb. 2012.

}


\examples{
data(HIVData)

res<-RchyOptimyx(PhenoCodes, -log10(LogRankPvals), '2000102001', 5,FALSE) 
plot(res,phenotypeCodes=PhenoCodes,phenotypeScores=-log10(LogRankPvals),marker.names=marker.names)
}

\keyword{ FlowCytData }
