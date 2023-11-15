\name{HIVData}
\alias{HIVData}
\alias{OverlapScores}
\alias{LogRankPvals}
\alias{PhenoCodes}
\alias{marker.names}
\docType{data}
\title{HIVData}
\description{
  A dataset of two sets of scores (particularly, correlation with protection against HIV and overlap with the Naive T-cell population) assigned to immunophenotypes measured by flow cytometry. 10 markers were measured: KI-67, CD28, CD45RO, CD8, CD4, CD57, CCR5, CD27, CCR7, and CD127.
}
\usage{data(HIVData)}
\details{
  This dataset consists of 2 numeric
  (\code{LogRankPvals} and \code{OverlapScores}) and 2 character (\code{Pheno.Codes} and \code{marker.names}). 

 \code{Pheno.Codes} vector is described in the RchyOptimyx function's manual.

  \code{LogRankPvals} is a vector of log-rank test P-values to
  determine the correlation between HIV's progression and each of the
  measured immunophenotypes in 466 HIV positive patients (lower values
  represent a stronger correlation). For more details see the Aghaeepour
  et.al. manuscript below.

  Ganesan et. al. define Naive T-cells as CD28+CD45RO-CD57-CCR5-CD27+
  CCR7+ within the CD3+CD14- compartment. The \code{OverlapScores}
  vector has the proportion of Naive T-cells (as defined above) to the
  total number of cells in any given immunophenotype (a higher value
  represents a larger overlap).  

\code{marker.names} is the name of all markers involved in the analysis.


}
%\source{
%}

\author{
  Nima Aghaeepour \email{<naghaeep@gmail.com>} and Adrin Jalali \email{<ajalali@bccrc.ca>}
}

\references{
  Nima Aghaeepour, Pratip K. Chattopadhyay , Anuradha Ganesan, Kieran O'Neill, Habil Zare, Adrin Jalali, Holger H. Hoos, Mario Roederer, and Ryan R. Brinkman. Early Immunologic Correlates of HIV Protection can be Identified from Computational Analysis of Complex Multivariate T-cell Flow Cytometry Assays. Bioinformatics, 2012.

  Anuradha Ganesan, Pratip K Chattopadhyay, Tess M. Brodie, Jing Qin, Wenjuan Gu, John R. Mascola, Nelson L. Michael, Dean A. Follmann, and Mario Roederer. Immunologic and virologic events in early HIV infection predict subsequent rate of progression. Journal of Infectious Diseases, 2010.
}

