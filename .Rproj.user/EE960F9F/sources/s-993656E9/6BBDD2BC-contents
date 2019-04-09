\name{generatedata}
\alias{generatedata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate Dataset for a Simulation Replicate
}
\description{
   Given a set of parameters, this function will simulate data from a single-cell sequencing experiment with batch effects. Called in intermediate steps to generate data for each simulation replicate.
}
\usage{
generatedata(nc,ng,xm,xsigma,ym
            ,ysigma,p1,p2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nc}{
     The number of cells in the sequencing experiment.
}

  \item{ng}{
     The number of genes per cell.
}

  \item{xm}{
      A vector of three numbers, representing the x-coordinates of the three cell types in a two-dimensional space before projection to the space of all genes.
}

  \item{xsigma}{
      The standard deviations of x-coordinates of the three cell types in a two-dimensional space before projection to the space of all genes.
}

  \item{ym}{
      A vector of three numbers, representing the y-coordinates of the three cell types in a two-dimensional space before projection to the space of all genes.
}

  \item{ysigma}{
      The standard deviations of y-coordinates of the three cell types in a two-dimensional space before projection to the space of all genes.
}

  \item{p1}{
      A vector of three numbers, representing the proportions of the three cell types in the first batch. The elements of the vector must sum to 1.
}

  \item{p2}{
      A vector of three numbers, representing the proportions of the three cell types in the second batch. The elements of the vector must sum to 1.
}

}

\value{
  A list of objects detailing the simulated data.
  \item{B1}{Simulated expression data for Batch 1}
  \item{B2}{Simulated expression data for Batch 2}
  \item{comp1}{Simulated cell type identities for Batch 1.}
  \item{comp2}{Simulated cell type identities for Batch 2.}
  \item{count1}{The number of cells for the cell type that least frequently occur in batch 1.}
  \item{count2}{The number of cells for the cell type that least frequently occur in batch 2.}
}

\examples{
    generatedata(
        nc=1000,ng=100,
        xm=c(0,5,5),xsigma=c(1,0.1,1),
        ym=c(5,5,0), ysigma=c(1,0.1,1),
        p1=c(0.3,0.5,0.2),p2=c(0.65,0.3,0.05)
    )
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
