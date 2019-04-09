\name{makeTSNE}
\alias{makeTSNE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a tSNE Plot to visually evaluate the quality of batch correction methods
}
\description{
Given a set of parameters, this function will simulate data and then output four tSNE plots: one corresponding to clustering after no batch correction, and another three, each corresponding to clustering after one of the three batch correction methods "limma," "combat," and "mnn."
}
\usage{
makeTSNE(
    seed,ncells,ngenes,
    xmus,xsds,
    ymus,ysds,
    prop1,prop2
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seed}{
    seed for reproducibility
}
  \item{ncells}{
     The number of cells in the sequencing experiment.
}

  \item{ngenes}{
     The number of genes per cell.
}

  \item{xmus}{
      A vector of three numbers, representing the x-coordinates of the three cell types in a two-dimensional space before projection to the space of all genes.
}

  \item{xsds}{
      The standard deviations of x-coordinates of the three cell types in a two-dimensional space before projection to the space of all genes.
}

  \item{ymus}{
      A vector of three numbers, representing the y-coordinates of the three cell types in a two-dimensional space before projection to the space of all genes.
}

  \item{ysds}{
      The standard deviations of y-coordinates of the three cell types in a two-dimensional space before projection to the space of all genes.
}

  \item{prop1}{
      A vector of three numbers, representing the proportions of the three cell types in the first batch. The elements of the vector must sum to 1.
}

  \item{prop2}{
      A vector of three numbers, representing the proportions of the three cell types in the second batch. The elements of the vector must sum to 1.
}

}
\value{
Returns four tSNE plots (on the same figure).
}
\examples{
makeTSNE(
    seed=0,ncells=1000,ngenes=100,
    xmus=c(0,5,5),xsds=c(1,0.1,1),
    ymus=c(5,5,0),ysds=c(1,0.1,1),
    prop1=c(0.3,0.5,0.2),prop2=c(0.65,0.3,0.05)
)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
