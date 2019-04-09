\name{prelimstudy}
\alias{prelimstudy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Perform a preliminary simulation study to estimate number of replicates
}
\description{
    This is a function that, given lists of parameters, performs a preliminary simulation study to estimate number of replicates. It functions similiarly to to the main simulation study function simstudy, but computes the standard deviations of the simulated quantities of interest to estimate the necessary number of replicates necessary for estimating these parameters to a certain precision.
}

\usage{
prelimstudy(
    tol=0.01,parameternames,nsims,
    seed,cellcounts,genecounts,
    xmeans,xsdss,ymeans,ysdss,
    propsbatch1,propsbatch2,
    pkeep=F,cutoff=5,mycore=1,
    dgeneratedata=generatedata,
    ddocluster=docluster)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tol}{
    The maximum standard deviation tolerable for all parameters of interest.
    }

  \item{parameternames}{
    The names of the parameter combinations supplied. Used only for identification.
    }

  \item{nsims}{
    A list giving the number of simulations for each parameter setting.
}

  \item{seed}{
    seed for reproducibility
}

  \item{cellcounts}{
     A list giving the number of cells for each parameter setting.
}

  \item{genecounts}{
     A list giving the number of genes for each parameter setting.
}

  \item{xmeans}{
      A list giving the vectors of mean x-coordinates of the three cell types for each parameter setting.
}

  \item{xsdss}{
      A list giving the vectors of the standard deviations of x-coordinates of the three cell types for each parameter setting.
}

  \item{ymeans}{
      A list giving the vectors of mean y-coordinates of the three cell types for each parameter setting.
}

  \item{ysdss}{
      A list giving the vectors of the standard deviations of y-coordinates of the three cell types for each parameter setting.
}

  \item{propsbatch1}{
        A list giving the vectors of proportions of the three cell types in the first batch for each parameter setting.
}

  \item{propsbatch2}{
        A list giving the vectors of proportions of the three cell types in the second batch for each parameter setting.
}

  \item{pkeep}{
    Whether or not to keep "bad" simulation replicates where a cell type is represented by less than a certain number of cells (see cutoff) in a batch. By default =F.
}
  \item{cutoff}{
     A number: if the number of cells for any cell type is represented by fewer than cutoff cells in a simulated batch, then this simulation replicate is deemed to be of bad quality. By default=5.
}

  \item{mycore}{
     The number of computing cores to use for parallelizing the simulation.
}

  \item{degeneratedata}{
      The function to use for generating data. By default this equals generatedata. Highly recommended not to touch this argument.
}

  \item{ddocluster}{
      The function to use for clustering data. By default this equals docluster. Highly recommended not to touch this argument.
}

  \item{degeneratedata}{
      The function to use for generating data. By default this equals generatedata. Highly recommended not to touch this argument.
}

  \item{ddocluster}{
      The function to use for clustering data. By default this equals docluster. Highly recommended not to touch this argument.
}

}
\value{
    A list of simulation replicate numbers calculated such that the standard deviation of each quantity (mean of silhoutte scores for all cells, mean of solhoutte scores for each cell type) for each batch correction method are all no greater than the tol argument for each parameter setting.
}

\examples{
parameternames=list('Set 1', 'Set 2')
nsims=list(10,10)
seed=0
cellcounts=list(500,500)
genecounts=list(100,100)
xmeans=list(c(0,5,5),c(0,2,2))
xsdss=list(c(1,0.1,1),c(1,0.1,1))
ymeans=list(c(5,5,0),c(2,2,0))
ysdss=list(c(1,0.1,1),c(1,0.1,1))
propsbatch1=list(c(0.3,0.5,0.2),c(0.3,0.5,0.2))
propsbatch2=list(c(0.65,0.3,0.05),c(0.65,0.3,0.05))

nsims=prelimstudy(
    tol=0.01,parameternames=parameternames,
    nsims=nsims,seed=seed,
    cellcounts=cellcounts,
    genecounts=genecounts,
    xmeans=xmeans,xsdss=xsdss,ymeans=ymeans,
    ysdss=ysdss,propsbatch1=propsbatch1,
    propsbatch2=propsbatch2,mycore=1
)

nsims

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
