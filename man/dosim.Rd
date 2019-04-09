\name{dosim}
\alias{dosim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    A function to perform a simulation for a set of parameters.
}
\description{
    Given a set of parameters, this function will generate a specified number simulation replicates. For each replicate, the docluster function is called 4 times: once with no batch corrections, once with mnn batch correction, once with limma batch correction, and once with combat batch corrction: outputting four mean silhoutte scores per batch correction type and replicate (see docluster). In addition, the runtime for each replicate is saved. Called in intermediate steps by other functions.
}
\usage{

dosim(
    nsim,ncells,ngenes,xmus,xsds,ymus,ysds,
    prop1,prop2,keep=F,cutoff=T,ncore,s.seed,
    dgeneratedata=generatedata,ddocluster=docluster
)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsim}{
    The number of simulation replicates
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

  \item{keep}{
    Whether or not to keep "bad" simulation replicates where a cell type is represented by less than a certain number of cells (see cutoff) in a batch. By default =F.
}
  \item{cutoff}{
     A number: if the number of cells for any cell type is represented by fewer than cutoff cells in a simulated batch, then this simulation replicate is deemed to be of bad quality. By default=5.
}

  \item{ncore}{
     The number of computing cores to use for parallelizing the simulation.
}

  \item{degeneratedata}{
      The function to use for generating data. By default this equals generatedata. Highly recommended not to touch this argument.
}

  \item{ddocluster}{
      The function to use for clustering data. By default this equals docluster. Highly recommended not to touch this argument.
}

}
\value{
    A list of simulation items
  \item{Uncorrected Silhoutte Scores}{
      A matrix of clustering results from simulation replicates for the uncorrected batch correction setting. Each row is a vector of four mean silhoutte scores from a single replicate. This vector consists of: mean silhoutte score of all cells, mean silhoutte score of cells from type 1, mean silhoutte score of cells from type 2, and mean silhoutte score of cells from type 3.
}

  \item{MNN Silhoutte Scores}{
      A matrix of clustering results from simulation replicates for the MNN batch correction setting. Each row is a vector of four mean silhoutte scores from a single replicate. This vector consists of: mean silhoutte score of all cells, mean silhoutte score of cells from type 1, mean silhoutte score of cells from type 2, and mean silhoutte score of cells from type 3.
}

  \item{limma Silhoutte Scores}{
      A matrix of clustering results from simulation replicates for the limma batch correction setting. Each row is a vector of four mean silhoutte scores from a single replicate. This vector consists of: mean silhoutte score of all cells, mean silhoutte score of cells from type 1, mean silhoutte score of cells from type 2, and mean silhoutte score of cells from type 3.
}

  \item{combat Silhoutte Scores}{
      A matrix of clustering results from simulation replicates for the combat batch correction setting. Each row is a vector of four mean silhoutte scores from a single replicate. This vector consists of: mean silhoutte score of all cells, mean silhoutte score of cells from type 1, mean silhoutte score of cells from type 2, and mean silhoutte score of cells from type 3.
}

 \item{Replicate Runtimes}{
    A vector of real-times (in seconds) necessary to complete each replicate.
 }

}
\examples{

dosim(
    nsim=30,ncells=500,ngenes=100,
    xmus=c(0,5,5),xsds=c(1,0.1,1),
    ymus=c(5,5,0),ysds=c(1,0.1,1),
    prop1=c(0.3,0.5,0.2),prop2=c(0.65,0.3,0.05),
    keep=F,cutoff=5,ncore=1,
    s.seed=0, dgeneratedata=generatedata,
    ddocluster=docluster
)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
