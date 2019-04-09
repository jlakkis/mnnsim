\name{docluster}
\alias{docluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Function to Evaluate Quality of Batch Correction via tSNE clustering.
}
\description{
    This is an intermediate function that will perform one of three types of batch correction to simulated data (or it can perform no batch correction). After batch correction, it clusters the data using Rtsne. To evaluate the quality of clustering numerically, a silhoutte coefficient is computed for each cingle cell. The mean of             silhoutte scores is computed for each of the three cell types, as well as for all cells combined and these means are returned. To evaluate the results of clustering visually, colors and shapes are returned for use in the makeTSNE function, which visualizes the results.
}
\usage{
docluster=function(mydata,type="uncorrected",silscores=TRUE,cosnorm=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{mydata}{
    A list object of simulated data items obtained from the "generatedata" function.
    }

    \item{type}{
    The type of batch correction to do before clustering. Options: "uncorrected" "mnn" "limma" "combat"
    }

    \item{silscores}{
    Whether to return silhoutte scores for quantitative evaluation of clustering results, or color and shape objects     for visual evaluation of clustering results.
    }

    \item{cosnorm}{
    Whether to perform cosine normalization for MNN batch correction.
    }
}

\value{
    If silscores=T, then a vector of four numebrs is returned: mean silhoutte score of all cells, mean silhoutte score of cells from type 1, mean silhoutte score of cells from type 2, and mean silhoutte score of cells from type 3.

    If silscores=F, then a list object is returned.
    \item{tsne}{
    The tsne clustering object.}
    \item{batch.id}{
    A vector of batch ids for all cells.}
    \item{clust.cols}{
    A vector of colors corresponding to true cell types.}
}


\examples{
    mydat=generatedata(
        nc=1000,ng=100,
        xm=c(0,5,5),xsigma=c(1,0.1,1),
        ym=c(5,5,0), ysigma=c(1,0.1,1),
        p1=c(0.3,0.5,0.2),p2=c(0.65,0.3,0.05)
    )
    doclust(mydat)
    doclust(mydat,type='mnn')
    doclust(mydat,type='limma')
    doclust(mydat,type='combat')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
