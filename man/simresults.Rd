\name{simresults}
\alias{simresults}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    A function to summarize the results of a simulation
}
\description{
    This function takes the complicated output of the simstudy function and condenses it into simple and interpretable tables from which conclusions can be made about the efficacy of various batch correction methods for various parameter settings.
}
\usage{
simresults(finishedstudy,sds=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{finishedstudy}{
    The output of the simstudy function.
    }
    \item{sds}{
    Whether to generate standard deviations for the various summary quantities in the       outputted tables.
    }
}

\value{
A list of five tables
  \item{All Cells Breakdown}{
  This table summarizes information from all replicates using all cells.
    For every parameter combination, this table provides the number of simulation replicates, the mean of mean silhouette coefficients across all replicates for uncorrected data, MNN corrected data, Limma corrected data, and ComBat data, as well as (optionally) the standard deviation of mean silhouette coefficients for each of these batch correction methods.
  }
  \item{Cell Type 1 Breakdown}{
  This table summarizes information from all replicates using cells of type 1.
    For every parameter combination, this table provides the number of simulation replicates, the mean of mean silhouette coefficients across all replicates for uncorrected data, MNN corrected data, Limma corrected data, and ComBat data, as well as (optionally) the standard deviation of mean silhouette coefficients for each of these batch correction methods.
  }
    \item{Cell Type 2 Breakdown}{
  This table summarizes information from all replicates using cells of type 2.
    For every parameter combination, this table provides the number of simulation replicates, the mean of mean silhouette coefficients across all replicates for uncorrected data, MNN corrected data, Limma corrected data, and ComBat data, as well as (optionally) the standard deviation of mean silhouette coefficients for each of these batch correction methods.
  }
    \item{Cell Type 3 Breakdown}{
  This table summarizes information from all replicates using cells of type 3.
    For every parameter combination, this table provides the number of simulation replicates, the mean of mean silhouette coefficients across all replicates for uncorrected data, MNN corrected data, Limma corrected data, and ComBat data, as well as (optionally) the standard deviation of mean silhouette coefficients for each of these batch correction methods.
  }

      \item{Replicate Runtime Analysis}{
  This table summarizes runtime information from all replicates.
    For every parameter combination, this table provides the mean and standard deviation of runtime per replicate.
  }

  }
\examples{
    parameternames=list('Set 1', 'Set 2')
    nsims=list(15,15)
    seed=0
    cellcounts=list(500,500)
    genecounts=list(100,100)
    xmeans=list(c(0,5,5),c(0,2,2))
    xsdss=list(c(1,0.1,1),c(1,0.1,1))
    ymeans=list(c(5,5,0),c(2,2,0))
    ysdss=list(c(1,0.1,1),c(1,0.1,1))
    propsbatch1=list(c(0.3,0.5,0.2),c(0.3,0.5,0.2))
    propsbatch2=list(c(0.65,0.3,0.05),c(0.65,0.3,0.05))

mystudy=simstudy(
    parameternames=parameternames,
    nsims=nsims,seed=seed,
    cellcounts=cellcounts,
    genecounts=genecounts,
    xmeans=xmeans,xsdss=xsdss,ymeans=ymeans,
    ysdss=ysdss,propsbatch1=propsbatch1,
    propsbatch2=propsbatch2,mycore=1
)

simresults(mystudy)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
