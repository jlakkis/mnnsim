% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dosim.R
\name{dosim}
\alias{dosim}
\title{generatedata}
\usage{
dosim(nsim, ncells, ngenes, xmus, xsds, ymus, ysds, prop1, prop2,
  keep = F, cutoff = 5, ncore, s.seed, dgeneratedata = generatedata,
  ddocluster = docluster)
}
\arguments{
\item{nsim}{number of simulations}

\item{ncells}{number of single-cells}

\item{ngenes}{number of genes}

\item{xmus}{low-d x-means}

\item{xsds}{low-d x-sds}

\item{ymus}{low-d y-means}

\item{ysds}{low-d y-sds}

\item{prop1}{batch 1 cell type proportions}

\item{prop2}{batch 2 cell type proportions}

\item{keep}{whether to keep "low-quality" replicates}

\item{cutoff}{criterion for determining whether a replicate is low-quality}

\item{ncore}{number of cores for parallel computing}

\item{s.seed}{rng seed}

\item{dgeneratedata}{technical argument}

\item{ddocluster}{technical argument}
}
\description{
generatedata
}
