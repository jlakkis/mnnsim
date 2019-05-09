#' @title Perform a simulation study
#' @description This is a function that, given lists of parameters, performs a simulation study. It essentially condenses the simulation study by running dosim multiple times, once per parameter combination, and then saving the output of each dosim in a different element of a list.
#' @param parameternames The names of the parameter combinations supplied. Used only for identification.
#' @param nsims A list giving the number of simulations for each parameter setting.
#' @param seed seed for reproducibility
#' @param cellcounts A list giving the number of cells for each parameter setting.
#' @param genecounts A list giving the number of genes for each parameter setting.
#' @param xmeans A list giving the vectors of mean x-coordinates of the three cell types for each parameter setting.
#' @param xsdss A list giving the vectors of the standard deviations of x-coordinates of the three cell types for each parameter setting
#' @param ymeans A list giving the vectors of mean y-coordinates of the three cell types for each parameter setting.
#' @param ysdss A list giving the vectors of the standard deviations of y-coordinates of the three cell types for each parameter setting
#' @param propsbatch1 A list giving the vectors of proportions of the three cell types in the first batch for each parameter setting.
#' @param propsbatch2 A list giving the vectors of proportions of the three cell types in the second batch for each parameter setting.
#' @param mykeep Whether or not to keep "bad" simulation replicates where a cell type is represented by less than a certain number of cells in a batch: see mycutoff. By default =F.
#' @param mycutoff A number: if the number of cells for any cell type is represented by fewer than cutoff cells in a simulated batch, then this simulation replicate is deemed to be of bad quality. By default=5.
#' @param mycore The number of computing cores to use for parallelizing the simulation.
#' @param dgeneratedata The function to use for generating data. By default this equals generatedata. Highly recommended not to modify this argument.
#' @param ddocluster The function to use for clustering data. By default this equals docluster. Highly recommended not to modify this argument.
#' @export
#' @return A list of simulation results. Each element of the list contains simulations results (output of dosim) for one of the supplied combinations of parameters.
#'
#' As mentioned in the documentation of dosim, the simulation results of a single simulation consist of the following list of simulation items:
#'
#'  \item{Uncorrected Silhoutte Scores}{
#'  A matrix of clustering results from simulation replicates for the uncorrected batch correction setting. Each row is a vector of four mean silhoutte scores from a single replicate. This vector consists of: mean silhoutte score of all cells, mean silhoutte score of cells from type 1, mean silhoutte score of cells from type 2, and mean silhoutte score of cells from type 3.
#'  }
#'  \item{MNN Silhoutte Scores}{
#'   A matrix of clustering results from simulation replicates for the MNN batch correction setting. Each row is a vector of four mean silhoutte scores from a single replicate. This vector consists of: mean silhoutte score of all cells, mean silhoutte score of cells from type 1, mean silhoutte score of cells from type 2, and mean silhoutte score of cells from type 3.
#'   }
#'
#'   \item{limma Silhoutte Scores}{
#'   A matrix of clustering results from simulation replicates for the limma batch correction setting. Each row is a vector of four mean silhoutte scores from a single replicate. This vector consists of: mean silhoutte score of all cells, mean silhoutte score of cells from type 1, mean silhoutte score of cells from type 2, and mean silhoutte score of cells from type 3.
#'   }
#'   \item{combat Silhoutte Scores}{
#'   A matrix of clustering results from simulation replicates for the combat batch correction setting. Each row is a vector of four mean silhoutte scores from a single replicate. This vector consists of: mean silhoutte score of all cells, mean silhoutte score of cells from type 1, mean silhoutte score of cells from type 2, and mean silhoutte score of cells from type 3.
#'   }
#'
#'  \item{Replicate Runtimes}{
#'  A vector of real-times (in seconds) necessary to complete each replicate.
#'  }
#' @examples \dontrun{
#' parameternames=list('Original', 'Smaller Differences', 'More Genes')
#' nsims=list(50,50,50)
#' seed=0
#' cellcounts=list(500,500,500)
#' genecounts=list(100,100,5000)
#' xmeans=list(c(0,5,5),c(0,2,2),c(0,5,5))
#' xsdss=list(c(1,0.1,1),c(1,0.1,1),c(1,0.1,1))
#' ymeans=list(c(5,5,0),c(2,2,0),c(5,5,0))
#' ysdss=list(c(1,0.1,1),c(1,0.1,1),c(1,0.1,1))
#' propsbatch1=list(c(0.3,0.5,0.2),c(0.3,0.5,0.2),c(0.3,0.5,0.2))
#' propsbatch2=list(c(0.65,0.3,0.05),c(0.65,0.3,0.05),c(0.65,0.3,0.05))
#'
#'simstudy(
#'    parameternames=parameternames,
#'    nsims=nsims,seed=seed,
#'    cellcounts=cellcounts,
#'    genecounts=genecounts,
#'    xmeans=xmeans,xsdss=xsdss,ymeans=ymeans,
#'    ysdss=ysdss,propsbatch1=propsbatch1,
#'    propsbatch2=propsbatch2,mycore=1
#')
#'}
#' @importFrom stats rnorm sd

simstudy=function(parameternames,nsims,seed,cellcounts,genecounts,xmeans,xsdss,ymeans,ysdss,propsbatch1,propsbatch2,mykeep=F,mycutoff=5,mycore=1,dgeneratedata=generatedata,ddocluster=docluster) {
    if(sd(sapply(list(parameternames,nsims,cellcounts,genecounts,xmeans,xsdss,ymeans,ysdss,propsbatch1,propsbatch2),length))!=0) {
        return(warning('All parameter lists should be the same length. Check arguments to make sure this is true.'))
    }

    set.seed(seed)
    seedset=abs(round(rnorm(length(nsims),1000,200)))

    results=as.list(vector(length=length(nsims)))
    names(results)=as.vector(parameternames)

    for(i in c(1:length(nsims))) {
        subseed=seedset[i]
        results[[i]]=dosim(nsim=nsims[[i]],
                           ncells=cellcounts[[i]],ngenes=genecounts[[i]],
                           xmus=xmeans[[i]],xsds=xsdss[[i]],
                           ymus=ymeans[[i]],ysds=ysdss[[i]],
                           prop1=propsbatch1[[i]],prop2=propsbatch2[[i]],
                           keep=mykeep,cutoff=mycutoff,ncore=mycore,s.seed=subseed,
                           dgeneratedata = dgeneratedata,ddocluster=ddocluster)
    }

    return(results)
}
