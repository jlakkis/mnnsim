#' @title Perform a preliminary simulation study to estimate number of replicates
#' @description This is a function that, given lists of parameters, performs a preliminary simulation study to estimate number of replicates. It functions similiarly to to the main simulation study function simstudy, except it computes the standard deviations of the simulated quantities of interest to approximate the number of replicates necessary for estimating these parameters to a certain precision.
#' @param tol The maximum standard deviation tolerable for all parameters of interest.
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
#' @param pkeep Whether or not to keep "bad" simulation replicates where a cell type is represented by less than a certain number of cells in a batch: see mycutoff. By default =F.
#' @param mycutoff A number: if the number of cells for any cell type is represented by fewer than cutoff cells in a simulated batch, then this simulation replicate is deemed to be of bad quality. By default=5.
#' @param mycore The number of computing cores to use for parallelizing the simulation.
#' @param dgeneratedata The function to use for generating data. By default this equals generatedata. Highly recommended not to modify this argument.
#' @param ddocluster The function to use for clustering data. By default this equals docluster. Highly recommended not to modify this argument.
#' @export
#' @return A list of simulation replicate numbers calculated such that the standard deviation of each quantity (mean of silhoutte scores for all cells, mean of solhoutte scores for each cell type) for each batch correction method are all no greater than the tol argument for each parameter setting.
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
#' nsims=prelimstudy(
#'     tol=0.01,parameternames=parameternames,
#'     nsims=nsims,seed=seed,
#'     cellcounts=cellcounts,
#'     genecounts=genecounts,
#'     xmeans=xmeans,xsdss=xsdss,ymeans=ymeans,
#'     ysdss=ysdss,propsbatch1=propsbatch1,
#'     propsbatch2=propsbatch2,mycore=1
#' )
#'
#' nsims
#' }
#' @importFrom stats rnorm sd
#'

prelimstudy=function(tol=0.01,parameternames,nsims,seed,cellcounts,genecounts,xmeans,xsdss,ymeans,ysdss,propsbatch1,propsbatch2,pkeep=F,mycutoff=5,mycore=1,dgeneratedata=generatedata,ddocluster=docluster) {
    if(sd(sapply(list(parameternames,nsims,cellcounts,genecounts,xmeans,xsdss,ymeans,ysdss,propsbatch1,propsbatch2),length))!=0) {
        return(warning('All parameter lists should be the same length. Check arguments to make sure this is true.'))
    }

    set.seed(seed)
    seedset=abs(round(rnorm(length(nsims),1000,200)))

    sduc=sdmnn=sdlm=sdcombat=vector(length=length(nsims))

    names(sduc)=names(sdmnn)=names(sdlm)=names(sdcombat)=as.vector(parameternames)
    failedsims=c(1:length(nsims))

    for(i in c(1:length(nsims))) {
        subseed=seedset[i]
        mysim=dosim(nsim=nsims[[i]],
                    ncells=cellcounts[[i]],ngenes=genecounts[[i]],
                    xmus=xmeans[[i]],xsds=xsdss[[i]],
                    ymus=ymeans[[i]],ysds=ysdss[[i]],
                    prop1=propsbatch1[[i]],prop2=propsbatch2[[i]],
                    keep=pkeep,cutoff=mycutoff,ncore=mycore,s.seed=subseed,
                    dgeneratedata = dgeneratedata,ddocluster=ddocluster)
        if(nrow(mysim[[1]])>1) {
            sduc[i]=max(apply(mysim[[1]],MARGIN=2,sd))
            sdmnn[i]=max(apply(mysim[[2]],MARGIN=2,sd))
            sdlm[i]=max(apply(mysim[[3]],MARGIN=2,sd))
            sdcombat[i]=max(apply(mysim[[4]],MARGIN=2,sd))
            failedsims=setdiff(failedsims,i)
        } else {
            warning(paste0('Preliminary Simulation Study Failed for Parameter Set: "', parameternames[[i]] ,'" Due to insufficient valid replicates \n Please Check Parameters \n Returning Supplied replicate number for this Parameter Set'))
        }
    }

    sds=data.frame(sduc,sdmnn,sdlm,sdcombat)
    maxsds=apply(sds,MARGIN=1,max)

    newnsims=ceiling((maxsds/tol)^2)

    if(length(failedsims>0)) {
        newnsims[failedsims]=unlist(nsims)[failedsims]
    }

    return(as.list(newnsims))

}
