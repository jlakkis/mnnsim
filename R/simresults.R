#' @title A function to summarize the results of a simulation
#' @description This function condenses the output of the simstudy function and into simple and interpretable tables from which conclusions can be made about the efficacy of various batch correction methods.
#' @param finishedstudy The output of the simstudy function.
#' @param sds Whether to generate standard deviations for the various summary quantities in the outputted tables.
#' @export
#' @return  A list of five tables
#' \item{All Cells Breakdown}{
#' This table summarizes information from all replicates using all cells.
#' For every parameter combination, this table provides the number of simulation replicates, the mean of mean silhouette coefficients across all replicates for uncorrected data, MNN corrected data, Limma corrected data, and ComBat data, as well as (optionally) the standard deviation of mean silhouette coefficients for each of these batch correction methods.
#' }
#' \item{Cell Type 1 Breakdown}{
#' This table summarizes information from all replicates using cells of type 1.
#' For every parameter combination, this table provides the number of simulation replicates, the mean of mean silhouette coefficients across all replicates for uncorrected data, MNN corrected data, Limma corrected data, and ComBat data, as well as (optionally) the standard deviation of mean silhouette coefficients for each of these batch correction methods.
#' }
#' \item{Cell Type 2 Breakdown}{
#' This table summarizes information from all replicates using cells of type 2.
#' For every parameter combination, this table provides the number of simulation replicates, the mean of mean silhouette coefficients across all replicates for uncorrected data, MNN corrected data, Limma corrected data, and ComBat data, as well as (optionally) the standard deviation of mean silhouette coefficients for each of these batch correction methods.
#' }
#' \item{Cell Type 3 Breakdown}{
#' This table summarizes information from all replicates using cells of type 3.
#' For every parameter combination, this table provides the number of simulation replicates, the mean of mean silhouette coefficients across all replicates for uncorrected data, MNN corrected data, Limma corrected data, and ComBat data, as well as (optionally) the standard deviation of mean silhouette coefficients for each of these batch correction methods.
#' }
#' \item{Replicate Runtime Analysis}{
#' This table summarizes runtime information from all replicates.
#' For every parameter combination, this table provides the mean and standard deviation of runtime per replicate.
#' }
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
#'parameternames=parameternames,
#'nsims=nsims,seed=seed,
#'cellcounts=cellcounts,
#'genecounts=genecounts,
#'xmeans=xmeans,xsdss=xsdss,ymeans=ymeans,
#'ysdss=ysdss,propsbatch1=propsbatch1,
#'propsbatch2=propsbatch2,mycore=1
#')
#'
#'simresults(mystudy)
#'
#'}
#' @importFrom stats sd
#' @importFrom tables tabular Factor

simresults=function(finishedstudy,sds=TRUE) {
    repcount=lapply(finishedstudy,function(x) nrow(x[[1]]))

    numberrows=sum(unlist(repcount))
    failedsims=repcount==0

    if(numberrows==0) {
        warning('Simulation Failed, Change Parameters')
        return(NA)
    } else if (sum(failedsims)>0) {
        warning('Some Simulation Settings Failed: Dropping Failed Simulation Settings')
        parameternames=parameternames[!failedsims]

        for(i in c(1:length(finishedstudy))){
            if(failedsims[i]) {
                finishedstudy[[i]]=NULL
                repcount[[i]]=NULL
            }
        }
    }

    results=resultstype1=resultstype2=resultstype3=as.data.frame(matrix(nrow=numberrows,ncol=5))
    times=as.data.frame(matrix(nrow=length(repcount),ncol=2))

    rownames(times)=as.vector(parameternames)
    colnames(results)=colnames(resultstype1)=colnames(resultstype2)=colnames(resultstype3)=c('Uncorrected', 'MNN', 'Limma', 'ComBat','Parameter.Set')
    colnames(times)=c('Mean Run Time', 'SD Run-Time')

    # typecheck=function(structure) {
    #     if(class(structure()))
    # }

    t=1
    for(i in c(1:length(repcount))) {
        results[c(t:(t+repcount[[i]]-1)),5]=rep(parameternames[[i]],repcount[[i]])
        results[c(t:(t+repcount[[i]]-1)),1]=finishedstudy[[i]][[1]][,1]
        results[c(t:(t+repcount[[i]]-1)),2]=finishedstudy[[i]][[2]][,1]
        results[c(t:(t+repcount[[i]]-1)),3]=finishedstudy[[i]][[3]][,1]
        results[c(t:(t+repcount[[i]]-1)),4]=finishedstudy[[i]][[4]][,1]
        times[i,]=c(mean(finishedstudy[[i]][[5]]),sd(finishedstudy[[i]][[5]]))
        t=t+repcount[[i]]
    }

    t=1
    for(i in c(1:length(repcount))) {
        resultstype1[c(t:(t+repcount[[i]]-1)),5]=rep(parameternames[[i]],repcount[[i]])
        resultstype1[c(t:(t+repcount[[i]]-1)),1]=finishedstudy[[i]][[1]][,2]
        resultstype1[c(t:(t+repcount[[i]]-1)),2]=finishedstudy[[i]][[2]][,2]
        resultstype1[c(t:(t+repcount[[i]]-1)),3]=finishedstudy[[i]][[3]][,2]
        resultstype1[c(t:(t+repcount[[i]]-1)),4]=finishedstudy[[i]][[4]][,2]
        t=t+repcount[[i]]
    }

    t=1
    for(i in c(1:length(repcount))) {
        resultstype2[c(t:(t+repcount[[i]]-1)),5]=rep(parameternames[[i]],repcount[[i]])
        resultstype2[c(t:(t+repcount[[i]]-1)),1]=finishedstudy[[i]][[1]][,3]
        resultstype2[c(t:(t+repcount[[i]]-1)),2]=finishedstudy[[i]][[2]][,3]
        resultstype2[c(t:(t+repcount[[i]]-1)),3]=finishedstudy[[i]][[3]][,3]
        resultstype2[c(t:(t+repcount[[i]]-1)),4]=finishedstudy[[i]][[4]][,3]
        t=t+repcount[[i]]
    }

    t=1
    for(i in c(1:length(repcount))) {
        resultstype3[c(t:(t+repcount[[i]]-1)),5]=rep(parameternames[[i]],repcount[[i]])
        resultstype3[c(t:(t+repcount[[i]]-1)),1]=finishedstudy[[i]][[1]][,4]
        resultstype3[c(t:(t+repcount[[i]]-1)),2]=finishedstudy[[i]][[2]][,4]
        resultstype3[c(t:(t+repcount[[i]]-1)),3]=finishedstudy[[i]][[3]][,4]
        resultstype3[c(t:(t+repcount[[i]]-1)),4]=finishedstudy[[i]][[4]][,4]
        t=t+repcount[[i]]
    }

    results[,5]=as.factor(results[,5])
    resultstype1[,5]=as.factor(resultstype1[,5])
    resultstype2[,5]=as.factor(resultstype2[,5])
    resultstype3[,5]=as.factor(resultstype3[,5])

    if(sds){
        overalltable= tables::tabular( ( Factor(Parameter.Set,'Overall Silhouette') ) ~ (n=1) + Format(digits=2)*(mean + sd)*(Uncorrected + MNN + Limma + ComBat), data=results)
        table1= tables::tabular( ( Factor(Parameter.Set,'Cell Type 1 Silhouette') ) ~ (n=1) + Format(digits=2)*(mean + sd)*(Uncorrected + MNN + Limma + ComBat), data=resultstype1)
        table2= tables::tabular( ( Factor(Parameter.Set,'Cell Type 2 Silhouette') ) ~ (n=1) + Format(digits=2)*(mean + sd)*(Uncorrected + MNN + Limma + ComBat), data=resultstype2)
        table3= tables::tabular( ( Factor(Parameter.Set,'Cell Type 3 Silhouette') ) ~ (n=1) + Format(digits=2)*(mean + sd)*(Uncorrected + MNN + Limma + ComBat), data=resultstype3)
    } else {
        overalltable= tables::tabular( ( Factor(Parameter.Set,'Overall Silhouette') ) ~ (n=1) + Format(digits=2)*mean*(Uncorrected + MNN + Limma + ComBat), data=results)
        table1= tables::tabular( ( Factor(Parameter.Set,'Cell Type 1 Silhouette') ) ~ (n=1) + Format(digits=2)*mean*(Uncorrected + MNN + Limma + ComBat), data=resultstype1)
        table2= tables::tabular( ( Factor(Parameter.Set,'Cell Type 2 Silhouette') ) ~ (n=1) + Format(digits=2)*mean*(Uncorrected + MNN + Limma + ComBat), data=resultstype2)
        table3= tables::tabular( ( Factor(Parameter.Set,'Cell Type 3 Silhouette') ) ~ (n=1) + Format(digits=2)*mean*(Uncorrected + MNN + Limma + ComBat), data=resultstype3)
    }

    myresults=list(overalltable,table1,table2,table3,times)

    names(myresults)=c('All Cells Breakdown', 'Cell Type 1 Breakdown',
                       'Cell Type 2 Breakdown', 'Cell Type 3 Breakdown',
                       'Replicate Runtime Analysis')

    return(myresults)
}
