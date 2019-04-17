prelimstudy=function(tol=0.01,parameternames,nsims,seed,cellcounts,genecounts,xmeans,xsdss,ymeans,ysdss,propsbatch1,propsbatch2,pkeep=F,cutoff=5,mycore=1,dgeneratedata=generatedata,ddocluster=docluster) {
    if(sd(sapply(list(parameternames,nsims,cellcounts,genecounts,xmeans,xsdss,ymeans,ysdss,propsbatch1,propsbatch2),length))!=0) {
        return(warning('All parameter lists should be the same length. Check arguments to make sure this is true.'))
    }

    if(!length(grep("scran",row.names(installed.packages())))){
        if(!length(grep("BiocManager",row.names(installed.packages())))){
            install.packages("BiocManager")
        }

        BiocManager::install("scran", version = "3.8",ask = F)
    }
    require("scran")

    set.seed(seed)

    sduc=sdmnn=sdlm=sdcombat=vector(length=length(nsims))

    names(sduc)=names(sdmnn)=names(sdlm)=names(sdcombat)=as.vector(parameternames)
    failedsims=c(1:length(nsims))

    for(i in c(1:length(nsims))) {
        subseed=abs(round(rnorm(1,1000,200)))
        mysim=dosim(nsim=nsims[[i]],
                    ncells=cellcounts[[i]],ngenes=genecounts[[i]],
                    xmus=xmeans[[i]],xsds=xsdss[[i]],
                    ymus=ymeans[[i]],ysds=ysdss[[i]],
                    prop1=propsbatch1[[i]],prop2=propsbatch2[[i]],
                    keep=pkeep,cut=cutoff,ncore=mycore,s.seed=subseed,
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
