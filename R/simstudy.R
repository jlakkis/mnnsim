simstudy=function(parameternames,nsims,seed,cellcounts,genecounts,xmeans,xsdss,ymeans,ysdss,propsbatch1,propsbatch2,mykeep=F,mycut=5,mycore=1,dgeneratedata=generatedata,ddocluster=docluster) {
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

    results=as.list(vector(length=length(nsims)))
    names(results)=as.vector(parameternames)

    for(i in c(1:length(nsims))) {
        subseed=abs(round(rnorm(1,1000,200)))
        results[[i]]=dosim(nsim=nsims[[i]],
                           ncells=cellcounts[[i]],ngenes=genecounts[[i]],
                           xmus=xmeans[[i]],xsds=xsdss[[i]],
                           ymus=ymeans[[i]],ysds=ysdss[[i]],
                           prop1=propsbatch1[[i]],prop2=propsbatch2[[i]],
                           keep=mykeep,cut=mycut,ncore=mycore,s.seed=subseed,
                           dgeneratedata = dgeneratedata,ddocluster=ddocluster)
    }

    return(results)
}
