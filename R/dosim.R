dosim=function(nsim,ncells,ngenes,xmus,xsds,ymus,ysds,prop1,prop2,keep=F,cutoff=5,ncore,s.seed,dgeneratedata=generatedata,ddocluster=docluster) {

    message('Starting a Simulation \n')

    lparallizer=function(iter,pgeneratedata,pdocluster,pncells,pngenes,pxmus,pxsds,pymus,pysds,pprop1,pprop2,pcutoff,pkeep) {

        tictoc::tic()
        mydat=pgeneratedata(pncells,pngenes,pxmus,pxsds,pymus,pysds,pprop1,pprop2)
        if(base::min(c(mydat[[5]],mydat[[6]]))<pcutoff) {
            myoutput=base::list(rep(NA,4),rep(NA,4),rep(NA,4),rep(NA,4))
            if(pkeep) {
                myoutput=base::list(pdocluster(mydat),pdocluster(mydat,type='mnn'),pdocluster(mydat,type='limma'),pdocluster(mydat,type='combat'))
            }
            myoutput[[5]]=T
        } else {
            myoutput=base::list(pdocluster(mydat),pdocluster(mydat,type='mnn'),pdocluster(mydat,type='limma'),pdocluster(mydat,type='combat'))
            myoutput[[5]]=F
        }

        info=tictoc::toc()

        myoutput[[6]]=info$toc-info$tic

        base::return(myoutput)
    }

    if(ncore>1) {
        RNGkind("L'Ecuyer-CMRG") # Set RNG
        cl = parallel::makeCluster(ncore)
        parallel::clusterSetRNGStream(cl, iseed=s.seed)
        output = parallel::parLapply(cl, as.list(c(1:nsim)), lparallizer,pgeneratedata=dgeneratedata,pdocluster=ddocluster,pncells=ncells,pngenes=ngenes,pxmus=xmus,pxsds=xsds,pymus=ymus,pysds=ysds,pprop1=prop1,pprop2=prop2,pcutoff=cutoff,pkeep=keep)
        stopCluster(cl)
    } else {
        set.seed(s.seed)
        lapply(cl, as.list(c(1:nsim)), lparallizer,pgeneratedata=dgeneratedata,pdocluster=ddocluster,pncells=ncells,pngenes=ngenes,pxmus=xmus,pxsds=xsds,pymus=ymus,pysds=ysds,pprop1=prop1,pprop2=prop2,pcutoff=cutoff,pkeep=keep)
    }

    coefficientsuc=t(sapply(output,function(r) r[[1]]))
    coefficientsmnn=t(sapply(output,function(r) r[[2]]))
    coefficientslm=t(sapply(output,function(r) r[[3]]))
    coefficientscombat=t(sapply(output,function(r) r[[4]]))
    warntrack=sapply(output,function(r) r[[5]])
    times=sapply(output,function(r) r[[6]])

    names=c('overall', 'celltype1', 'celltype2', 'celltype3')
    colnames(coefficientsuc) = colnames(coefficientsmnn) = colnames(coefficientslm) = colnames(coefficientscombat) = names

    times=times[!is.na(coefficientsuc[,1])]
    coefficientsuc=coefficientsuc[!is.na(coefficientsuc[,1]),]
    coefficientsmnn=coefficientsmnn[!is.na(coefficientsmnn[,1]),]
    coefficientslm=coefficientslm[!is.na(coefficientslm[,1]),]
    coefficientscombat=coefficientscombat[!is.na(coefficientscombat[,1]),]

    numwarn=sum(warntrack==T)

    if(numwarn>0) {
        if(keep == T) {
            warning(paste0('Warning: A cell type is underrepresented in a batch for '), numwarn, ' Simulation Replicates \n Check Parameters')
        } else{
            warning(paste0('Warning: A cell type is underrepresented in a batch for '), numwarn, ' Simulation Replicates \n Dropping These Replicates')
        }
    }

    message('\n')
    myoutput=list(coefficientsuc,coefficientsmnn,coefficientslm,coefficientscombat,times)
    names(myoutput)=c('Uncorrected Silhoutte Scores', 'MNN Silhoutte Scores', 'limma Silhoutte Scores', 'combat Silhoutte Scores', 'Replicate Runtimes')

    return(myoutput)

}
