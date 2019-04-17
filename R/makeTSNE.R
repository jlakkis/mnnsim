makeTSNE=function(seed=0,ncells=1000,ngenes=100,xmus=c(0,5,5),xsds=c(1,0.1,1),ymus=c(5,5,0),ysds=c(1,0.1,1),prop1=c(0.3,0.5,0.2),prop2=c(0.65,0.3,0.05)) {
    set.seed(seed)

    if(!length(grep("scran",row.names(installed.packages())))){
        if(!length(grep("BiocManager",row.names(installed.packages())))){
            install.packages("BiocManager")
        }

        BiocManager::install("scran", version = "3.8",ask = F)
    }
    require("scran")

    plotFUN <- function(fname, Y, batch.id, cols, xlab="tSNE 1", ylab="tSNE 2", ...) {
        plot(Y[,1], Y[,2],
             pch=c(16, 2)[batch.id],
             cex=c(2.5, 3.5)[batch.id],
             col=alpha(cols, 0.6),
             xlab=xlab, ylab=ylab, ...)
    }

    mydat=generatedata(ncells,ngenes,xmus,xsds,ymus,ysds,prop1,prop2)

    unc=docluster(mydat,silscores = FALSE)
    mnn=docluster(mydat,type = 'mnn', silscores = FALSE,cosnorm = FALSE)
    lm=docluster(mydat,type = 'limma', silscores = FALSE)
    combat=docluster(mydat,type = 'combat', silscores = FALSE)

    layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))
    plotFUN(paste0("figs/unc.png"), Y=unc[[1]]$Y, batch.id=unc[[2]], cols=unc[[3]], main="Uncorrected")
    plotFUN(paste0("figs/mnn.png"), Y=mnn[[1]]$Y, batch.id=mnn[[2]], col=mnn[[3]], main="MNN")
    plotFUN(paste0("figs/lmfit.png"), Y=lm[[1]]$Y, batch.id=lm[[2]], col=lm[[3]], main="limma")
    plotFUN(paste0("figs/combat.png"), Y=combat[[1]]$Y, batch.id=combat[[2]], col=combat[[3]], main="ComBat")

    plot.new()
    legend('center', legend = c("Cell type 1", "Cell type 2", "Cell type 3", "Batch 1", "Batch 2"),
           col = c("brown1", "gold2", "blue", "black", "black"),
           pch = c(15, 15, 15, 16, 2),
           cex = 1.75,y.intersp = 0.15, bty = "n")
}
