#' @title docluster
#' @importFrom stats dist
#' @importFrom batchelor mnnCorrect
#' @importFrom SummarizedExperiment assays
#' @importFrom limma removeBatchEffect
#' @importFrom sva ComBat
#' @importFrom Rtsne Rtsne
#' @importFrom cluster silhouette
#' @param mydata output from generatedata function
#' @param type batch correction method to use before clustering
#' @param silscores whether to produce silhoutte scores or not
#' @param cosnorm whether to perform cosine normalization

docluster=function(mydata,type="uncorrected",silscores=TRUE,cosnorm=TRUE) {
    if(!(type %in% c("uncorrected","mnn","limma","combat"))){
        return(warning('Not a Valid Method: Please set argument "type" to one of: "uncorrected","mnn","limma","combat"'))
    }

    quiet <- function(x) {
        sink(tempfile())
        on.exit(sink())
        invisible(force(x))
    }

    raw.all <- cbind(mydata[[1]], mydata[[2]])
    batch.id <- rep(1:2, c(ncol(mydata[[1]]), ncol(mydata[[2]])))

    if(type=='uncorrected') {
        all.dists <- as.matrix(dist(t(raw.all)))
    } else if (type =='mnn') {
        Xmnn <- batchelor::mnnCorrect(raw.all, batch = batch.id, k=20, sigma=1, cos.norm.in=cosnorm, cos.norm.out=cosnorm, var.adj=TRUE)
        Xmnn = SummarizedExperiment::assays(Xmnn)
        Xmnn = Xmnn$'corrected'
        all.dists <- as.matrix(dist(t(Xmnn)))
    } else if (type =='limma') {
        Xlm <- limma::removeBatchEffect(raw.all, factor(batch.id))
        all.dists <- as.matrix(dist(t(Xlm)))
    } else {
        cleandat <- quiet(suppressMessages( sva::ComBat(raw.all, factor(batch.id), mod=NULL, prior.plots = FALSE) ))
        all.dists <- as.matrix(dist(t(cleandat)))
    }

    tsne = tryCatch ({
        Rtsne::Rtsne(all.dists, is_distance=TRUE)
    },
    error = function (cond){
        Rtsne::Rtsne(all.dists, is_distance=TRUE,perplexity = 10)
    })

    if(!silscores) {
        ref.cols <- c("blue", "brown1", "gold2")
        clust1 <- ref.cols[mydata[[3]]]
        clust2 = ref.cols[mydata[[4]]]
        clust.cols <- c(clust1, clust2)
        return(list(tsne,batch.id,clust.cols))
    }

    dd <- as.matrix(dist(tsne$Y))

    ct <- c(mydata[[3]],mydata[[4]])

    sil <- (cluster::silhouette(as.numeric(ct), dd))

    coefficients=c(mean(sil[,3]),
                   mean(sil[sil[,1]==1,3]),
                   mean(sil[sil[,1]==2,3]),
                   mean(sil[sil[,1]==3,3]))

    return(coefficients)
}
