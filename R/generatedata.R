generatedata=function(nc,ng,xm,xsigma,ym,ysigma,p1,p2) {
    checkquality=function(types) {
        return(min(sapply(c(1:3),function(c) sum(types==c))))
    }

    comp1 <- sample(1:3, prob=p1, size=nc, replace=TRUE)

    # Sampling locations for cells in each component.

    samples1 <- cbind(rnorm(n=nc, mean=xm[comp1],sd=xsigma[comp1]),
                      rnorm(n=nc, mean=ym[comp1],sd=ysigma[comp1]))

    # Random projection to D dimensional space, to mimic high-dimensional expression data.
    proj <- matrix(rnorm(ng*nc), nrow=ng, ncol=2)
    A1 <- samples1 %*% t(proj)

    # Add normally distributed noise.
    A1 <- A1 + rnorm(ng*nc)
    rownames(A1) <- paste0("Cell", seq_len(nc), "-1")
    colnames(A1) <- paste0("Gene", seq_len(ng))

    # Setting proportions of each of the three cell types in batch 2.
    comp2 <- sample(1:3, prob=p2, size=nc, replace=TRUE)

    # Sampling locations for cells in each component.
    samples2 <- cbind(rnorm(n=nc, mean=xm[comp2], sd=xsigma[comp2]),
                      rnorm(n=nc, mean=ym[comp2], sd=ysigma[comp2]))

    # Random projection, followed by adding batch effects and random noise.
    A2 <- samples2 %*% t(proj)
    A2 <- A2 + matrix(rep(rnorm(ng), each=nc), ncol=ng) # gene-specific batch effect (genes are columns)
    A2 <- A2 + rnorm(ng*nc) # noise
    rownames(A2) <- paste0("Cell", seq_len(nc), "-2")
    colnames(A2) <- paste0("Gene", seq_len(ng))

    B1 <- t(A1)
    B2 <- t(A2)

    count1=checkquality(comp1)
    count2=checkquality(comp2)

    return(list(B1,B2,comp1,comp2,count1,count2))
}
