# mnnsim

# mnnsim 1.0.0

mnnsim (Mutual Nearest Neighbor Simulation) implements a simulation to compare the quality of various batch correction methods for single-cell data.

## Installation

You can install the most recent updates of mnnsim from github with:

```R
if("mnnsim" %in% rownames(installed.packages()) == FALSE) {
    if("devtools" %in% rownames(installed.packages()) == FALSE) 
        install.packages("devtools")
    
    devtools::install_github('jlakkis/mnnsim')
}
```

## Reference
For more information about the package, check out the [vignette](http://rpubs.com/jlakkis/495004). For code assessing the parallel code, check out the [parallel code script](http://rpubs.com/jlakkis/495061).
