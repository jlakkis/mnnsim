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
For a detailed walkthrough of this package, please check out this [vignette](http://rpubs.com/jlakkis/495004). For a quantitative assessment of the benefits of parallelizing this simulation, check out the [parallel code script](http://rpubs.com/jlakkis/495061).
