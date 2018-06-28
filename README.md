# (expertly named) CPPCxN
This version of PCxN is optimized for speed by adjusting code structure and using C++(Rcpp package). The changes do not alter the idea of the original PCxN method found [here](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006042). 

## Rcpp (C++)
The Rcpp package is widely used to speed-up R scripts by using lower level C++ calculations, majorly decreasing looping overheads (among other things). A large number of R functions are supported by the package (e.g. `colMeans()`) but inevitably a few have to be implemented manually, or worst case imported from R (e.g. `cor.shrink()`).

## R code
The original version is based on two nested loops (through tissues and process pair elements). The latter one is distributed between cores. Turning these 2 loops into Rcpp would be hard as we would need to implement all fucntions included to C++. The main problem would be the `mcapply()` function, which would require us to find a C++ that provides the same functionality.

What we can do is restructure the code in the following ways:
1. Precalculate the Summary matrices (joint and disjoint)
2. Subset the pair elements. That means only calculate interesting relationships (e.g. pathways to drugs)

## Current performance report
We are tackling performance issues from two sides

### R code adjustments/additions
A more detailed report of the runtime comparisons can be found [here](https://docs.google.com/spreadsheets/d/1359vW0Rua5wTmuHGkCloJ8ft-8A-TltPQdUCIocctBE/edit?usp=sharing) In short, the steps implemented:

1. Pre-calculate the matrix that holds the disjoint summaries(estimated represent 40% of the geneset pairs).
2. Pre-calculate the matrix that holds the joint summaries

### New C++ code
Implementing current PCxN functions in Rcpp. Four functions have been translated to C++ but at the moment they don't offer a speed advantage (yet). More effort will go towards implementing them as effieciently as possible. A new function has been created (`precalculate_matrices.cpp`) that pre-calculates both joint and disjoint matrices. To gain speed with this function, C++ and multicores have to be used.

### Memory Issues
The original PCxN code was built and run for 1,330 gene sets/pathways. When we increased the pathways to 5000, memory issues occcured during execution on sharc. Using rmem=48G and mem=48G scripts that handled tissues with more than one GSEs got stopped shortly after the first completed GSE. One thing we can do is adjust the scripts to clean-up after each loop but still the memory needed at any time will definitely be considerably more than before. Only considering a subset of pairs would do wonders for that problem. 

### Comparisons
1. [R/C++](https://docs.google.com/spreadsheets/d/18Z3dXQc22dZ0K_BdF5zopYjVm5YCLSJ_TPR49G47ap0/edit?usp=sharing)
2. [250gene sets,500 gene sets, subpairs](https://docs.google.com/spreadsheets/d/1359vW0Rua5wTmuHGkCloJ8ft-8A-TltPQdUCIocctBE/edit?usp=sharing)

## Suggestions
1. [C++] Try and avoid ALL imports from R. Especially in repeated parts of the code it results in way slower execution. Specifically important case: `cor.shrink()` in *ShrinkCor.cpp* and *ShrinkPCor.cpp* scripts. This is the current reason we can't gain speed from these two scripts.
2. The code is structured to run per tissue. By splitting the batches(each one having a few tissues) we can cut down the total runtime.  

## Notes about the code
1. In Shrink(P)Cor : currently I am importing the cor.shrink function from R and run in C++. Altought the results are identical, this is not speed optimized. The optimal solution would be using a pure C++ function 

