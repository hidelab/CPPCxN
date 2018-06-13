# (expertly named) CPPCxN
This version of PCxN is optimized for speed by adjusting code structure and using C++(Rcpp package). The changes do not alter the idea or results of the original PCxN method found [here](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006042). 

## Rcpp (C++)
The Rcpp package is widely used to speed-up R scripts by using lower level C++ calculations, majorly decreasing looping overheads. A large number of R functions are supported by the package (e.g. `colMeans()`) but inevitably a few have to be implemented manually, or worst case imported from R (e.g. `cor.shrink()`).



## Notes about the code
1. In Shrink(P)Cor : currently I am importing the cor.shrink function from R and run in C++. Altought the results are identical, this is not speed optimized. The optimal solution would be using a pure C++ function 
