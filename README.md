# (expertly named) CPPCxN
This version of PCxN is optimized for speed by adjusting code structure and using C++. It does not change the idea or results of the original PCxN method found [here](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006042).



## Notes about the code
1. In Shrink(P)Cor : currently I am importing the cor.shrink function from R and run in C++. Altought the results are identical, this is not speed optimized. The optimal solution would be using a pure C++ function 
