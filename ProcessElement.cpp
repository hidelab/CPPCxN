// [[Rcpp::depends(corpcor)]]

#include <cstdio>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <String.h>
#include <Rcpp.h> 
#include <Rmath.h>

using namespace Rcpp;
using namespace std;

double GetStatistic_cpp(double r,double n) {
    double result = r*(sqrt((n-2)/(1-pow(r,2))));
    return result;
}

// [[Rcpp::export]]
NumericVector ShrinkCor_cpp(NumericVector ic){
    // Helper function to get the experiment-level estimates for a 
    // gene-set pair
    
    
    
    
    
 return 0;
}