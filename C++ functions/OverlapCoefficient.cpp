#include <cstdio>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <String.h>
#include <Rcpp.h> 
#include <Rmath.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
float OverlapCoefficient_cpp(NumericVector ar1, NumericVector ar2){
    // Function to calculate the overlap coefficient between x and y
    // which is defined as the size of the intersection divided by the
    // size of the smaller set
    
    // Args
    // x: a vector
    // y: a vector
    
    // Returns
    // The overlap coefficient, a number between 0 and 1
    
    // Use the length of any of the inputs arrays since the intersection cannot be larger than any of the sets
    int ar1_length = ar1.size();
    int ar2_length = ar2.size();
    
    // Initialize vector that holds the intersection elements
    NumericVector intersection_vector(ar1_length);
    
    // Initialize loop variables
    int i, j, k = 0;
    
    // Using C++ intersect() function
    int intersection_size = (intersect(ar1, ar2)).size();
    
    // Size of smaller set
    int min_length = 0;
    if (ar1_length <= ar2_length) {
        min_length = ar1_length;
    } else {
        min_length = ar2_length;	
    }
    
    float intersection_result = float(intersection_size)/float(min_length);
    
    return intersection_result;
}