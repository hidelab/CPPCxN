
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
NumericVector ProcessElement_cpp(int ic){
    // Helper function to get the experiment-level estimates for a gene-set pair
    
    // Obtain environment containing function
    //Rcpp::Environment base("package:base");
    //Function rank = base["choose"];
    
    int i = ceil((sqrt(8*(ic+1)-7)+1)/2);
    
    NumericVector j_inner = NumericVector::create(floor(1/2+sqrt(2*ic)));
    NumericVector j_outer = choose(j_inner,2);
    
    int j = ic - j_outer[0];

    ///*** R
    //    # Pathway gene sets
    //    gsA=gs_lst[[i]]
    //    gsB=gs_lst[[j]]
    //*/
    
    Environment R_Env = Environment::global_env();
    List gs_lst = R_Env["gs_lst"];
    List gsA = gs_lst[i];
    List gsB = gs_lst[j];
    NumericVector gsA_C = wrap(gsA);
    NumericVector gsB_C = wrap(gsB);
    
    Rf_PrintValue(gsA_C);
    
    //List gsB = gs_lst[j];
    
    // Shared genes
    //List gsAB = intersect(gsA,gsB);
    
    //std::string m1 = "i_C++ : ";
    //Rcout << m1 << i << std::endl;
    //std::string m2 = "j_C++ : ";
    //Rcout << m2 << j << std::endl;
    
 return 0;
}