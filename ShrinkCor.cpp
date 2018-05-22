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
NumericVector ShrinkCor_cpp(NumericVector x, NumericVector y, int method, Function cshrink){
    // Wrapper to estimate the correlation coefficient between x and y using the
    // shrinkage estimator from Schaffer & Strimmer (2005) [corpcor package]
    // and the corresponding t-statistic and p-value
    // Args
    //  x: a vector with n observations
    //  y: a vector with n observations
    //  method: 1 = Pearson and 2 = Spearman correlation coefficient
    // Returns
    //  A named vector with the correlation estimate, the sample size n, 
    //  the t-statistic and its corresponding p-value
    
// String messages
//    Rcpp::StringVector myvector(3);
     //std::string my_string = "Sample size: ";
//    Rcout << my_string << std::endl;
    
    // Loading R environment to retrieve some variables 
    Environment R_Env = Environment::global_env();
    
    // Printing arguments - CHECK
    //std::string my_string = "Arg 1_C++ : ";
    //Rcout << my_string << std::endl;
    //Rf_PrintValue(x);
    //Rf_PrintValue(y);
    
    // Get sample size and exit if x and y are not of same length
    int n;
    if(x.size() == y.size()){
        n = x.size();
    }else{
        stop("Doh! x and y have different lengths!");
    }
    // Pearson
    if(method == 1){
        
        // Running cor.shrink() using embeddedR
        ///*** R
        //    estimate_1 <- cor.shrink(cbind(x2,y2),verbose=F)
        //*/
        
        NumericMatrix estimate_1 = cshrink(cbind(x,y));
        
        // Fetching the result of cor.shrink
        //NumericMatrix estimate_1 = R_Env["estimate_1"];
        
        //std::string my_string = "Grab_C++ : ";
        //Rcout << my_string << estimate_1(1,0) << std::endl;
        
        double statistic_1 = GetStatistic_cpp(estimate_1(1,0),n);
        double p_value_1 = 2*R::pt(-abs(statistic_1),(n-2),1,0);
        
        //std::string my_1 = "Results 1: ";
        //Rcout << my_1 << std::endl;
        //Rcout << statistic_1 << std::endl;
        //Rcout << p_value_1 << std::endl;
        
        // prepare results
        NumericVector res_1 = NumericVector::create(estimate_1(1,0),n,statistic_1,p_value_1);
        return res_1;
    }
    // Spearman
//    else if (method == 2) {
//        
        // Running cor.shrink() using embeddedR
//        /*** R
//            estimate_2 <- cor.shrink(cbind(rank(x2),rank(y2)),verbose=F)
//        */
//        NumericMatrix estimate_2 = cshrink(cbind(match(x, sort(x)),match(x, sort(y)));
        
        // Fetching the result of cor.shrink
        //NumericMatrix estimate_2 = R_Env["estimate_2"];
//        double statistic_2 = GetStatistic_cpp(estimate_2(1,0),n);
//        double p_value_2 = 2*R::pt(-abs(statistic_2),(n-2),1,0);
        
        // prepare results
//        NumericVector res_2 = NumericVector::create(estimate_2(1,0),n,statistic_2,p_value_2);
        
//        return res_2;
//    }
    // Wrong method input
    else {
        stop("Doh! Invalid method input!");
    }
}