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
    double result = r*(sqrt((n-3)/(1-pow(r,2))));
    return result;
}

// [[Rcpp::export]]
NumericVector ShrinkPCor_cpp(NumericVector x, NumericVector y, 
                            NumericVector z, int method, Function cshrink, Function cor2pcor){
    // Wrapper to estimate the partial correlation coefficient x,y|z using the 
    // shrinkage estimator from Schaffer & Strimmer (2005) [corpcor package]
    // and the corresponding t-statistic and p-value
    // Args
    //  x: a vector with n observations
    //  y: a vector with n observations
    //  z: a vector with n observations
    //  method: 1 = Pearson and 2 = Spearman correlation coefficient
    // Returns
    //  A named vector with the correlation estimate, the sample size n, 
    //  the t-statistic and its corresponding p-value
    
    // Loading R environment to retrieve some variables 
    Environment R_Env = Environment::global_env();
    
    // Get sample size and exit if x and y are not of same length
    int n;
    if((x.size() == y.size()) && (x.size() == z.size())){
        n = x.size();
        //Rcout << my_string << n << std::endl;
    }else{
        stop("Doh! x,y and z have different lengths!");
    }
    // Pearson
    if(method == 1){
        
        // Running cor.shrink() using embeddedR
        ///*** R
        //    estimate_P_1 <- 0
        //    cor.xyz <- cor.shrink(cbind(x3,y3,z3),verbose=F)
        //    estimate_P_1 <- cor2pcor(cor.xyz)[1,2] 
        //   statistic <- GetStatistic(estimate_P_1,n)
        //    p.value <- 2*pt(-abs(statistic),n-3)
        //*/
        
        NumericMatrix cor_xyz = cshrink(cbind(x,y,z));
        NumericMatrix estimate_P_1 = cor2pcor(cor_xyz) ;
        double statistic_P_1 = GetStatistic_cpp(estimate_P_1(0,1),n);
        double p_value_P_1 = 2*R::pt(-abs(statistic_P_1),(n-3),1,0);
        
        // prepare results
        NumericVector res_P_1 = NumericVector::create(estimate_P_1(0,1),n,statistic_P_1,p_value_P_1);
        return res_P_1;
    }
    // Spearman
    else if (method == 2) {
        
        // Obtain environment containing function
        Rcpp::Environment base("package:base");
        Function rank = base["rank"];
        
        NumericMatrix cor_xyz = cshrink(cbind(rank(x),rank(y),rank(z)));
        NumericMatrix estimate_P_2 = cor2pcor(cor_xyz) ;
        double statistic_P_2 = GetStatistic_cpp(estimate_P_2(0,1),n);
        double p_value_P_2 = 2*R::pt(-abs(statistic_P_2),(n-3),1,0);
        
        // prepare results
        NumericVector res_P_2 = NumericVector::create(estimate_P_2(0,1),n,statistic_P_2,p_value_P_2);
        return res_P_2;
    }
    // Wrong method input
    else {
        stop("Doh! Invalid method input!");
    }
}