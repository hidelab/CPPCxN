#Rcpp::sourceCpp("step1.cpp")

require(Rcpp)
require(RcppGSL)
library(RcppArmadillo)

library('inline')

rcpp_inc <- '
using namespace Rcpp;
using namespace arma;
'

# Test data (vectors etc)
y <- runif(10, min=0, max=200)
x <- runif(10, min=0, max=200) 

x_matrix <- matrix(1:100, nrow = 10, dimnames = list(c("gene1","gene2","gene3","gene4","gene5","gene6","gene7","gene8","gene9","gene10"), c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10")))
x_gs <- c("gene1","gene3","gene5","gene9")
    
# Benchmarking
library("microbenchmark")

results <- microbenchmark(OverlapCoefficient = OverlapCoefficient(x,y), OverlapCoefficient_cpp = OverlapCoefficient_cpp(x,y))
print(summary(results)[, c(1:7)],digits=1)

results <- microbenchmark(GetSummary = GetSummary(x_matrix,x_gs,sum_fun=colMeans), GetSummary_cpp = GetSummary_cpp(x_matrix,x_gs,1))
print(summary(results)[, c(1:7)],digits=1)


# --- OverlapCoefficient  START ---
cppFunction("
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
")

OverlapCoefficient <- function(x,y){
    # function to calculate the overlap coefficient between x and y
    # which is defined as the size of the intersection divided by the
    # size of the smaller set
    #
    # Args
    #   x: a vector
    #   y: a vector
    #
    # Returns
    #   the overlap coefficient, a number between 0 and 1
    
    length(intersect(x,y))/min(length(unique(x)),length(unique(y)))
}
# --- OverlapCoefficient  END ---

# --- GetSummary START ---
cppFunction("
    #include <cstdio>
    #include <algorithm>
    #include <iostream>
    #include <String.h>

    NumericMatrix GetSummary_cpp(NumericMatrix dat, CharacterVector gs, int sum_fun){
        // Function to calculate the summary statistic for the pathway
        // Args.
        //  dat: genes by samples matrix
        //  gs: vector with the names of the genes in the gene set
        //  sum_fun: function to calculate the summary
        // Returns
        //  A 1 by samples vector with the summary statistic for the pathway
        
        // rownames(genes) of dat
        CharacterVector dat_rownames = rownames(dat);
        int rows = dat_rownames.size();

        // colnames(samples) of dat
        CharacterVector dat_colnames = colnames(dat);
        int cols = dat_colnames.size();

        // dat genes that are also in gs => dat_gs[]
        CharacterVector dat_gs(gs.size());
        NumericVector dat_gs_pos(gs.size());
        int in_gs_counter = 0;
        for (int j=0; j < rows; j++) {
            for (int i=0; i < gs.size(); i++) {
                // checking if matrix gene is in gs
                if ((dat_rownames[j] == gs[i]) == 1) {
                    dat_gs[in_gs_counter] = dat_rownames[j];
                    dat_gs_pos[in_gs_counter] = j;
                    in_gs_counter++;
                    break;
                }
            }
        }

        // pre-allocate result matrix
        NumericMatrix mat(in_gs_counter,cols);
        rownames(mat) = dat_gs;
        colnames(mat) = dat_colnames;

        CharacterVector result_colnames;

        // construct result matrix
        for (int m=0; m < in_gs_counter; m++) {
            mat(m,_) = dat( dat_gs_pos[m],_);
        }

        //Rf_PrintValue(dat_gs);
        //Rf_PrintValue(dat_gs_pos);

        // pre-allocate result 1-dim means matrix
        NumericMatrix mean_mat(1,cols);
        colnames(mean_mat) = dat_colnames;
        //NumericMatrix median_mat(1,cols);
        //colnames(median_mat) = dat_colnames;

        // Non-empty input matrix
        if (gs.size() > 1) {

            // Case: colMeans
            if(sum_fun == 1) {
                
                mean_mat(0,_) = colMeans(mat);
                return mean_mat;
            } 
            // Case: any other 
            else {
                return mat;
            }
        }

        else {
            return mat;
        }
    }

")

GetSummary = function(dat,gs,sum_fun){
    # function to calculate the summary statistic for the pathway
    #
    # Args.
    #   dat: genes by samples matrix
    #   gs: vector with the names of the genes in the gene set
    #   sum_fun: function to calculate the summary
    #
    # Returns
    #   a 1 by samples vector with the summary statistic for the pathway
    
    if(length(gs) > 1){
        # calculate summary for pathways with more than 1 element
        return(sum_fun(dat[rownames(dat) %in% gs,]))
    }else{
        # return actual value for pathways with a single element
        return(dat[rownames(dat) %in% gs,])
    }
}
# --- GetSummary END ---


print("C++ runtime:")
system.time(OverlapCoefficient_cpp(y,x))
print("--------------------------------")
print("R runtime:")
system.time(OverlapCoefficient(y,x))
print("--------------------------------")

#OverlapCoefficient_cpp(length(x),length(y))