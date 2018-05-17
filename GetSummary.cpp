#include <cstdio>
#include <algorithm>
#include <iostream>
#include <String.h>
#include <cmath>
#include <String.h>
#include <Rcpp.h> 
#include <Rmath.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
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