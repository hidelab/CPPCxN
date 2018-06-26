// [[Rcpp::depends(RcppParallel)]]

#include <cstdio>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <String.h>
#include <Rcpp.h> 
#include <Rmath.h>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace RcppParallel;
using namespace std;

NumericMatrix GetSummary_cpp(NumericMatrix dat, CharacterVector gs){

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
            mean_mat(0,_) = colMeans(mat);
            return mean_mat;
    }
    else {
        return mat;
    }
}
/*
struct Summarizer : public Worker {
    // source matrix (dat) and the gene-sets set
    const RMatrix<double> input;
    const RVector<char> jl;
    
    // destination matrix (=summary_joint_matrix_cpp)
    RMatrix<double> summary_joint_matrix_cpp;
    
    
    // initialize with source and destination
    Summarizer(const NumericMatrix input, 
               const List jl, 
               NumericMatrix summary_joint_matrix_cpp) 
        : input(input), jl(jl) , summary_joint_matrix_cpp(summary_joint_matrix_cpp) {}
    
    
    
    // take the summary for the specific range of elements passed down from the call
    void operator()(std::size_t begin, std::size_t end) {
        
        for (std::size_t i = begin; i < end; i++) {
            CharacterVector gs = jl[i];
            Rcout << gs << std::endl;
            
        }
            //summary_joint_matrix_cpp(it,_) = GetSummary_cpp(input, gs);
        
        
        // std::transform(jl.begin() + begin, 
        //                jl.begin() + end, 
        //                summary_joint_matrix_cpp.begin() + begin, 
        //                ::GetSummary_cpp);
    }
};

*/


// [[Rcpp::export]]
NumericMatrix create_disjoint_matrix(List jl, NumericMatrix dat) {
    
    //dat = dat;
    
    // colnames(samples) of dat
    CharacterVector dat_colnames = colnames(dat);
    int cols = dat_colnames.size();
    
    List names = jl.names();
    const int n = names.size();
    
    //CharacterVector cv = Rcpp::as< CharacterVector >(jl);
    
    //Rcout << cv[1] << std::endl;
    
    // Instead of list
    // NumericMatrix summary_joint_matrix_cpp(n, cols);
    // rownames(summary_joint_matrix_cpp) = names;
    // colnames(summary_joint_matrix_cpp) = dat_colnames;
    
    // allocate the output matrix
    NumericMatrix summary_joint_matrix_cpp(n, cols);
    rownames(summary_joint_matrix_cpp) = names;
    colnames(summary_joint_matrix_cpp) = dat_colnames;
    
    // Summarizer functor (pass input and output matrixes)
    //Summarizer summarizer(dat, jl, summary_joint_matrix_cpp);
    
    // call parallelFor to do the work
    //parallelFor(0, n, summarizer);
    
    for(int i=0; i < n; i++) {
        CharacterVector gs = jl[i];
        summary_joint_matrix_cpp(i,_) = GetSummary_cpp(dat=dat,gs=gs);
        double perc = i/n;
        Rcout << "Done: " << perc << "\r";
    }
    
    return summary_joint_matrix_cpp;  
}
