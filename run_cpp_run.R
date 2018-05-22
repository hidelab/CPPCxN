require(Rcpp)
require(RcppGSL)
library(RcppArmadillo)
library('inline')
library("microbenchmark")
library(corpcor)

# Example data for OverlapCoefficient
y1 <- runif(10, min=0, max=200)
x1 <- runif(10, min=0, max=200) 

# Example data for GetSummary
x_matrix <- matrix(1:100, nrow = 10, dimnames = list(c("gene1","gene2","gene3","gene4","gene5","gene6","gene7","gene8","gene9","gene10"), c("s1","s2","s3","s4","s5","s6","s7","s8","s9","s10")))
x_gs <- c("gene1","gene3","gene5","gene9")

# Example data for ShrinkCor !!! these variables are imported to ShrinkCor.cpp
# DO NOT CHANGE NAMES
y2 <- runif(1500, min=0, max=2000)
x2 <- runif(1500, min=0, max=2000) 

# Example data for ShrinkPCor !!! these variables are imported to ShrinkPCor.cpp
# DO NOT CHANGE NAMES
y3 <- runif(1500, min=0, max=2000)
x3 <- runif(1500, min=0, max=2000) 
z3 <- runif(1500, min=0, max=2000) 



# R functions
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

ShrinkCor = function(x,y,method="pearson"){
    # wrapper to estimate the correlation coefficient between x and y using the 
    # shrinkage estimator from Schaffer & Strimmer (2005) [corpcor package]
    # and the corresponding t-statistic and p-value
    #
    # Args
    #   x: a vector with n observations
    #   y: a vector with n observations
    #   method: character to pick either the Pearson or Spearman correlation coefficient
    #
    # Returns
    #   a named vector with the correlation estimate, the sample size n, the t-statistic
    #   and its corresponding p-value
    
    # function to get the t-statistic
    GetStatistic <- function(r,n){
        r*sqrt((n-2)/(1-r^2))
    }
    # get sample size
    if(length(x) == length(y)){
        n <- length(x)
    }else{
        cat("\n x and y have different lengths! >=( \n")
        return(NA)
    }
    # determine method
    selected_method <- match(method,c("pearson","spearman"))
    # Pearson correlation
    if(selected_method == 1){
        estimate <- cor.shrink(cbind(x,y),verbose=F)
        #print(paste0("Statistic_R : ", r*sqrt((n-2)/(1-r^2))))
        
        statistic <- GetStatistic(estimate[2,1],n)
        p.value <- 2*pt(-abs(statistic),n-2)
    }else if(selected_method == 2){
        estimate <- cor.shrink(cbind(rank(x),rank(y)),verbose=F)
        statistic <- GetStatistic(estimate[2,1],n)
        p.value <- 2*pt(-abs(statistic),n-2)
    }else{
        cat("invalid method! >=( \n")
        return(NA)
    }
    # prepare results
    res <- c(estimate[2,1],n,statistic,p.value)
    names(res) <- c("estimate","n","statistic","p.value")
    return(res)
}

ShrinkPCor <- function(x,y,z,method="pearson"){
    # wrapper to estimate the partial correlation coefficient x,y|z using the 
    # shrinkage estimator from Schaffer & Strimmer (2005) [corpcor package]
    # and the corresponding t-statistic and p-value
    #
    # Args
    #   x: a vector with n observations
    #   y: a vector with n observations
    #   z: a vector with n observations
    #   method: character to pick either the Pearson or Spearman partial correlation coefficient
    #
    # Returns
    #   a named vector with the partial correlation estimate, the sample size n, the t-statistic
    #   and its corresponding p-value
    
    # function to get the t-statistic
    GetStatistic <- function(r,n){
        #print("t-statistic(R): ")
        #print(r*sqrt((n-3)/(1-r^2)))
        r*sqrt((n-3)/(1-r^2))
        }
    # get sample size
    if(length(x) == length(y) & length(z) == length(x)){
        n <- length(x)
    }else{
        cat("x,y and z have different lengths! >=( \n")
        return(NA)
    }
    # determine method
    selected_method <- match(method,c("pearson","spearman"))
    # Pearson correlation
    if(selected_method == 1){
        cor.xyz <- cor.shrink(cbind(x,y,z),verbose=F)
        estimate <- cor2pcor(cor.xyz)[1,2] 
        statistic <- GetStatistic(estimate,n)
        p.value <- 2*pt(-abs(statistic),n-3)
    }else if(selected_method == 2){
        cor.xyz <- cor.shrink(cbind(rank(x),rank(y),rank(z)),verbose=F)
        estimate <- cor2pcor(cor.xyz)[1,2] 
        statistic <- GetStatistic(estimate,n)
        p.value <- 2*pt(-abs(statistic),n-3)
    }else{
        cat("invalid method! >=( \n")
        return(NA)
    }
    # prepare results
    res <- c(estimate,n,statistic,p.value)
    names(res) <- c("estimate","n","statistic","p.value")
    return(res)
}

ProcessElement = function(ic){
    # helper function to get the experiment-level estimates for a 
    # gene-set pair
    i = ceiling((sqrt(8*(ic+1)-7)+1)/2)
    j = ic-choose(floor(1/2+sqrt(2*ic)),2)
    
    # pathway gene sets
    gsA=gs_lst[[i]]
    gsB=gs_lst[[j]]
    
    # shared genes
    gsAB <- intersect(gsA,gsB)
    
    
    # get correlation between the summaries for the unique genes
    tmp = data.frame(Pathway.A=names(gs_lst)[i],Pathway.B=names(gs_lst)[j])
    
    # get pathway summaries for disjoint gene sets
    summaryA = GetSummary(dat=exprs_rnk,gs=gsA,sum_fun=colMeans)
    summaryB = GetSummary(dat=exprs_rnk,gs=gsB,sum_fun=colMeans)
    
    if(length(gsAB) > 0){
        # if pathways share genes, estimate conditional correlation (on shared genes)
        summaryAB = GetSummary(dat=exprs_rnk,gs=gsAB,sum_fun=colMeans)
        tmp = c(tmp,ShrinkPCor(
            x=summaryA,
            y=summaryB,
            z=summaryAB,
            method = "spearman"
        ))
    }else{
        # otherwise, estimate correlation between gene sets
        tmp = c(tmp,ShrinkCor(
            x=summaryA,
            y=summaryB,
            method = "spearman"
        ))
    }
    
    # calculate overlap coefficient
    tmp$Overlap.Coeff= OverlapCoefficient(gs_lst[[i]],gs_lst[[j]])
    
    setTxtProgressBar(pb,ic)
    return(tmp)
}


Rcpp::sourceCpp("OverlapCoefficient.cpp")
Rcpp::sourceCpp("GetSummary.cpp")
Rcpp::sourceCpp("ShrinkCor.cpp",embeddedR = TRUE)
Rcpp::sourceCpp("ShrinkPCor.cpp",embeddedR = TRUE)
Rcpp::sourceCpp("ProcessElement.cpp")

# Benchmarking
results_OverlapCoefficient <- microbenchmark(OverlapCoefficient = OverlapCoefficient(x1,y1), OverlapCoefficient_cpp = OverlapCoefficient_cpp(x1,y1))
results_summary_OverlapCoefficient <- summary(results_OverlapCoefficient)[, c(1:7)]
print(summary(results_OverlapCoefficient)[, c(1:7)],digits=1)

results_GetSummary <- microbenchmark(GetSummary = GetSummary(x_matrix,x_gs,sum_fun=colMeans), GetSummary_cpp = GetSummary_cpp(x_matrix,x_gs,1))
results_summary_GetSummary <- summary(results_GetSummary)[, c(1:7)]
print(summary(results_GetSummary)[, c(1:7)],digits=1)

results_ShrinkCor <- microbenchmark(ShrinkCor = ShrinkCor(x2,y2), ShrinkCor_cpp = ShrinkCor_cpp(x2,y2,1,cor.shrink))
results_summary_ShrinkCor <- summary(results_ShrinkCor)[, c(1:7)]
print(summary(results_ShrinkCor)[, c(1:7)],digits=1)

results_ShrinkPCor <- microbenchmark(ShrinkPCor = ShrinkPCor(x3,y3,z3), ShrinkPCor_cpp = ShrinkPCor_cpp(x3,y3,z3,1,cor.shrink,cor2pcor))
results_summary_ShrinkPCor <- summary(results_ShrinkPCor)[, c(1:7)]
print(summary(results_ShrinkPCor)[, c(1:7)],digits=1)

