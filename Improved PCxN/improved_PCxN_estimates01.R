    # rels: array of desired pairs based on the list below
    # 1. pathway-CMAP
    # 2. pathway-CTD
    # 3. pathway-PharmGKB
    # 4. CMAP-CTD
    # 5. CMAP-PharmGKB
    # 6. CMAP.up-CMAP.down
    
    rm(list = ls(all=TRUE)) 
    gc()
    options(stringsAsFactors = F)
    
    # ==== Arguments ====
    # command line arguments from submission (sharc) script
    # 1 - job index, pick tissue type
    # 2 - number of cores
    # 3 - desired relationships/pairs
    # 4 - genesets file
    cmd_args <- commandArgs(trailingOnly = T)
    
    # ================== libs =================
    library(svd)
    library(corpcor)
    library(parallel)
    #require(Rcpp)
    #require(RcppGSL)
    #library(RcppArmadillo)
    #library('inline')
    #library("microbenchmark")
    
    # ============ C++ functions =============
    #Rcpp::sourceCpp("OverlapCoefficient.cpp")
    #Rcpp::sourceCpp("GetSummary.cpp")
    #Rcpp::sourceCpp("ShrinkCor.cpp",embeddedR = TRUE)
    #Rcpp::sourceCpp("ShrinkPCor.cpp",embeddedR = TRUE)
    # ================ INPUTS ================
    rels_char <- cmd_args[3]
    geneset_file <- cmd_args[4]
    output_folder <- "output_improved_PCxN"

    # directory with gene expression background
    barcode_dir <- "data/HGU133plus2/"
    
    # ==== PCxN Functions ====
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
        GetStatistic <- function(r,n){r*sqrt((n-2)/(1-r^2))}
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
        GetStatistic <- function(r,n){r*sqrt((n-3)/(1-r^2))}
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
    
    # ==== Pathway Annotation ====
    # Filtered gene set annotation
    gs_lst = readRDS(paste("data/",geneset_file, sep="")) 
    
    
    # ==== Barcode Annotation ====
    # Sample annotation for the gene expression background
    tissue_annot <- readRDS( "data/Barcode3.tissue.RDS" )
    
    # ==== GSE Series ====
    # GSE series per tissue
    tmp <- subset(tissue_annot,select=c(tissue,series))
    gse_lst <- lapply(split(tmp, tmp$tissue),function(x){table(x$series)})
    # order by number of samples
    gse_lst <- gse_lst[order(sapply(gse_lst,sum),decreasing=T)]
    rm(tmp)
    # filter series with at least 5 samples
    res <- lapply(gse_lst,function(x){x[which(x >= 10)]})
    res <- res[lapply(res,length)>0]
    # order by number of samples
    res <- res[order(sapply(res,sum),decreasing=T)]
    
    # ==== Gene Expression Data ====
    # get fRMA normalized values
    getExprs <- function(x){
        tissue_fn <- gsub("[#,%:]",".",x)
        # get path to sample of a given tissue
        tissue_rds <- paste0(barcode_dir,tissue_fn,"/",tissue_fn,".collapse.RDS")
        # load normalized expression values
        tissue_exprs <- readRDS(tissue_rds)$datETcollapsed
        return(tissue_exprs)
    }
    
    # select tissue type
    tissue_select <- names(res)[ as.numeric(cmd_args[1]) ]
    # load normalized expression values
    tissue_exprs <- getExprs( tissue_select )
    # tissue rank expression values
    tissue_rnk <- apply(tissue_exprs,2,rank)
    # tissue meta-data
    tissue_meta <- subset(tissue_annot[tissue_annot$tissue == tissue_select,],select=c(sample,series))
    # get tissue GSE series (experiment IDs)
    tissue_series <- names(res[[ tissue_select ]])
    tissue_seriesn <- tissue_series
    tissue_series <- gsub(";","_",tissue_series)
    
    # keep only genes present in given tissue type
    gs_lst = lapply(gs_lst,function(x){x[x %in% rownames(tissue_exprs)]})
    gs_lst = gs_lst[ sapply(gs_lst,length) > 0 ]
    
    # ==== Experiment-level Estimates ====
    
    # Checks if we have a specific interaction: TRUE if relationships is allowed, otherwise FALSE
    check_rel <- function(n1,n2,rel){
        
        switch(rel,
               "1" = if((startsWith(n1, "Pathway.") & startsWith(n2, "CMAP.")) | (startsWith(n2, "Pathway.") & startsWith(n1, "CMAP."))){return(TRUE)}else {return(FALSE)},
               "2" = if((startsWith(n1, "Pathway.") & startsWith(n2, "CTD.")) | (startsWith(n2, "Pathway.") & startsWith(n1, "CTD."))){return(TRUE)}else {return(FALSE)},
               "3" = if((startsWith(n1, "Pathway.") & startsWith(n2, "PharmGKB.")) | (startsWith(n2, "Pathway.") & startsWith(n1, "PharmGKB."))){return(TRUE)}else {return(FALSE)},
               "4" = if((startsWith(n1, "CMAP.") & startsWith(n2, "CTD.")) | (startsWith(n2, "CMAP.") & startsWith(n1, "CTD."))){return(TRUE)}else {return(FALSE)},
               "5" = if((startsWith(n1, "CMAP.") & startsWith(n2, "PharmGKB.")) | (startsWith(n2, "CMAP.") & startsWith(n1, "PharmGKB."))){return(TRUE)}else {return(FALSE)},
               "6" = if((startsWith(n1, "CMAP.up.") & startsWith(n2, "CMAP.down")) | (startsWith(n2, "CMAP.up") & startsWith(n1, "CMAP.down"))){return(TRUE)}else {return(FALSE)}
        )
    }
    
    # helper function to get the experiment-level estimates for a gene-set pair
    ProcessElement = function(ic){
        i = ceiling((sqrt(8*(ic+1)-7)+1)/2)
        j = ic-choose(floor(1/2+sqrt(2*ic)),2)
        
        # pathway gene sets
        gsA=gs_lst[[i]]
        gsB=gs_lst[[j]]

        n1 <- names(gs_lst[i])
        n2 <- names(gs_lst[j])
        
        pass <- TRUE
        
        for (r in 1:length(rels)) {
            if(!check_rel(n1,n2,rels[r])) pass <- FALSE
        }
        
        # Check if this pairs passes all relationship checks
        if(!pass) return(NULL)
        
        # shared genes
        gsAB <- intersect(gsA,gsB)
        
        # get correlation between the summaries for the unique genes
        tmp = data.frame(Pathway.A=names(gs_lst)[i],Pathway.B=names(gs_lst)[j])
        
        if(length(gsAB) > 0){
            # if pathways share genes, estimate conditional correlation (on shared genes)
            summaryAB = GetSummary(dat=exprs_rnk,gs=gsAB,colMeans)
            
            tmp = c(tmp,ShrinkPCor(
                x=unlist(summary_dis_list[[i]]),
                y=unlist(summary_dis_list[[j]]),
                z=summaryAB,
                method = "pearson"
            ))
        }else{
            # otherwise, estimate correlation between gene sets
            tmp = c(tmp,ShrinkCor(
                x=unlist(summary_dis_list[[i]]),
                y=unlist(summary_dis_list[[j]]),
                method = "pearson"
            ))
        }
        
        # calculate overlap coefficient
        tmp$Overlap.Coeff= OverlapCoefficient(as.numeric(gs_lst[[i]]),as.numeric(gs_lst[[j]]))
        
        setTxtProgressBar(pb,ic)
        return(tmp)
        
    }
    
    rels <- as.numeric(unlist(strsplit(rels_char, ",")))
    
    # indices for pathway pairs 
    number_of_pathways = choose(length(gs_lst),2)
    input = 1:(number_of_pathways)
    
    #summary_joint_list <- readRDS("~/TheSandbox/C++ examples/summary_joint_list.RDS")
    
    # Checks if we have a pathway - pathway interaction: TRUE or FALSE
    check_path_path <- function(x){
        str1 <- unlist(strsplit(unlist(x), split = " - "))[1]
        str2 <- unlist(strsplit(unlist(x), split = " - "))[2]
        
        if(startsWith(str1, "Pathway.") & startsWith(str2, "Pathway.")){
            return(TRUE)
        }else {
            return(FALSE)
        }
    }
    
    # loop thru each experiment (GSE series) for a given tissue type
    for(j in seq_along(tissue_series)){
        # subset expression ranks
        seriesn <- tissue_seriesn[j]
        series <- tissue_series[j]
        
        ind_series <- colnames(tissue_rnk) %in% tissue_meta$sample[tissue_meta$series == seriesn]
        exprs_rnk <- tissue_rnk[ ,ind_series ]
        cat(
            "+==========================+\n",
            "Tissue:",tissue_select,"\n",
            "Series:",series,"\n",
            "+==========================+\n\n\n"
        )
        
        # ---- Get pathway summaries for disjoint gene sets ----
        summary_dis_list <-list()
        
        for (gs in 1:length(gs_lst)) {
            nm <- names(gs_lst[gs])
            tmp <- list(GetSummary(dat=exprs_rnk,gs=gs_lst[[gs]],colMeans))
            summary_dis_list[[nm]] <- tmp
        }
        
        # get experiment-level estimates (parallel loop for pathway pairs)
        pb = txtProgressBar(min=0,max=number_of_pathways,style=3,initial=0)
        cat("\n")
        res = mclapply(input,ProcessElement,mc.cores=as.numeric(cmd_args[2]))
        # Remove NULLs
        res[sapply(res, is.null)] <- NULL
        close(pb)
        
        # save experiment-level estimates
        fname = paste0(make.names(tissue_select),"_",series)
        saveRDS(res, paste0(output_folder,"/mean_pcor2_barcode/",fname,"_cpad_pathcor.RDS"))
    }