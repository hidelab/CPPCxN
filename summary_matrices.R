# Run this every time you change genesets (e.g. the file DPD.Hs.gs.xxx.PDN.RDS)

# ================== libs =================
library(svd)
library(corpcor)
library(parallel)
require(Rcpp)
require(RcppGSL)
library(RcppArmadillo)
library('inline')
library("microbenchmark")

# ================ INPUTS ================
geneset_file <- "DPD.Hs.gs.mini.PDN.RDS"
output_folder <- "CPP_output_mini_PDN"

# directory with gene expression background
barcode_dir <- "data/HGU133plus2/"

# ==== PCxN C++ Functions ====
Rcpp::sourceCpp("create_disjoint_matrix.cpp")
# cppFunction("NumericMatrix create_disjoint_matrix(List jl, NumericMatrix dat) {
#             
#             // colnames(samples) of dat
#             CharacterVector dat_colnames = colnames(dat);
#             int cols = dat_colnames.size();
#             
#             List names = jl.names();
#             const int n = names.size();
#             
#             // Instead of list
#             NumericMatrix summary_joint_matrix_cpp(n, cols);
#             rownames(summary_joint_matrix_cpp) = names;
#             colnames(summary_joint_matrix_cpp) = dat_colnames;
#             
#             for(int i=0; i < n; i++) { 
#             CharacterVector gs = jl[i];
#             summary_joint_matrix_cpp(i,_) = GetSummary_cpp(dat=dat,gs=gs,1);
#             }
#             
#             return summary_joint_matrix_cpp;  
#             }"
#         ,includes='NumericMatrix GetSummary_cpp(NumericMatrix dat, CharacterVector gs, int sum_fun){
#         // Args.
#         //  dat: genes by samples matrix
#         //  gs: vector with the names of the genes in the gene set
#         //  sum_fun: function to calculate the summary
#         // Returns
#         //  A 1 by samples vector with the summary statistic for the pathway
#         
#         // rownames(genes) of dat
#         CharacterVector dat_rownames = rownames(dat);
#         int rows = dat_rownames.size();
#         
#         // colnames(samples) of dat
#         CharacterVector dat_colnames = colnames(dat);
#         int cols = dat_colnames.size();
#         
#         // dat genes that are also in gs => dat_gs[]
#         CharacterVector dat_gs(gs.size());
#         NumericVector dat_gs_pos(gs.size());
#         int in_gs_counter = 0;
#         for (int j=0; j < rows; j++) {
#         for (int i=0; i < gs.size(); i++) {
#         // checking if matrix gene is in gs
#         if ((dat_rownames[j] == gs[i]) == 1) {
#         dat_gs[in_gs_counter] = dat_rownames[j];
#         dat_gs_pos[in_gs_counter] = j;
#         in_gs_counter++;
#         break;
#         }
#         }
#         }
#         
#         // pre-allocate result matrix
#         NumericMatrix mat(in_gs_counter,cols);
#         rownames(mat) = dat_gs;
#         colnames(mat) = dat_colnames;
#         
#         CharacterVector result_colnames;
#         
#         // construct result matrix
#         for (int m=0; m < in_gs_counter; m++) {
#         mat(m,_) = dat( dat_gs_pos[m],_);
#         }
#         
#         //Rf_PrintValue(dat_gs);
#         //Rf_PrintValue(dat_gs_pos);
#         
#         // pre-allocate result 1-dim means matrix
#         NumericMatrix mean_mat(1,cols);
#         colnames(mean_mat) = dat_colnames;
#         //NumericMatrix median_mat(1,cols);
#         //colnames(median_mat) = dat_colnames;
#         
#         // Non-empty input matrix
#         if (gs.size() > 1) {
#         // Case: colMeans
#         if(sum_fun == 1) {
#         
#         mean_mat(0,_) = colMeans(mat);
#         return mean_mat;
#         } 
#         // Case: any other 
#         else {
#         return mat;
#         }
#         }
#         else {
#         return mat;
#         }
#         }')
    
# ==== PCxN Functions ====
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

for (current_tissue in 1:length(res)){
    # ---- tissue specifics ----
    # select tissue type
    tissue_select <- names(res)[ current_tissue ]
    print(paste("Starting ", tissue_select, " tissue...", sep = ""))
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
    
    # ---- custom ----
    #gs_lst <- gs_lst[1:100]
    
    # indices for pathway pairs 
    number_of_pathways = choose(length(gs_lst),2)
    input = 1:(number_of_pathways)
    
    # loop thru each experiment (GSE series) for a given tissue type
    for(j in seq_along(tissue_series)){
        
        # subset expression ranks
        seriesn <- tissue_seriesn[j]
        series <- tissue_series[j]
        
        print(paste("> Calculating ", series, " tissue...", sep = ""))
        
        ind_series <- colnames(tissue_rnk) %in% tissue_meta$sample[tissue_meta$series == seriesn]
        exprs_rnk <- tissue_rnk[ ,ind_series ]
        
        # ---- Get pathway summaries for joint gene sets ----
        nms <- combn( names(gs_lst) , 2 , FUN = paste0 , collapse = " - " , simplify = FALSE )
        ll <- combn( gs_lst , 2 , simplify = FALSE )
        out <- lapply( ll , function(x) intersect( x[[1]] , x[[2]] ) )
        joint_gs_lst <- setNames( out , nms )
        joint_gs_lst <- joint_gs_lst[lapply(joint_gs_lst,length)>0]


        summary_joint_list <- create_disjoint_matrix(joint_gs_lst,exprs_rnk)
        summary_joint_list<- apply(summary_joint_list, 1, as.list)
        out.name <- paste(output_folder,"/Summaries/summary_joint_list_",tissue_select, "_", series, ".RData",sep = "")
        save(summary_joint_list, file=out.name)
        
        # ---- Get pathway summaries for disjoint gene sets ----
        # summary_dis_list <-list()
        # 
        # for (gs in 1:length(gs_lst)) {
        #     nm <- names(gs_lst[gs])
        #     tmp <- list(GetSummary(dat=exprs_rnk,gs=gs_lst[[gs]],colMeans))
        #     summary_dis_list[[nm]] <- tmp
        # }
        # 
        # out.name.dis <- paste(output_folder,"/Summaries/summary_dis_list_",tissue_select, "_", series, ".RData",sep = "")
        # save(summary_dis_list, file=out.name.dis)
    }
    
}

# ==== Get pathway summaries for joint gene sets ====
# nms <- combn( names(gs_lst) , 2 , FUN = paste0 , collapse = " - " , simplify = FALSE )
# ll <- combn( gs_lst , 2 , simplify = FALSE )
# out <- lapply( ll , function(x) intersect( x[[1]] , x[[2]] ) )
# joint_gs_lst <- setNames( out , nms )
# joint_gs_lst <- joint_gs_lst[lapply(joint_gs_lst,length)>0]
# 
# summary_joint_list <- list()
# 
# for (gs in 1:length(joint_gs_lst)) {
#     nm <- names(joint_gs_lst[gs])
#     tmp <- list(GetSummary(dat=exprs_rnk,gs=joint_gs_lst[[gs]],colMeans))
#     summary_joint_list[[nm]] <- tmp
# }


# ---- Get pathway summaries for disjoint gene sets ----
# summary_dis_list <-list()
# 
# for (gs in 1:length(gs_lst)) {
#     nm <- names(gs_lst[gs])
#     tmp <- list(GetSummary(dat=exprs_rnk,gs=gs_lst[[gs]],colMeans))
#     summary_dis_list[[nm]] <- tmp
# }