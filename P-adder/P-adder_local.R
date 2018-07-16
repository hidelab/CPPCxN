p_adder <- function (base,base_gs_n,new_gs_n) {
    # We need a base pcxn result/matrix
    
    options(stringsAsFactors = F)

    # ==== example inputs ====
    # base <- "pcxn_base.RDS"
    #base_gs_n <- "DPD.Hs.gs.mini.PDN.CMAP.RDS"
    #new_gs_n <- "new_gs.RDS"
    
    # Load base matrix and new geneset file
    base_matrix <- readRDS(base)
    base_gs <- readRDS(base_gs_n)
    new_gs <- readRDS(new_gs_n)
    
    # Changing names for testing
    names(new_gs)[1]<-"Random name"
    names(new_gs)[12]<-"Random name 2"
    
    # ==== Find common elements to check whether the genesets have identical members
    common_names <- intersect(names(base_gs),names(new_gs))
    diffs <- c()
    identical_genesets <- TRUE
    for(n in common_names) {
        if(new_gs[[n]] != base_gs[[n]]) {
            identical_genesets <- FALSE
            #print(paste(n , " has different gene members in the two gene set files!",sep = ""))
            diffs <- c(diffs,n)
            print("BOM")
        }
    }
    
    # Checking if identical
    if (!identical_genesets) {
        cl <- paste(c(diffs), collapse=",")
        stop(paste("The following gene sets have different gene members in the two gene set files:", cl, sep = ""))
        #print(diffs)

    }
    
    
    
 
}