convert_pcxn_result <- function(pcxn_res) {
    # Original PCxN's result
    original_matrix <- readRDS(pcxn_res)
    
    unique_pathways <- unique(c(original_matrix$Pathway.A,original_matrix$Pathway.B))
    pathway_number <- length(unique_pathways)
    rns <- unique_pathways
    cns <- unique_pathways
    
    square_matrix <- matrix(, pathway_number, pathway_number, dimnames = list(rns, cns))
    
    for(row in 1:nrow(original_matrix)) {
        tp1 <- original_matrix[row, 1]
        tp2 <- original_matrix[row, 2]
        square_matrix[tp1,tp2] <- original_matrix[row, 4]
        square_matrix[tp2,tp1] <- original_matrix[row, 4]
    }
    
    diag(square_matrix) <- "1"
    square_matrix <- `class<-`(square_matrix, 'numeric')
    
    saveRDS(square_matrix, paste("square_", pcxn_res, sep = "") )
}
