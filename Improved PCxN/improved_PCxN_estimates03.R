# Aggregate results into a single data frame

rm(list = ls(all=TRUE)) 
options(stringsAsFactors = F)

# ==== Arguments ====
# command line arguments from submission (sharc) script
# 1 - genesets file
cmd_args <- commandArgs(trailingOnly = T)

# ==== INPUTS ====
geneset_file <- cmd_args[1]
output_folder <- "output_improved_PCxN"

# empty variable to store results
pcxn = c()
pb = txtProgressBar(min=0,max=2,initial=0,style=3)
for(k in 1:2){
  pcxn = rbind(pcxn, readRDS(paste0("../",output_folder,"/mean_pcor2_barcode/res/pcxn_mean_pcor2_barcode_part",k,".RDS")))
  setTxtProgressBar(pb,k)
}
close(pb)

# adjust p-values for multiple comparison
pcxn$p.Adjust = p.adjust(p = pcxn$p.value, method = "fdr")

# save results
saveRDS(pcxn, paste0("../",output_folder,"/improved_PCxN_",geneset_file))

original_matrix <- pcxn

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

saveRDS(pcxn, paste0("../",output_folder,"/square_improved_PCxN_",geneset_file))
