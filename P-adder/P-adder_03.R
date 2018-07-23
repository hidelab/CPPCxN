# Aggregate results into a single data frame

rm(list = ls(all=TRUE)) 
options(stringsAsFactors = F)

# ==== Arguments ====
# command line arguments from submission (sharc) script
# 1 - genesets file
# 2 - old genesets folder (include "/" in the end)
cmd_args <- commandArgs(trailingOnly = T)

# ==== INPUTS ====
geneset_file <- cmd_args[1]
output_folder <- "output_adder_PCxN"
old_output_folder <- cmd_args[2]

new_files <- list.files(path = paste0("../",output_folder,"/mean_pcor2_barcode/res/") ,pattern="pcxn_mean_pcor2")
old_files <- list.files(path = old_output_folder ,pattern="pcxn_mean_pcor2")

#length(number_of_files)

# empty variable to store results
pcxn = c()

# Going through new files
for(n in new_files){
  pcxn = rbind(pcxn, readRDS(paste0("../",output_folder,"/mean_pcor2_barcode/res/",n)))
}

# Going through new files
for(o in old_files){
    pcxn = rbind(pcxn, readRDS(paste0(old_output_folder,o)))
}

# adjust p-values for multiple comparison
pcxn$p.Adjust = p.adjust(p = pcxn$p.value, method = "fdr")

# save results
saveRDS(pcxn, paste0("../",output_folder,"/adder_PCxN_",geneset_file))

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

saveRDS(pcxn, paste0("../",output_folder,"/square_adder_PCxN_",geneset_file))
