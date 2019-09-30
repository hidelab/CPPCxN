# Aggregate results into a single data frame

rm(list = ls(all=TRUE)) 
options(stringsAsFactors = F)

# ==== Arguments ====
# command line arguments from submission (sharc) script
# 1 - genesets file
# 2 - output folder

setwd("/Users/kkoler/Downloads/")

# ==== INPUTS ====
geneset_file <- "MSigDB.L1000CDS2.staticmod.6.subset.RDS"

# empty variable to store results
pcxn = c()
pb = txtProgressBar(min=0,max=2,initial=0,style=3)

parts <- readRDS(paste0("res/parts.RDS"))

for(k in 1:parts[[1]]){
    pcxn = rbind(pcxn, readRDS(paste0("res/pcxn_mean_pcor2_barcode_part",k,".RDS")))
    setTxtProgressBar(pb,k)
}
close(pb)
# adjust p-values for multiple comparison
pcxn$p.Adjust = p.adjust(p = pcxn$p.value, method = "fdr")

# save results
saveRDS(pcxn, paste0("improved_PCxN_",geneset_file))
