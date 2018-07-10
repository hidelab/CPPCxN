# Aggregate results into a single data frame
# We need estimates02 to be completed for both runs


rm(list = ls(all=TRUE)) 
options(stringsAsFactors = F)

# ==== Arguments ====
# command line arguments from submission (sharc) script
# 1 - genesets file
cmd_args <- commandArgs(trailingOnly = T)

# ==== INPUTS ====
m1 <- "pcxn_conc_base"
m2 <- "pcxn_conc_plus_10"
output_folder <- "conc_output_mini_PDN"

# empty variable to store results
pcxn = c()
pb = txtProgressBar(min=0,max=2,initial=0,style=3)
for(k in 1:2){
    # read each part of each
    tmp1 <- readRDS(paste0(output_folder, "/mean_pcor2_barcode/res/",m1 ,"_part", k, ".RDS"))
    tmp2 <- readRDS(paste0(output_folder, "/mean_pcor2_barcode/res/",m2 ,"_part", k, ".RDS"))
    conc_tmp <- rbind(tmp1,tmp2)
    conc_tmp <- unique(conc_tmp)
    pcxn = rbind(pcxn, conc_tmp)
    setTxtProgressBar(pb,k)
}
close(pb)

# adjust p-values for multiple comparison
pcxn$p.Adjust = p.adjust(p = pcxn$p.value, method = "fdr")

# save results
saveRDS(pcxn, paste0(output_folder,"/",m1,"_",m2,".RDS"))