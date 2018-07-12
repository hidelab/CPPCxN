# Aggregate results into a single data frame
# We need estimates02 to be completed for both runs


rm(list = ls(all=TRUE)) 
options(stringsAsFactors = F)

# ==== Arguments ====
# command line arguments from submission (sharc) script
# 1 - genesets file
cmd_args <- commandArgs(trailingOnly = T)

# ==== INPUTS ====
inputs <- c("pcxn_conc_base","pcxn_conc_plus_10","pcxn_conc_plus_20")
output_name <- "pcxn_conc_base_plus_10_plus_20"

# Check if all necessary files exists in this directory
for(i in inputs) {
    tname1 <- paste(i,"_part1.RDS",sep = "")
    tname2 <- paste(i,"_part2.RDS",sep = "")
    
    if(!file_test("-f", tname1)) stop(paste(tname1, ": No such file found in this directory!",sep = ""))
    if(!file_test("-f", tname2)) stop(paste(tname2, ": No such file found in this directory!",sep = ""))
}

# empty variable to store results
pcxn = c()
pb = txtProgressBar(min=0,max=2,initial=0,style=3)

for(k in 1:2){
    # read each part of each
    conc_tmp <- c()
    for (m in inputs) {
        tmp <- readRDS(paste0(m ,"_part", k, ".RDS"))
        conc_tmp <- rbind(conc_tmp,tmp)
    }
    
    conc_tmp <- unique(conc_tmp)
    pcxn = rbind(pcxn, conc_tmp)
    setTxtProgressBar(pb,k)
}
close(pb)

# adjust p-values for multiple comparison
pcxn$p.Adjust = p.adjust(p = pcxn$p.value, method = "fdr")

# save results
saveRDS(pcxn, paste0(output_name,".RDS"))