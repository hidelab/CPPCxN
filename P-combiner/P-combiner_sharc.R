# Aggregate results into a single data frame
# We need estimates02 to be completed for both matrices

# How to use
# Place parts1 and 2 (output of estimates 2) of both matrices in the same folder as this script.
# As the first argument use the matrix names (the part before "_part1.RDS") seperated by commas. 
# As the second argument use the entire name of the combined matrix (e.g. matri1_matrix2.RDS)
# Example (files already in folder):
# Rscript P-combiner_sharc.R pcxn_conc_base,pcxn_conc_plus_10,pcxn_conc_plus_20 1,10,20 all_concs

rm(list = ls(all=TRUE)) 
options(stringsAsFactors = F)

# ==== Arguments ====
# command line arguments from submission (sharc) script
# 1 - genesets file
cmd_args <- commandArgs(trailingOnly = T)

# ==== INPUTS ====
inputs <- unlist(strsplit(cmd_args[1], ","))
#inputs <- c("pcxn_conc_base","pcxn_conc_plus_10","pcxn_conc_plus_20")
nparts <- as.numeric(unlist(strsplit(cmd_args[2], ",")))

output_name <- cmd_args[3]


# Check if all necessary files exists in this directory
for(i in 1:length(inputs)) {
	for (k in 1:nparts[i]){
		tname1 <- paste0(inputs[i],"_part", k, ".RDS")
		if(!file_test("-f", tname1)) stop(paste(tname1, ": No such file found in this directory!",sep = ""))
	}
}

# empty variable to store results
pcxn = c()
pb = txtProgressBar(min=0,max=max(nparts),initial=0,style=3)


# read each part of each 
for (m in 1:length(inputs)) {
	conc_tmp <- c()
	for (k in nparts[m]) {
		tmp <- readRDS(paste0(inputs[m] ,"_part", k, ".RDS"))
    	conc_tmp <- rbind(conc_tmp,tmp)
	}
	conc_tmp <- unique(conc_tmp)
    pcxn = rbind(pcxn, conc_tmp)
    setTxtProgressBar(pb,k)
}

close(pb)
pcxn <- unique(pcxn)

# adjust p-values for multiple comparison
pcxn$p.Adjust = p.adjust(p = pcxn$p.value, method = "fdr")

# save results
saveRDS(pcxn, paste0(output_name,".RDS"))
