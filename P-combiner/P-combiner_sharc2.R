# Aggregate results into a single data frame
# We need estimates02 to be completed for both matrices

# How to use
# Place parts1 and 2 (output of estimates 2) of both matrices in the same folder as this script.
# As the first argument use the matrix names (the part before "_part1.RDS") seperated by commas. 
# As the second argument use the entire name of the combined matrix (e.g. matri1_matrix2.RDS)
# Example (files already in folder):
# Rscript P-combiner_sharc.R pcxn_conc_base.RDS,pcxn_conc_plus_10.RDS,pcxn_conc_plus_20.RDS all_concs

rm(list = ls(all=TRUE)) 
options(stringsAsFactors = F)

# ==== Arguments ====
# command line arguments from submission (sharc) script
# 1 - genesets file
cmd_args <- commandArgs(trailingOnly = T)

# ==== INPUTS ====
inputs <- unlist(strsplit(cmd_args[1], ","))
#inputs <- c("pcxn_conc_base","pcxn_conc_plus_10","pcxn_conc_plus_20")

output_name <- cmd_args[2]

# empty variable to store results
pcxn = c()
conc_tmp <- c()
pb = txtProgressBar(min=0,max=length(inputs),initial=0,style=3)

# read each part of each 
for (m in 1:length(inputs)) {
    print(paste0("adding ", inputs[m]))
    tmp <- readRDS(paste0(inputs[m]))
    print(nrow(tmp))
    # conc_tmp <- unique(conc_tmp)
    
    #test for sameness here
    #order by drug? 
    # same <- tmp[conc_tmp$Pathway.A == conc_tmp]
    
    conc_tmp <- rbind(conc_tmp,tmp)
    print(nrow(conc_tmp))
    setTxtProgressBar(pb,m)
    rm(tmp)
}

pcxn = conc_tmp
rm(conc_tmp)

close(pb)
print(nrow(pcxn))
print("about to remove duplicates")

pcxn$p.Adjust <- NULL

print(Sys.time())
pcxn <- dplyr::distinct(pcxn, Pathway.A, Pathway.B, PathCor, Overlap.Coefficient, .keep_all = T)
print(Sys.time())

# print(Sys.time())
# pcxn <- pcxn[!duplicated(pcxn[,c(1:4)] ),]
# print(Sys.time())


print(nrow(pcxn))

# adjust p-values for multiple comparison
pcxn$p.Adjust = p.adjust(p = pcxn$p.value, method = "fdr")

# save results
saveRDS(pcxn[pcxn$p.Adjust < 1, ], paste0(output_name,"p.adj<1.RDS"))
saveRDS(pcxn, paste0(output_name,".RDS"))
print("done")
