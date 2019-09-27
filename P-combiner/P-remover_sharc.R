# Remove genesets from tissue level correlations
# We need estimates01 to be completed for all parts

# How to use
# Specify path to output of estimates 1
# Specify a path to where the finished product should go, same folder structure as before
# Specfy name of genesets to be removed. 
# 
# Example (files already in folder):
# Rscript P-remover_sharc.R ../output_v2_rel7/ ../output_v2_rel7_core/ seed_drugs.RDS 1:134

rm(list = ls(all=TRUE)) 
options(stringsAsFactors = F)

# remove_genesets_name <- "seed_drugs.RDS"
# remove_gs <- readRDS(paste0("/Users/kkoler/Google Drive/Documents/new_PDN/LINCSCDgenesets/", remove_genesets_name))

# ==== Arguments ====
# command line arguments from submission (sharc) script
# 1 - genesets file
cmd_args <- commandArgs(trailingOnly = T)

# ==== INPUTS ====
dir_to_be_changed <- cmd_args[1]
#inputs <- c("pcxn_conc_base","pcxn_conc_plus_10","pcxn_conc_plus_20")
dir_changed <- cmd_args[2]

remove_genesets_name <- cmd_args[3]

range_to_change <- unlist(strsplit(cmd_args[4], ":"))



remove_gs <- readRDS(paste0("../data/", remove_genesets_name))
remove_gs_names <- names(remove_gs)

files_to_change <- list.files(path=paste0(dir_to_be_changed, "mean_pcor2_barcode"))
pb = txtProgressBar(min=0,max=length(files_to_change),initial=0,style=3)

# Read in file and subset it, then save
for(i in files_to_change[range_to_change[1]:range_to_change[2]]) {
  print(Sys.time())
  path_to_file <- paste0(dir_to_be_changed, "mean_pcor2_barcode/", i)
  print(i)
  print(path_to_file)
  file_to_change <- readRDS(path_to_file)
  orig_size <- length(file_to_change)
  file_to_change <- file_to_change[which(!sapply(file_to_change, `[[`, 1) %in% remove_gs_names)]
  new_path_to_file <- sub(dir_to_be_changed, dir_changed, path_to_file)
  new_size <- length(file_to_change)
  print(paste0("saving to ", new_path_to_file))
  saveRDS(file_to_change, new_path_to_file)
  
  rm(file_to_change)
  print(paste0("done ", which(files_to_change %in% i)))
  setTxtProgressBar(pb,which(files_to_change %in% i))
}

close(pb)

print(paste0(new_size/orig_size*100, "% of the old size"))

accepted_pairs <- new_size
actual_pairs <- data.frame("actual_number_of_pairs" = accepted_pairs)
saveRDS(actual_pairs, paste0(dir_changed,"/mean_pcor2_barcode/res/pairs.RDS"))

print("done all of the files")