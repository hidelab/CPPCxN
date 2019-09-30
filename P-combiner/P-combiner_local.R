p_combiner <- function (inputs,nparts,output_name) {
    # We need estimates02 to be completed for both runs
    
    options(stringsAsFactors = F)

    # ==== example inputs ====
    # inputs <- c("pcxn_conc_base","pcxn_conc_plus_10","pcxn_conc_plus_20")
    # nparts <- c(2,2,2)
    # output_name <- "pcxn_conc_base_plus_10_plus_20"
    
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
        print(paste0("adding ", inputs[m]))
        conc_tmp <- c()
        for (k in 1:nparts[m]) {
            print(k)
            tmp <- readRDS(paste0(inputs[m] ,"_part", k, ".RDS"))
            conc_tmp <- rbind(conc_tmp,tmp)
        }
        # conc_tmp <- unique(conc_tmp)
        pcxn = rbind(pcxn, conc_tmp)
        setTxtProgressBar(pb,k)
    }
    
    close(pb)
    
    # pcxn <- unique(pcxn)
    
    # adjust p-values for multiple comparison
    pcxn$p.Adjust = p.adjust(p = pcxn$p.value, method = "fdr")
    
    # save results
    saveRDS(pcxn, paste0(output_name,".RDS"))
}




