p_combiner <- function (inputs,output_name) {
    # We need estimates02 to be completed for both runs
    
    options(stringsAsFactors = F)
    
    # ==== example inputs ====
    # inputs <- c("path_to_pdxn_1","path_to_pdxn_2","path_to_pdxn_3")
    # output_name <- "pcxn_conc_base_plus_10_plus_20"
    
    
    # empty variable to store results
    pcxn = c()
    conc_tmp <- c()
    pb = txtProgressBar(min=0,max=max(length(inputs)),initial=0,style=3)
    
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
    }
    
    # xxx <- conc_tmp[grepl("SKMEL1", conc_tmp$Pathway.A),]
    # table(duplicated(xxx[,1:4]))
    # nrow(conc_tmp)
    pcxn = conc_tmp
    
    close(pb)
    print(nrow(pcxn))
    print("about to remove duplicates")
    print(Sys.time())
    pcxn <- pcxn[!duplicated(pcxn[,c(1:4)] ),]
    print(Sys.time())
    print(nrow(pcxn))
    
    # adjust p-values for multiple comparison
    pcxn$p.Adjust = p.adjust(p = pcxn$p.value, method = "fdr")
    
    # save results
    saveRDS(pcxn[pcxn$p.Adjust < 1, ], paste0(output_name,"p.adj<1.RDS"))
    saveRDS(pcxn, paste0(output_name,".RDS"))
    print("done")
}


# inputs <- c("/Users/kkoler/Google Drive/Documents/new_PDN/miniPDN/improved_PCxN_MSigDB.L1000CDS2.staticmod.6.subset.RDS", "/Users/kkoler/Google Drive/Documents/new_PDN/miniPDN/improved_PCxN_MSigDB.L1000CDS2.staticmod.24.subset.RDS")
# 
# output_name <- c("test_combined")
# 
# p_combiner(inputs, output_name)
