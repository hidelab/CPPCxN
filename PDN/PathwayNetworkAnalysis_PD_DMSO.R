# Network based analysis of alzheimers pathways from mt.saina data from Sandeep
# Version: Alzhemier's Disease
# Init date: 20171220

# Guide
# 1. Place the pathway data file in the same dir as this script

#### Loading libraries ####
library(readxl)
library(pathprint)

#### Inputs ####
# Project files (change those lines if you are starting a new project)
project_outer_dir <- "Parkinsons"
project_inner_folder <- "patient_control_DMSO" 
pathwayDataFilename <- "pathway_Fibroblast_RNAseq_Patient_DMSO_vs_Control_DMSO_limma sigpathways 0.05p_1.2FC_modified.xls"
humanFilename <-"square_data_50000.RDS"

# Set desired column names
pathway_col <- "pathway_id"
foldChange_col <- "logFC"

#### Loading resources ####

# Creating the project directory if needed(new project)
if (!file.exists(project_outer_dir)) {
    dir.create(project_outer_dir, showWarnings = FALSE)
    dir.create(paste(project_outer_dir, "/", project_inner_folder, sep = ""), showWarnings = FALSE)
}

# load data from file
# extract significant pathways based on cutoffs
pathwayData <- read_excel(paste("data/", pathwayDataFilename, sep = ""))

setwd(paste(project_outer_dir, "/", project_inner_folder, "/", sep = ""))

# Output dir setup
OutputDataDir <- paste("/", project_outer_dir, "/", project_inner_folder, sep = "")

#### Already ordered but re-do just to check #### 
# convert list to vector to be able to use order() and have a variable column name
lname <- pathwayData[,foldChange_col]
vname <- unlist(lname, use.names=FALSE)
pathwayData <- pathwayData[order(vname, decreasing = T),]

#### NEEDS TO ADJUST FOR THE CASE THERE ARE LESS THAN 5 pathways #### 
# Cluster A) target top 5 pathways
ClusterA <- head(pathwayData[,pathway_col],5)
# Cluster B) target top 10 pathways
ClusterB <- head(pathwayData[,pathway_col],10)
# Cluster C) target bottom 5 pathways
ClusterC <- tail(pathwayData[,pathway_col],5)
# Cluster D) target bottom 10 pathways
ClusterD <- tail(pathwayData[,pathway_col],10)

print("ClusterA")
ClusterA
print("ClusterB")
ClusterB
print("ClusterC")
ClusterC
print("ClusterD")
ClusterD

#### Loading human data from pathprint ####
data(pathprint.Hs.gs)

#### Check that all pathways exist and are named correctly #### 
ClusterA[!(ClusterA %in% names(pathprint.Hs.gs))]
ClusterB[!(ClusterB %in% names(pathprint.Hs.gs))]
ClusterC[!(ClusterC %in% names(pathprint.Hs.gs))]
ClusterD[!(ClusterD %in% names(pathprint.Hs.gs))]

length(ClusterA) == sum(names(pathprint.Hs.gs) %in% ClusterA)
length(ClusterB) == sum(names(pathprint.Hs.gs) %in% ClusterB)
length(ClusterC) == sum(names(pathprint.Hs.gs) %in% ClusterC)
length(ClusterD) == sum(names(pathprint.Hs.gs) %in% ClusterD)

#### Strategy ####
# 1) extract a network with links to at least 4 nodes from each cluster.
# 2) Find subnetworks within these connected clusters

#### Load human data and replacing NAs with zeros ####
cor.matrix = readRDS(paste("../../data/", humanFilename, sep = ""))
cor.matrix[is.na(cor.matrix)] <- 0

#### Just use highest concentration of CMPA drug ####	 
CMAPnames <- rownames(cor.matrix)[grepl("CMAP.", rownames(cor.matrix))]
CMAPdrug <- sapply(CMAPnames, function(x){unlist(strsplit(x, "_"))[2]})
CMAPconcentration <- sapply(CMAPnames, function(x){unlist(strsplit(x, "_"))[3]})
CMAPcell <- sapply(CMAPnames, function(x){unlist(strsplit(x, "_"))[4]})
CMAPdirection <- sapply(CMAPnames, function(x){unlist(strsplit(x, "_"))[5]})

#### Define CMAP set using the highest drug concentration ####
CMAPnamesUse <- c()
for (i in unique(CMAPdirection)){
    for (j in unique(CMAPcell[CMAPdirection == i])){
        for (k in unique(CMAPdrug[CMAPdirection == i & CMAPcell == j])){
            CMAPnamesSub <- CMAPnames[CMAPdirection == i & CMAPcell == j & CMAPdrug == k]
            CMAPnamesUse <- c(CMAPnamesUse, CMAPnamesSub[which.max(CMAPconcentration[CMAPdirection == i & CMAPcell == j & CMAPdrug == k])])
        }
    }
}

#### NEEDS TO RELOCATE OR SPLIT. Create CMAP network with all Pathways ####
library(GeneNet)
tempm <- cor.matrix[
    (grepl("CTD.disease.", rownames(cor.matrix)) |
         grepl("CTD.chem.", rownames(cor.matrix)) |
         grepl("PharmGKB.", rownames(cor.matrix)) |
         grepl("Pathway.", rownames(cor.matrix)) |
         grepl("Hide.", rownames(cor.matrix)) |
         rownames(cor.matrix) %in% CMAPnamesUse
    ),
    (grepl("CTD.disease.", rownames(cor.matrix)) |
         grepl("CTD.chem.", rownames(cor.matrix)) |
         grepl("PharmGKB.", rownames(cor.matrix)) |
         grepl("Pathway.", rownames(cor.matrix)) |
         grepl("Hide.", rownames(cor.matrix)) |
         rownames(cor.matrix) %in% CMAPnamesUse
    )
]

class(tempm) <- "numeric"

DPD.pcor.est <- ggm.estimate.pcor(tempm)

save(DPD.pcor.est, file = paste("DPD.pcor.est.RData", sep = "/"))

#### Loading GeneNet object and create network ####
load("DPD.pcor.est.RData")
DPD.pcor.est.results <- network.test.edges(DPD.pcor.est, direct = F, plot=F)

nodeNames <- rownames(DPD.pcor.est)
testNode <- 1:length(nodeNames)
CMAPtestnodes <- match(CMAPnamesUse, nodeNames)

#### Create p-value for positive correlation ####
DPD.pcor.est.results$pValPos <- NULL
DPD.pcor.est.results$pValPos[DPD.pcor.est.results$pcor >= 0] <- DPD.pcor.est.results$pval[DPD.pcor.est.results$pcor >= 0]
#### If correlation is negative, set the p-value for +ve correlation to 1 ####
DPD.pcor.est.results$pValPos[DPD.pcor.est.results$pcor < 0] <- 1

#### Create p-value for negative correlation ####
DPD.pcor.est.results$pValNeg <- NULL
DPD.pcor.est.results$pValNeg[DPD.pcor.est.results$pcor <= 0] <- DPD.pcor.est.results$pval[DPD.pcor.est.results$pcor <= 0]
#### If correlation is positive, set the p-value for +ve correlation to 1 ####
DPD.pcor.est.results$pValNeg[DPD.pcor.est.results$pcor > 0] <- 1

#### Ensure that CMAP and target reference is consistent whether node1 or node2	####
indx1a <- DPD.pcor.est.results$node1 %in% CMAPtestnodes & !(DPD.pcor.est.results$node2 %in% CMAPtestnodes)
indx2a <- DPD.pcor.est.results$node2 %in% CMAPtestnodes & !(DPD.pcor.est.results$node1 %in% CMAPtestnodes)
indx3a <- DPD.pcor.est.results$node1 %in% CMAPtestnodes & DPD.pcor.est.results$node2 %in% CMAPtestnodes

indx1b <- DPD.pcor.est.results$node1 %in% testNode & !(DPD.pcor.est.results$node1 %in% CMAPtestnodes)
indx2b <- DPD.pcor.est.results$node2 %in% CMAPtestnodes & !(DPD.pcor.est.results$node1 %in% CMAPtestnodes)
indx3b <- DPD.pcor.est.results$node1 %in% CMAPtestnodes & DPD.pcor.est.results$node2 %in% CMAPtestnodes

DPD.pcor.est.results$CMAPnode <- NULL
DPD.pcor.est.results$testNode <- NULL

DPD.pcor.est.results$CMAPnode[indx1a] <- DPD.pcor.est.results$node1[indx1a]
DPD.pcor.est.results$CMAPnode[indx2a] <- DPD.pcor.est.results$node2[indx2a]
DPD.pcor.est.results$testNode[indx1a] <- DPD.pcor.est.results$node2[indx1a]
DPD.pcor.est.results$testNode[indx2a] <- DPD.pcor.est.results$node1[indx2a]

#### For CMAP-CMAP nodes just use node1 and node2 ####
#DPD.pcor.est.results$testNode[indx3a] <- DPD.pcor.est.results$node1[indx3a]
#DPD.pcor.est.results$CMAPnode[indx3a] <- DPD.pcor.est.results$node2[indx3a]

save(DPD.pcor.est.results, file = "DPD.pcor.est.results.RData")
# load("DPD.pcor.est.results.RData)

#### Want to know ####
# 1 how pathways relate to drugs
# 2 how pathways relate to each other
# 3 how drugs relate to each other

testClusters <- list(ClusterA = paste("Pathway", ClusterA, sep ="."),
                     ClusterB = paste("Pathway", ClusterB, sep ="."),
                     ClusterC = paste("Pathway", ClusterC, sep ="."),
                     ClusterD = paste("Pathway", ClusterD, sep ="."))


#### loop ####
if (2 == 2){
    # extract network at high cutoff
    DPD.pcor.est.net.0.9 <- extract.network(DPD.pcor.est.results, cutoff.ggm=0.9)
    
    # Now want to extract nodes that are connected to 4 out of 5 of the cluster nodes
    nodeNames <- rownames(DPD.pcor.est)
    
    testNetworkList <- vector("list", length = length(testClusters)) 
    names(testNetworkList) <- names(testClusters)
    for (i in 1:length(testClusters))
    {
        testNodes <- match(testClusters[[i]], nodeNames)
        includeNodeNames <- testClusters[[i]]
        # Set threshold number of nodes in the cluster that must be connected for a node to be included in the network
        threshold <- (length(testNodes) - 1)
        scoreData = as.data.frame(matrix(nrow = length(nodeNames),ncol = 3))
        rownames(scoreData) <- nodeNames
        colnames(scoreData) <- c("Score", "SumPos", "SumNeg")
        tempTable <- DPD.pcor.est.net.0.9[DPD.pcor.est.net.0.9$node1 %in% testNodes | DPD.pcor.est.net.0.9$node2 %in% testNodes,]
        for (j in 1:length(nodeNames))
        {
            
            # Calculate the number of edges with the testCluster group
            score <- sum(
                (tempTable$node1 %in% testNodes & tempTable$node2 == j) |
                    (tempTable$node2 %in% testNodes & tempTable$node1 == j)
            )
            sumPos <- sum(tempTable$pcor[
                (tempTable$node1 %in% testNodes & tempTable$node2 == j) |
                    (tempTable$node2 %in% testNodes & tempTable$node1 == j)
                ] > 0)
            sumNeg <- sum(tempTable$pcor[
                (tempTable$node1 %in% testNodes & tempTable$node2 == j) |
                    (tempTable$node2 %in% testNodes & tempTable$node1 == j)
                ] < 0)											
            #print(j)
            #print(score)
            scoreData[j,1] <- score
            scoreData[j,2] <- sumPos
            scoreData[j,3] <- sumNeg
            # if (score >= threshold)
            # {
            # # if connected to at least threshold members then add to network
            # includeNodeNames <- union(includeNodeNames, nodeNames[j])
            # }
        }
        includeNodeNames <- union(includeNodeNames, rownames(scoreData)[scoreData[,1] >= 3])	
        includeNodeIndx = match(includeNodeNames, nodeNames)
        # Create network that uses only the included nodes	
        testNetwork.0.9 <- DPD.pcor.est.net.0.9[
            DPD.pcor.est.net.0.9$node1 %in% includeNodeIndx & DPD.pcor.est.net.0.9$node2 %in% includeNodeIndx,
            ]
        testNetwork.0.9 <- testNetwork.0.9[order(testNetwork.0.9$prob,decreasing = F),]
        #  	testNetwork <- head(testNetwork,20)
        testNetwork.0.9$node1Name = nodeNames[testNetwork.0.9$node1]
        testNetwork.0.9$node2Name = nodeNames[testNetwork.0.9$node2]
        
        # assign to list
        testNetworkList[[i]] <- list(high = testNetwork.0.9, scoreData = scoreData) 
    }
    
    # Save testNetwork list
    testNetworkListFile <- "PD.DMSO.Pathway.human.network.RData"
    save(testNetworkList, file = testNetworkListFile)
    
    for (i in 1:length(testClusters))
    { 
        testNetwork.0.9 <- testNetworkList[[i]]$high
        scoreData <- testNetworkList[[i]]$scoreData
        print(head(scoreData))
        
        # Produce networks with annotations
        testNetwork.0.9$NodeType1 <- NA
        testNetwork.0.9$NodeType1[grep("CMAP.up", testNetwork.0.9$node1Name)] <- "CMAP_UP"
        testNetwork.0.9$NodeType1[grep("CMAP.down", testNetwork.0.9$node1Name)] <- "CMAP_DOWN"
        testNetwork.0.9$NodeType1[grep("Pathway.", testNetwork.0.9$node1Name)] <- "Pathway"
        testNetwork.0.9$NodeType1[grep("PharmGKB.disease", testNetwork.0.9$node1Name)] <- "PharmGKB_Disease"
        testNetwork.0.9$NodeType1[grep("PharmGKB.drug", testNetwork.0.9$node1Name)] <- "PharmGKB_drug"
        testNetwork.0.9$NodeType1[grep("CTD.disease", testNetwork.0.9$node1Name)] <- "CTD_Disease"
        testNetwork.0.9$NodeType1[grep("CTD.chem", testNetwork.0.9$node1Name)] <- "CTD_Chem"
        
        testNetwork.0.9$NodeType2 <- NA
        testNetwork.0.9$NodeType2[grep("CMAP.up", testNetwork.0.9$node2Name)] <- "CMAP_UP"
        testNetwork.0.9$NodeType2[grep("CMAP.down", testNetwork.0.9$node2Name)] <- "CMAP_DOWN"
        testNetwork.0.9$NodeType2[grep("Pathway.", testNetwork.0.9$node2Name)] <- "Pathway"
        testNetwork.0.9$NodeType2[grep("PharmGKB.disease", testNetwork.0.9$node2Name)] <- "PharmGKB_Disease"
        testNetwork.0.9$NodeType2[grep("PharmGKB.drug", testNetwork.0.9$node2Name)] <- "PharmGKB_drug"
        testNetwork.0.9$NodeType2[grep("CTD.disease", testNetwork.0.9$node2Name)] <- "CTD_Disease"
        testNetwork.0.9$NodeType2[grep("CTD.chem", testNetwork.0.9$node2Name)] <- "CTD_Chem"
        
        # Remove commas for output
        testNetwork.0.9.Output <- testNetwork.0.9
        testNetwork.0.9.Output$node1Name <- gsub(",", "", testNetwork.0.9.Output$node1Name)
        testNetwork.0.9.Output$node2Name <- gsub(",", "", testNetwork.0.9.Output$node2Name)
        rownames(scoreData) <- gsub(",", "", rownames(scoreData))
        
        # Save network
        networkOutputFile <- paste("PDDMSOPathwayNetwork", names(testClusters)[i], "csv", sep = ".")
        write.csv(testNetwork.0.9.Output[,c(12,13,14,15,1,4,5,6)],
                  file = networkOutputFile, quote = F)
        
        # Save attributes
        attribute.table <- data.frame(NodeName = unique(c(testNetwork.0.9.Output$node1Name,
                                                          testNetwork.0.9.Output$node2Name)))
        attribute.table$Score <- scoreData$Score[match(attribute.table$NodeName, rownames(scoreData))]
        attribute.table$Score[
            attribute.table$NodeName %in% gsub(",", "", testClusters[[i]])] <- 99
        attribute.table$SumPos <- scoreData$SumPos[match(attribute.table$NodeName, rownames(scoreData))]		
        attribute.table$SumNeg <- scoreData$SumNeg[match(attribute.table$NodeName, rownames(scoreData))]			
        
        attribute.table$NodeType <- NA
        attribute.table$NodeType[grep("CMAP.up", attribute.table$NodeName)] <- "CMAP_UP"
        attribute.table$NodeType[grep("CMAP.down", attribute.table$NodeName)] <- "CMAP_DOWN"
        attribute.table$NodeType[grep("Pathway.", attribute.table$NodeName)] <- "Pathway"
        attribute.table$NodeType[grep("PharmGKB.disease", attribute.table$NodeName)] <- "PharmGKB_Disease"
        attribute.table$NodeType[grep("PharmGKB.drug", attribute.table$NodeName)] <- "PharmGKB_drug"
        attribute.table$NodeType[grep("CTD.disease", attribute.table$NodeName)] <- "CTD_Disease"
        attribute.table$NodeType[grep("CTD.chem", attribute.table$NodeName)] <- "CTD_Chem"
        
        networkAttributesOutputFile <- paste("ADPathwayNetworkAttributes", names(testClusters)[i], "csv", sep = ".")
        write.csv(attribute.table, file = networkAttributesOutputFile, quote = F)
    }
}

#### ??ry to create composite CMAP score by combining p-value for positive and negative correlation ####
# based on ideas proposed here: http://www.stat.wisc.edu/~wardrop/courses/meta2.pdf
# assume that the p-values have a uniform distribution under the null hypothesis
# sum of -natural log of p-values has chi-squared distribution

#### ??eed to split into up and down ####
CMAPtestnodes <- match(CMAPnamesUse, nodeNames)
clusterTestnames <- unlist(testClusters)
clusterTestnodes <- match(unlist(testClusters), nodeNames)
CMAP.results <- DPD.pcor.est.results[(DPD.pcor.est.results$node1 %in% CMAPtestnodes & DPD.pcor.est.results$node2 %in% clusterTestnodes) | 
                                         (DPD.pcor.est.results$node2 %in% CMAPtestnodes & DPD.pcor.est.results$node1 %in% clusterTestnodes),
                                     ]

#### Create p-value for positive correlation ####
CMAP.results$pValPos <- NULL
CMAP.results$pValPos[CMAP.results$pcor >= 0] <- CMAP.results$pval[CMAP.results$pcor >= 0]
##### If correlation is negative, set the p-value for +ve correlation to 1 ####
# CMAP.results$pValPos[CMAP.results$pcor < 0] <- (1-CMAP.results$pval[CMAP.results$pcor < 0])
CMAP.results$pValPos[CMAP.results$pcor < 0] <- 1

#### Create p-value for negative correlation ####
CMAP.results$pValNeg <- NULL
CMAP.results$pValNeg[CMAP.results$pcor <= 0] <- CMAP.results$pval[CMAP.results$pcor <= 0]
#### If correlation is positive, set the p-value for +ve correlation to 1 ####
# CMAP.results$pValNeg[CMAP.results$pcor > 0] <- (1-CMAP.results$pval[CMAP.results$pcor >= 0])
CMAP.results$pValNeg[CMAP.results$pcor > 0] <- 1

#### Now combine pos and negative p-values across CMAP_UP and CMAP_DOWN ####
#### Create CMAP lookup table ####
# This table contains p-vales for CMap negative or positive association combined across the CMAP_UP and CMAP_DOWN signatures.
CMAPlookup <- as.data.frame(matrix(ncol = 3, nrow = length(CMAPnamesUse)/2))
colnames(CMAPlookup) <- c("UPindx", "DOWNindx", "MasterIndx")
rownames(CMAPlookup) <-gsub("CMAP.up.", "", gsub("_up", "", CMAPnamesUse[grep("CMAP.up", CMAPnamesUse)]))
CMAPlookup$MasterIndx <- 1:nrow(CMAPlookup)
CMAPlookup$UPindx <- match(paste("CMAP.up.", paste(rownames(CMAPlookup), "_up", sep = ""), sep = ""), nodeNames)
CMAPlookup$DOWNindx <- match(paste("CMAP.down.", paste(rownames(CMAPlookup), "_down", sep = ""), sep = ""), nodeNames)

CMAPlookup <- cbind(CMAPlookup, matrix(nrow = nrow(CMAPlookup), ncol = length(clusterTestnodes)*2))
colnames(CMAPlookup)[4:(3+length(clusterTestnames))]<-paste(clusterTestnames, "pos", sep = "_")
colnames(CMAPlookup)[(4+length(clusterTestnames)):(3+length(clusterTestnames)*2)]<-paste(clusterTestnames, "neg", sep = "_")

for (i in 1:length(clusterTestnodes)){
    temp.Table <- CMAP.results[CMAP.results$testNode %in% clusterTestnodes[i],]
    temp.upPpos <- temp.Table$pValPos[match(CMAPlookup$UPindx, temp.Table$CMAPnode)]
    temp.downPneg <- temp.Table$pValNeg[match(CMAPlookup$DOWNindx, temp.Table$CMAPnode)]
    temp.combinedPpos <- pchisq(-log(temp.upPpos) - log(temp.downPneg),2, lower.tail = F)	
    temp.upPneg <- temp.Table$pValNeg[match(CMAPlookup$UPindx, temp.Table$CMAPnode)]
    temp.downPpos <- temp.Table$pValPos[match(CMAPlookup$DOWNindx, temp.Table$CMAPnode)]
    temp.combinedPneg <- pchisq(-log(temp.upPneg) - log(temp.downPpos),2, lower.tail = F)	
    CMAPlookup[,(3+i)] <- temp.combinedPpos
    CMAPlookup[,(length(clusterTestnames)+3+i)] <- temp.combinedPneg
    print(i)
}

#### Examine sample node ####
# combined
CMAPlookup[grep("nicardipine_7.80e-06_MCF7", rownames(CMAPlookup)),
           colnames(CMAPlookup) %in% paste(testClusters[[4]],"_neg",sep = "")]
CMAPlookup[grep("nicardipine_7.80e-06_MCF7", rownames(CMAPlookup)),
           colnames(CMAPlookup) %in% paste(testClusters[[4]],"_pos",sep = "")]

# original
temp1 <- CMAP.results[CMAP.results$node2 %in% grep("nicardipine_7.80e-06_MCF7",nodeNames) & 
                          CMAP.results$node1 %in% match(unlist(testClusters[[4]]), nodeNames),]
temp1$node2Name <- nodeNames[temp1$node2]
temp1$node1Name <- nodeNames[temp1$node1]

# Now combine across clusters
ClusterCMAPpvals <- as.data.frame(matrix(nrow = nrow(CMAPlookup), ncol = 2*length(testClusters)))
rownames(ClusterCMAPpvals) <- rownames(CMAPlookup)
colnames(ClusterCMAPpvals) <- c(paste(names(testClusters), ".Pos", sep = ""),
                                paste(names(testClusters), ".Neg", sep = ""))

for (i in 1:length(testClusters))
{ 
    CMAPlookupSub_pos <- CMAPlookup[,colnames(CMAPlookup) %in% paste(testClusters[[i]],"_pos",sep = "")]
    CMAPlookupSub_neg <- CMAPlookup[,colnames(CMAPlookup) %in% paste(testClusters[[i]],"_neg",sep = "")]
    ClusterCMAPpvals[,i] <- apply(CMAPlookupSub_pos,1,function(x){
        pchisq(sum(-log(x)),length(x),lower.tail = F)
    })
    ClusterCMAPpvals[,(i+length(testClusters))] <- apply(CMAPlookupSub_neg,1,function(x){
        pchisq(sum(-log(x)),length(x),lower.tail = F)
    })
}

#### Create single combined score from the 3 p-values ####
ClusterCMAPpvals$AposBnegCpos <- pchisq( (-log(ClusterCMAPpvals$ClusterA.Pos)) + 
                                             (-log(ClusterCMAPpvals$ClusterB.Neg)) +
                                             (-log(ClusterCMAPpvals$ClusterC.Pos)),3,lower.tail = F)
ClusterCMAPpvals$AnegBposCneg <- pchisq( (-log(ClusterCMAPpvals$ClusterA.Neg)) + 
                                             (-log(ClusterCMAPpvals$ClusterB.Pos)) +
                                             (-log(ClusterCMAPpvals$ClusterC.Neg)),3,lower.tail = F)	

# now order by negative correlation to each of the pathways
#write.csv(ClusterCMAPpvals[order(ClusterCMAPpvals$AnegBposCneg, decreasing = F),], file = "ClusterCMAPpvals.csv")
#ClusterCMAPpvals[grep("fulvestrant", rownames(ClusterCMAPpvals)),]
#ClusterCMAPpvals <- read.csv("ClusterCMAPpvals.csv", row.names = 1, stringsAsFactors = F)

#### Combine p-values to create scored output by rank ####
ClusterCMAPScores <- as.data.frame(matrix(nrow = nrow(ClusterCMAPpvals), ncol = 3*length(testClusters)))
rownames(ClusterCMAPScores) <- rownames(ClusterCMAPpvals)
colnames(ClusterCMAPScores) <- c(paste(names(testClusters), ".CombinedScore", sep = ""),
                                 paste(names(testClusters), ".PosRank", sep = ""),
                                 paste(names(testClusters), ".NegRank", sep = ""))

for (i in 1:length(testClusters))
{ 
    # create rank for positive and negative association
    tempPos <- ClusterCMAPpvals[,paste(names(testClusters)[i], ".Pos", sep = "")]
    tempNeg <- ClusterCMAPpvals[,paste(names(testClusters)[i], ".Neg", sep = "")] 
    ClusterCMAPScores[,paste(names(testClusters)[i], ".PosRank", sep = "")] <- rank(tempPos)
    ClusterCMAPScores[,paste(names(testClusters)[i], ".NegRank", sep = "")] <- rank(tempNeg)
    # combine positive and negative
    tempCombined <-  rank(rank(tempNeg) - rank(tempPos))/(nrow(ClusterCMAPScores)/2) - 1
    tempCombined[tempPos > 0.1 & tempNeg > 0.1] <- 0
    ClusterCMAPScores[,paste(names(testClusters)[i], ".CombinedScore", sep = "")] <- tempCombined
}  

ClusterCMAPpvals <- cbind(ClusterCMAPScores, ClusterCMAPpvals)
write.csv(ClusterCMAPpvals[order(ClusterCMAPpvals[,paste(names(testClusters)[1], ".CombinedScore", sep = "")], decreasing = F),],
          file = "ClusterCMAPpvals.csv")
