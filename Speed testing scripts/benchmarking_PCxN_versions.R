library(microbenchmark)

source("og_01.R")
source("og_subpairs_01.R")
source("JD_01.R")
source("JD_subpairs_01.R")

# # Checking identical results
# check_identical_output <- function(a,b,c,d) {
#     
# }

lower <- 600
upper <- 699

mbm <- microbenchmark("original" = { og_01(lower,upper) },
                      "original_subpairs" = { og_subpairs_01(lower,upper) },
                      "joint_disjoint" = { JD_01(lower,upper) },
                      "joint_disjoint_subpairs" = { JD_subpairs_01(lower,upper) },
                      times = 10)

original_performance <-subset(mbm, mbm$expr == 'original')
original_sub_performance <-subset(mbm, mbm$expr == 'original_subpairs')
joint_disjoint_performance <-subset(mbm, mbm$expr == 'joint_disjoint')
joint_disjoint_sub_performance <-subset(mbm, mbm$expr == 'joint_disjoint_subpairs')