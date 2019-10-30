## A few commands for how to read json files into R

# Set working directory
setwd("C:/Users/kewo0002/gitstuff/motif_role_stability/data/Julia")

# load jsonlite
library(jsonlite)

# read in files
sim_ext_0 <- fromJSON("fw_dstrbnc/fw_1/dstrbnc_0.json")
sim_ext_1 <- fromJSON("fw_dstrbnc/fw_1/dstrbnc_1.json")
sim_ext_2 <- fromJSON("fw_dstrbnc/fw_1/dstrbnc_2.json")
sim_ext_3 <- fromJSON("fw_dstrbnc/fw_1/dstrbnc_3.json")
sim_ext_4 <- fromJSON("fw_dstrbnc/fw_1/dstrbnc_4.json")
sim_ext_5 <- fromJSON("fw_dstrbnc/fw_1/dstrbnc_5.json")

# plot biomasses
plot(sim_ext_0$B[1,],type = "l",ylim = c(0,max(sim_ext_0$B)))
for(sp in 2:20){
  lines(sim_ext_0$B[sp,],col = sp)
}


# Tip:
## sim_ext_1$B has all the biomasses (rows as species, columns as time steps)
## sim_ext_1$p has all the parameters used in the simulation
## sim_ext_1$t has time steps
