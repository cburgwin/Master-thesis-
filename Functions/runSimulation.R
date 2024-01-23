#### Libraries
library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx) # to save csv file
library(ggplot2) # to plot
library(RColorBrewer) # colour palette
library(future.apply)
library(parallel)
library(emmeans)

rm(list=ls())

# use source function to use functions that are created in another R script
source("callSimulation.R")

runSimulation <- function(t, # number of time periods
                          n, # number of patients for each time period
                          armsNumb_vec = 1:3, # number of arms
                          prevalence, # allocation probability to the different types of groups
                          deltaT1, # treatment effect for treatment 1
                          deltaT2, # treatment effect for treatment 2
                          meanCX, # mean for control in group X
                          meanCY, # mean for control in group Y
                          meanCZ, # mean for control in group Z
                          sdT1X, # standard deviation for treatment 1 in group X
                          sdT2Y, # standard deviation for treatment 2 in group Y
                          sdT1Z, # standard deviation for treatment 1 in group Z
                          sdT2Z, # standard deviation for treatment 2 in group Z
                          sdCX, # standard deviation for control in group X
                          sdCY, # standard deviation for control in group Y
                          sdCZ, # standard deviation for control in group 
                          blocksizeXY, # block size for block randomization for groups X and Y
                          blocksizeZ, # block size for block randomization for group Z
                          alpha, # alpha level for t-test
                          allocProbRandXY_vec, # allocation probability for groups X & Y to treatment or control
                          allocProbRandZ_vec, # allocation probability for group Z to different arms
                          complete, # boolean for randomization method
                          random, # boolean for deterministic or random allocation to groups
                          full,  # boolean for returning input parameters or just the output
                          iterations # number of iterations
) {
  
  # calculate the sample size per group
  n = (n * 3) / t
  
  # determine the grid over all combinations
  grid1_mat <- expand.grid("Treatment effect T1" = deltaT1,
                           "Treatment effect T2" = deltaT2, 
                           "Mean Control X" = meanCX,
                           "Mean Control Y" = meanCY,
                           "Mean Control Z" = meanCZ,
                           "Std Control X" = sdCX,
                           "Std Control Y" = sdCY,
                           "Std Control Z" = sdCZ,
                           "Std T1 X" = sdT1X,
                           "Std T1 Z" = sdT1Z,
                           "Std T2 Y" = sdT2Y,
                           "Std T2 Z" = sdT2Z,
                           "Time periods" = t,
                           "Sample Size" = n,
                           "Prevalence" = prevalence,
                           "Alloc Prob Rand XY" = allocProbRandXY_vec,
                           "Alloc Prob Rand Z" = allocProbRandZ_vec,
                           "Blocksize XY" = blocksizeXY,
                           "Blocksize Z" = blocksizeZ,
                           "Number of Iterations" = iterations,
                           "Signficance level" = alpha,
                           "Complete Randomization" = complete,
                           "Random Allocation" = random,
                           "Full Output" = full)
  
 # only keep the combinations of the three control means that are of interest
   grid1_mat <- grid1_mat %>%
     filter((`Mean Control X` == 0 & `Mean Control Y` == 0 & `Mean Control Z` == 0) |
              (`Mean Control X` == -0.5 & `Mean Control Y` == 0.5 & `Mean Control Z` == 0) |
              (`Mean Control X` == -0.5 & `Mean Control Y` == 0 & `Mean Control Z` == 0.5))
  
  
  # only keep combinations of randomization ratios in groups which are of interest
  grid1_mat <- grid1_mat %>%
    filter(!(`Alloc Prob Rand XY` == "c(0.333333333333333, 0.333333333333333, 0.333333333333333)" & `Alloc Prob Rand Z` == "c(0.25, 0.25, 0.25, 0.25)") &
             !(`Blocksize XY` == 2 & `Alloc Prob Rand XY` == "c(0.333333333333333, 0.333333333333333, 0.333333333333333)") &
             !(`Blocksize XY` == 3 & `Alloc Prob Rand XY` == "c(0.5, 0.5)") &   
             !(`Blocksize Z` == 4 & `Alloc Prob Rand Z` == "c(0.333333333333333, 0.333333333333333, 0.333333333333333)") &
             !(`Blocksize Z` == 3 & `Alloc Prob Rand Z` == "c(0.25, 0.25, 0.25, 0.25)"))


  # only keep the unqiue combinations
  grid1_mat <- unique(grid1_mat)
  
  
  # iterate over the grid and call the function for all combinations save the output in a list
  results_list <- list()
  for(i in 1:nrow(grid1_mat)) {
    results_list[[i]] <- callSimulation(n = grid1_mat$`Sample Size`[i],
                                        allocProb_vec = grid1_mat$`Prevalence`[[i]],
                                        armsNumb_vec = 1:3,
                                        t = grid1_mat$`Time periods`[i],
                                        deltaT1 = grid1_mat$`Treatment effect T1`[i],
                                        deltaT2 = grid1_mat$`Treatment effect T2`[i], 
                                        meanCX = grid1_mat$`Mean Control X`[i], 
                                        meanCY = grid1_mat$`Mean Control Y`[i],
                                        meanCZ = grid1_mat$`Mean Control Z`[i],
                                        sdCX = grid1_mat$`Std Control X`[i], 
                                        sdCY = grid1_mat$`Std Control Y`[i], 
                                        sdCZ = grid1_mat$`Std Control Z`[i], 
                                        sdT1X = grid1_mat$`Std T1 X`[i],
                                        sdT1Z = grid1_mat$`Std T1 Z`[i],
                                        sdT2Y = grid1_mat$`Std T2 Y`[i],
                                        sdT2Z = grid1_mat$`Std T2 Z`[i],
                                        allocProbRandXY_vec = grid1_mat$`Alloc Prob Rand XY`[[i]],
                                        allocProbRandZ_vec = grid1_mat$`Alloc Prob Rand Z`[[i]],
                                        blocksizeXY = grid1_mat$`Blocksize XY`[i],
                                        blocksizeZ = grid1_mat$`Blocksize Z`[i],
                                        complete = grid1_mat$`Complete Randomization`[i],
                                        random = grid1_mat$`Random Allocation`[i],
                                        full = grid1_mat$`Full Output`[i],
                                        iterations = grid1_mat$`Number of Iterations`[i],
                                        alpha = grid1_mat$`Signficance level`[i])
    
  }
  
  # turn the list into a data frame
  data1 <- do.call(rbind,
                    results_list)
  
  # save data frame as csv file with the date of the respective day
  today <- format(Sys.Date(), "%Y_%m_%d")
  write.csv2(data1,
             paste0("data_EqualPrevalence", today, ".csv"),
             row.names = TRUE)

}

#####################################################################################################################
# Example to call the simulation for the following setting:
# for an equal prevalence to the groups, four different sample sizes, different treatment effects for T1 and T2, 5000
# iterations, block and complete randomization and the three different randomization strategies

data <- runSimulation(t = 1,
                      n = c(10, 25, 50, 100),
                      deltaT1 = seq(from = 0, to = 0.5, by = 0.1),
                      deltaT2 = seq(from = 0, to = 0.5, by = 0.1),
                      meanCX =  c(0, -0.5),
                      meanCY =  c(0, 0.5),
                      meanCZ = c(0, 0.5),
                      sdCX = 1,
                      sdCY = 1,
                      sdCZ = 1,
                      sdT1X = 1,
                      sdT1Z = 1,
                      sdT2Y = 1,
                      sdT2Z = 1,
                      iterations = 5000,
                      alpha = 0.025,
                      complete =  c("True", "False"),
                      full = FALSE,
                      random = FALSE,
                      prevalence = list(c((1/3), (1/3), (1/3))),
                      allocProbRandXY_vec = list(c((1/3), (1/3), (1/3)), c((1/2), (1/2))),
                      allocProbRandZ_vec = list(c((1/3), (1/3), (1/3)), c((1/4), (1/4), (1/4), (1/4))),
                      blocksizeXY = c(2, 3),
                      blocksizeZ = c(3, 4)
                      )
