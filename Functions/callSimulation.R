### Libraries
library(parallel) # to detect number of cores
library(foreach) # to parallelize simulation
library(doParallel) # to parallelize simulation

### Use source function to use functions that are created in another R script
source("determineDataset.R")

#### Function for calling the other function
callSimulation <-
  function(
    t, # number of time periods
    n, # number of patients for each time period
    armsNumb_vec = 1:3, # number of arms
    allocProb_vec = c((1/3), (1/3), (1/3)), # allocation probability to the different types of groups
    deltaT1, # treatment effect for treatment 1
    deltaT2, # treatment effect for treatment 2
    meanCX, # mean for control in group X
    meanCY, # mean for control in group Y
    meanCZ, # mean for control in group Z
    sdT1X = 1, # standard deviation for treatment 1 in group X
    sdT2Y = 1, # standard deviation for treatment 2 in group Y
    sdT1Z = 1, # standard deviation for treatment 1 in group Z
    sdT2Z = 1, # standard deviation for treatment 2 in group Z
    sdCX = 1, # standard deviation for control in group X
    sdCY = 1, # standard deviation for control in group Y
    sdCZ = 1, # standard deviation for control in group 
    blocksizeXY, # block size for block randomization for groups X and Y
    blocksizeZ, # block size for block randomization for group Z
    alpha = 0.025, # alpha level for t-test
    allocProbRandXY_vec, # allocation probability for groups X & Y to treatment or control
    allocProbRandZ_vec, # allocation probability for group Z to different arms
    complete, # boolean for randomization method
    random, # boolean for deterministic or random allocation to groups
    full,  # boolean for returning input parameters or just the output
    iterations # number of iterations
  ) {
    
    # determine how many cores to use
    nCores <- parallel::detectCores() - 1
    
    #create the cluster
    myCluster <- parallel::makeCluster(
      nCores, 
      type = "PSOCK")
    
    #register the cluster to be used by %dopar%
    doParallel::registerDoParallel(cl = myCluster)
    
    # create empty list
    opChar_list <- list()
    
    # run the simulation for a given number of iterations  over the grid of all combinations
    opChar_list <- foreach(i = 1:iterations, .packages = c('emmeans', 'dplyr', 'stats'),
                           .export = "fnSimulation") %dopar% {
      
      # call the function 
      opChar_list[[i]] <- fnSimulation(t = t, # number of time periods
                                       n = n, # number of patients for each time period
                                       armsNumb_vec = armsNumb_vec, # number of arms
                                       allocProb_vec = allocProb_vec, # allocation probability to the different types of groups
                                       deltaT1 = deltaT1, # treatment effect for treatment 1
                                       deltaT2 = deltaT2, # treatment effect for treatment 2
                                       meanCX = meanCX, # mean for control in group X
                                       meanCY = meanCY, # mean for control in group Y
                                       meanCZ = meanCZ, # mean for control in group Z
                                       sdT1X = 1, # standard deviation for treatment 1 in group X
                                       sdT2Y = 1, # standard deviation for treatment 2 in group Y
                                       sdT1Z = 1, # standard deviation for treatment 1 in group Z
                                       sdT2Z = 1, # standard deviation for treatment 2 in group Z
                                       sdCX = 1, # standard deviation for control in group X
                                       sdCY = 1, # standard deviation for control in group Y
                                       sdCZ = 1, # standard deviation for control in group 
                                       blocksizeXY = blocksizeXY, # block size for block randomization for groups X and Y
                                       blocksizeZ = blocksizeZ, # block size for block randomization for group Z
                                       alpha = 0.025, # alpha level for t-test
                                       allocProbRandXY_vec = allocProbRandXY_vec, # allocation probability for groups X & Y to treatment or control
                                       allocProbRandZ_vec, # allocation probability for group Z to different arms
                                       complete = complete, # boolean for randomization method
                                       random = random, # boolean for deterministic or random allocation to groups
                                       full = full # boolean for returning input parameters or just the output
      )$list2 
    }
      # stop cluster when no longer needed
      parallel::stopCluster(cl = myCluster)
      
      # calculate the means for treatment 1 and treatment 2 in the different groups to make them part of the output and show them in the csv
      meanT1X <- meanCX + deltaT1
      meanT2Y <- meanCY + deltaT2
      meanT1Z <- meanCZ + deltaT1
      meanT2Z <- meanCZ + deltaT2
      
      # get the mean for the prevalence for the different types of groups over all iterations
      listAllocAll_mat <- do.call(rbind,
                                  lapply(opChar_list,
                                         function(x) table(x[["AllocationList"]])))
      groupMean <- colMeans(listAllocAll_mat)
      groupMean <- unname(groupMean)
      
      # save the mean for each group separately
      groupMeanX <- groupMean[1]
      groupMeanY <- groupMean[2]
      groupMeanZ <- groupMean[3]
      
      # get the mean for the number of patients per arm over all iterations
      tablePatientsPerArm <- do.call(rbind,
                                     lapply(opChar_list,
                                            function(x) table(x[["PatientsperArm"]])))
      armMean <- colMeans(tablePatientsPerArm)
      armMean <- unname(armMean)
      
      # save the mean for each arm separately
      meanControl <- armMean[1]
      meanT1 <- armMean[2]
      meanT2 <- armMean[3]
      
      ### calculate power, bias, type I error and CI for the different scenarios ###
      
      ## Analysis Scenario 1
      ## one-way ANOVA control vs treatment 1 using all control data
      pValCT1All <- do.call(rbind,
                            lapply(opChar_list, 
                                   function(x) x[["PValueControlT1All"]]))
      
      # save which p-values are smaller than alpha 
      rejectOneCT1All_mat <- pValCT1All < alpha
      
      # calculate the rejection probability
      rejProbCT1All <- mean(rejectOneCT1All_mat)
      
      # get estimates and calculate bias
      estimateCT1All <- do.call(rbind,
                                lapply(opChar_list,
                                       function(x) x[["estimateCT1All"]]))
      biasCT1All <- mean(estimateCT1All - deltaT1)
      
      # get lower and upper bound of confidence intervals and calculate the average
      CICT1All <- do.call(rbind,
                          lapply(opChar_list,
                                 function(x) x[["CICT1All"]]))
      CICT1All <- colMeans(CICT1All)
      
      
      
      ## Analysis Scenario 1  
      ## one-way ANOVA control vs treatment 2 using all control data
      pValCT2All <- do.call(rbind,
                            lapply(opChar_list,
                                   function(x) x[["PValueControlT2All"]]))
      
      # save which p-values are smaller than alpha 
      rejectOneCT2All_mat <- pValCT2All < alpha
      
      # calculate the rejection probability
      rejProbCT2All <- mean(rejectOneCT2All_mat)
      
      # check which p-values are smaller than alpha for both T1 and T2 and T1 or T2 
      rejectAtLeastOneCT12All_mat <- (pValCT1All < alpha) | (pValCT2All < alpha)
      rejectBothCT12All_mat <- (pValCT1All < alpha) & (pValCT2All < alpha)
      
      # get estimates and calculate bias
      estimateCT2All <- do.call(rbind,
                                lapply(opChar_list,
                                       function(x) x[["estimateCT2All"]]))
      biasCT2All <- mean(estimateCT2All - deltaT2)
      
      # get lower and upper bound of confidence intervals and calculate the average
      CICT2All <- do.call(rbind,
                          lapply(opChar_list,
                                 function(x) x[["CICT2All"]]))
      CICT2All <- colMeans(CICT2All)   
      
      
      
      ## Analysis Scenario 3  
      ## one-way ANOVA control vs treatment 1 using only suitable control data of groups X and Z
      pValCT1Fit <- do.call(rbind,
                            lapply(opChar_list, 
                                   function(x) x[["pValCT1Fit"]]))
      
      # save which p-values are smaller than alpha 
      rejectOneCT1Fit_mat <- pValCT1Fit < alpha
      
      # calculate the rejection probability
      rejProbCT1Fit <- mean(rejectOneCT1Fit_mat)
      
      # get estimates and calculate bias
      estimateCT1Fit <- do.call(rbind,
                                lapply(opChar_list,
                                       function(x) x[["estimateCT1Fit"]]))
      biasCT1Fit <- mean(estimateCT1Fit - deltaT1)
      
      # get lower and upper bound of confidence intervals and calculate the average
      CICT1Fit <- do.call(rbind,
                          lapply(opChar_list,
                                 function(x) x[["CICT1Fit"]]))
      CICT1Fit1 <- colMeans(CICT1Fit)  
      
      
      
      ## Analysis Scenario 3  
      ## one-way ANOVA control vs treatment 2 using only suitable control data of groups Y and Z
      pValCT2Fit <- do.call(rbind,
                            lapply(opChar_list, 
                                   function(x) x[["pValCT2Fit"]]))
      
      # save which p-values are smaller than alpha 
      rejectOneCT2Fit_mat <- pValCT2Fit < alpha
      
      # calculate the rejection probability
      rejProbCT2Fit <- mean(rejectOneCT2Fit_mat)
      
      # check which p-values are smaller than alpha for both T1 and T2 and T1 or T2  
      rejectAtLeastOneCT12Fit_mat <- (pValCT1Fit < alpha) |(pValCT2Fit < alpha)
      rejectBothCT12Fit_mat <- (pValCT1Fit < alpha) & (pValCT2Fit < alpha)
      
      # get estimates and calculate bias
      estimateCT2Fit <- do.call(rbind,
                                lapply(opChar_list,
                                       function(x) x[["estimateCT2Fit"]]))
      biasCT2Fit <- mean(estimateCT2Fit - deltaT2)
      
      # get lower and upper bound of confidence intervals and calculate the average
      CICT2Fit <- do.call(rbind,
                          lapply(opChar_list,
                                 function(x) x[["CICT2Fit"]]))
      CICT2Fit <- colMeans(CICT2Fit)  
      
      
      
      
      ## Analysis Scenario 2  
      ## one-way ANOVA control vs treatment 1 using all control data of the data set which includes T1 and C
      pValCT1FitAllC <- do.call(rbind,
                                lapply(opChar_list, 
                                       function(x) x[["pValCT1FitAllC"]]))
      
      # save which p-values are smaller than alpha 
      rejectOneCT1FitAllC_mat <- pValCT1FitAllC < alpha
      
      # calculate the rejection probability
      rejProbCT1FitAllC <- mean(rejectOneCT1FitAllC_mat)
      
      # get estimates and calculate bias
      estimateCT1FitAllC <- do.call(rbind,
                                    lapply(opChar_list,
                                           function(x) x[["estimateCT1FitAllC"]]))
      biasCT1FitAllC <- mean(estimateCT1FitAllC - deltaT1)
      
      # get lower and upper bound of confidence intervals and calculate the average
      CICT1FitAllC <- do.call(rbind,
                              lapply(opChar_list,
                                     function(x) x[["CICT1FitAllC"]]))
      CICT1FitAllC1 <- colMeans(CICT1FitAllC) 
      
      
      
      ## Analysis Scenario 2  
      ## one-way ANOVA control vs treatment 2 using all control data of the data set which includes T2 and C
      pValCT2FitAllC <- do.call(rbind,
                                lapply(opChar_list, 
                                       function(x) x[["pValCT2FitAllC"]]))
      
      # save which p-values are smaller than alpha 
      rejectOneCT2FitAllC_mat <- pValCT2FitAllC < alpha
      
      # calculate the rejection probability
      rejProbCT2FitAllC <- mean(rejectOneCT2FitAllC_mat)
      
      # check which p-values are smaller than alpha for both T1 and T2 and T1 or T2 
      rejectAtLeastOneCT12FitAllC_mat <- (pValCT1FitAllC < alpha) |(pValCT2FitAllC < alpha)
      rejectBothCT12FitAllC_mat <- (pValCT1FitAllC < alpha) & (pValCT2FitAllC < alpha)
      
      # get estimates and calculate bias
      estimateCT2FitAllC <- do.call(rbind,
                                    lapply(opChar_list,
                                           function(x) x[["estimateCT2FitAllC"]]))
      biasCT2FitAllC <- mean(estimateCT2FitAllC - deltaT2)
      
      # get lower and upper bound of confidence intervals and calculate the average
      CICT2FitAllC <- do.call(rbind,
                              lapply(opChar_list,
                                     function(x) x[["CICT2FitAllC"]]))
      CICT2FitAllC <- colMeans(CICT2FitAllC)  
      
      
      ## Analysis Scenario 4  
      ## two-way ANOVA control vs treatment 1 for all data
      pValCT1Two <- do.call(rbind,
                            lapply(opChar_list,
                                   function(x) x[["pValCT1Two"]]))
      
      # save which p-values are smaller than alpha 
      rejectOne2CT1_mat <- pValCT1Two < alpha
      
      # calculate the rejection probability
      rejProbCT1Two <- mean(rejectOne2CT1_mat)
      
      # get estimates and calculate bias
      estimateCT1Two <- do.call(rbind,
                                lapply(opChar_list,
                                       function(x) x[["estimateCT1Two"]]))
      biasCT1Two <- mean(estimateCT1Two - deltaT1)
      
      # get lower and upper bound of confidence intervals and calculate the average
      CICT1Two <- do.call(rbind,
                          lapply(opChar_list,
                                 function(x) x[["CICT1Two"]]))
      CICT1Two <- colMeans(CICT1Two)
      
      
      ## Analysis Scenario 4  
      ## two-way ANOVA control vs treatment 2 for all data
      pValCT2Two <- do.call(rbind,
                            lapply(opChar_list, 
                                   function(x) x[["pValCT2Two"]]))
      
      # save which p-values are smaller than alpha 
      rejectOne2CT2_mat <- pValCT2Two < alpha
      
      # calculate the rejection probability
      rejProbCT2Two <- mean(rejectOne2CT2_mat)
      
      # check which p-values are smaller than alpha for both T1 and T2 and T1 or T2 and take their mean 
      rejectAtLeastOne2CT12_mat <- (pValCT1Two < alpha) | (pValCT2Two < alpha)
      rejectBoth2CT12_mat <- (pValCT1Two < alpha) & (pValCT2Two < alpha)
      
      # get estimates and calculate bias
      estimateCT2Two <- do.call(rbind, 
                                lapply(opChar_list,
                                       function(x) x[["estimateCT2Two"]]))
      biasCT2Two <- mean(estimateCT2Two - deltaT2)
      
      # get lower and upper bound of confidence intervals and calculate the average
      CICT2Two <- do.call(rbind,
                          lapply(opChar_list,
                                 function(x) x[["CICT2Two"]]))
      CICT2Two <- colMeans(CICT2Two)
      
      
      
      ## Analysis Scenario 6  
      # two-way ANOVA control vs treatment 1 using the data set which only includes T1 and C
      pValCT1XZ <- do.call(rbind,
                           lapply(opChar_list,
                                  function(x) x[["pValCT1XZ"]]))
      
      # save which p-vales are smaller than alpha 
      rejectOne2CT1XZ_mat <- pValCT1XZ < alpha
      
      # calculate the rejection probability
      rejProbCT1XZ <- mean(rejectOne2CT1XZ_mat)
      
      # get estimates and calculate bias
      estimateCT1XZ <- do.call(rbind,
                               lapply(opChar_list,
                                      function(x) x[["estimateCT1XZ"]]))
      biasCT1XZ <- mean(estimateCT1XZ - deltaT1)
      
      # get lower and upper bound of confidence intervals and calculate the average
      CICT1XZ <- do.call(rbind,
                         lapply(opChar_list,
                                function(x) x[["CICT1XZ"]]))
      
      CICT1XZ <- colMeans(CICT1XZ)  
      
      
      
      ## Analysis Scenario 6
      # two-way ANOVA control vs. treatment 2  using the data set which only includes T2 and C
      pValCT2YZ <- do.call(rbind,
                           lapply(opChar_list,
                                  function(x) x[["pValCT2YZ"]]))
      
      # save which p-vales are smaller than alpha 
      rejectOne2CT2YZ_mat <- pValCT2YZ < alpha
      
      # calculate the rejection probability
      rejProbCT2YZ <- mean(rejectOne2CT2YZ_mat)
      
      # check which p-values are smaller than alpha for both T1 and T2 and T1 or T2 and take their mean
      rejectAtLeastOne2CT12Fit_mat <- (pValCT1XZ < alpha) | (pValCT2YZ < alpha)
      rejectBoth2CT12Fit_mat <- (pValCT1XZ < alpha) & (pValCT2YZ < alpha)
      
      # get estimates and calculate bias
      estimateCT2YZ <- do.call(rbind,
                               lapply(opChar_list,
                                      function(x) x[["estimateCT2YZ"]]))
      biasCT2YZ <- mean(estimateCT2YZ - deltaT2)
      
      # get lower and upper bound of confidence intervals and calculate the average
      CICT2YZ <- do.call(rbind,
                         lapply(opChar_list,
                                function(x) x[["CICT2YZ"]]))
      CICT2YZ <- colMeans(CICT2YZ)
      
      
      
      ## Analysis Scenario 5
      # two-way ANOVA control vs. treatment 1  using the data set which only includes all control data but not T2
      pValCT1AllC <- do.call(rbind,
                             lapply(opChar_list,
                                    function(x) x[["pValCT1AllC"]]))
      
      # save which p-vales are smaller than alpha 
      rejectOne2CT1AllC_mat <- pValCT1AllC < alpha
      
      # calculate the rejection probability
      rejProbCT1AllC <- mean(rejectOne2CT1AllC_mat)
      
      # get estimates and calculate bias
      estimateCT1AllC <- do.call(rbind,
                                 lapply(opChar_list,
                                        function(x) x[["estimateCT1AllC"]]))
      biasCT1AllC <- mean(estimateCT1AllC - deltaT1)
      
      # get lower and upper bound of confidence intervals and calculate the average
      CICT1AllC <- do.call(rbind,
                           lapply(opChar_list,
                                  function(x) x[["CICT1AllC"]]))
      CICT1AllC <- colMeans(CICT1AllC)  
      
      
      
      ## Analysis Scenario 5
      # two-way ANOVA control vs. treatment 2 using the data set which only includes all control data but not T1
      pValCT2AllC <- do.call(rbind,
                             lapply(opChar_list,
                                    function(x) x[["pValCT2AllC"]]))
      
      # save which p-vales are smaller than alpha 
      rejectOne2CT2AllC_mat <- pValCT2AllC < alpha
      
      # calculate the rejection probability
      rejProbCT2AllC <- mean(rejectOne2CT2AllC_mat)
      
      # check which p-values are smaller than alpha for both T1 and T2 and T1 or T2 and take their mean
      rejectAtLeastOne2CT12AllC_mat <- (pValCT1AllC < alpha) | (pValCT2AllC < alpha)
      rejectBoth2CT12AllC_mat <- (pValCT1AllC < alpha) & (pValCT2AllC < alpha)
      
      # get estimates and calculate bias
      estimateCT2AllC <- do.call(rbind,
                                 lapply(opChar_list,
                                        function(x) x[["estimateCT2AllC"]]))
      biasCT2AllC <- mean(estimateCT2AllC - deltaT2)
      
      # get lower and upper bound of confidence intervals and calculate the average
      CICT2AllC <- do.call(rbind,
                           lapply(opChar_list,
                                  function(x) x[["CICT2AllC"]]))
      CICT2AllC <- colMeans(CICT2AllC)  
      
      
      # calculate the conjunctive and disjunctive power for all analysis scenarios
      # the conjunctive and disjunctive power are calculated when the treatment effects are different from 0
      if(deltaT1 != 0 & deltaT2 != 0) {
        conjPower1All <- mean(rejectBothCT12All_mat)
        disjPower1All <- mean(rejectAtLeastOneCT12All_mat)
        conjPower1Fit <- mean(rejectBothCT12Fit_mat)
        disjPower1Fit <- mean(rejectAtLeastOneCT12Fit_mat)
        conjPower1AllC <- mean(rejectBothCT12FitAllC_mat)
        disjPower1AllC <- mean(rejectAtLeastOneCT12FitAllC_mat)
        
        conjPower2All <- mean(rejectBoth2CT12_mat)
        disjPower2All <- mean(rejectAtLeastOne2CT12_mat)
        conjPower2Fit <- mean(rejectBoth2CT12Fit_mat)
        disjPower2Fit <- mean(rejectAtLeastOne2CT12Fit_mat)
        conjPower2AllC <- mean(rejectBoth2CT12AllC_mat)
        disjPower2AllC <- mean(rejectAtLeastOne2CT12AllC_mat)
        
      }  else {
        conjPower1All <- NA
        disjPower1All <- NA
        conjPower1Fit <- NA
        disjPower1Fit <- NA
        conjPower1AllC <- NA
        disjPower1AllC <- NA
        
        conjPower2All <- NA
        disjPower2All <- NA
        conjPower2Fit <- NA
        disjPower2Fit <- NA
        conjPower2AllC <- NA
        disjPower2AllC <- NA
      }
      
      # return a list of the operating characteristics of interest
      results_list <- list(#opChar_list = opChar_list,
                          Timepoints = t,
                          SampleSize = n*t,
                          NumbArms = length(armsNumb_vec),
                          Prevalence = allocProb_vec, 
                          MeanCX = meanCX,
                          MeanCY = meanCY,
                          MeanCZ = meanCZ,
                          MeanT1X = meanT1X, 
                          MeanT2Y = meanT2Y, 
                          MeanT1Z = meanT1Z,
                          MeanT2Z = meanT2Z,
                          Iterations = iterations,
                          TreatmentEffectT1 = deltaT1,
                          TreatmentEffectT2 = deltaT2,
                          blocksizeXY = blocksizeXY, 
                          blocksizeZ = blocksizeZ,
                          StandDevT1X = sdT1X, 
                          StandDevT2Y = sdT2Y,
                          StandDevT1Z = sdT1Z,
                          StandDevT2Z = sdT2Z,
                          StandDevCX = sdCX,
                          StandDevCY = sdCY,
                          StandDevCZ = sdCZ,
                          AllocProbXY = allocProbRandXY_vec, 
                          AllocProbZ = allocProbRandZ_vec, 
                          CompleteRand = complete,
                          RandomAlloc = random, 
                          MeanGroupX = groupMeanX,
                          MeanGroupY = groupMeanY,
                          MeanGroupZ = groupMeanZ,
                          MeanControl = meanControl,
                          MeanT1 = meanT1,
                          MeanT2 = meanT2,
                          rejProbCT1All = rejProbCT1All,
                          biasCT1All = biasCT1All,
                          CICT1All = CICT1All,
                          rejProbCT2All = rejProbCT2All,
                          biasCT2All = biasCT2All,
                          CICT2All = CICT2All,
                          rejProbCT1Fit = rejProbCT1Fit,
                          biasCT1Fit = biasCT1Fit,
                          CICT1Fit = CICT1Fit1,
                          rejProbCT2Fit = rejProbCT2Fit,
                          biasCT2Fit = biasCT2Fit,
                          CICT2Fit = CICT2Fit,
                          rejProbCT1FitAllC = rejProbCT1FitAllC,
                          biasCT1FitAllC = biasCT1FitAllC,
                          CICT1FitAllC = CICT1FitAllC1,
                          rejProbCT2FitAllC = rejProbCT2FitAllC,
                          biasCT2FitAllC = biasCT2FitAllC,
                          CICT2FitAllC = CICT2FitAllC,
                          rejProbCT1Two = rejProbCT1Two,
                          biasCT1Two = biasCT1Two,
                          CICT1Two = CICT1Two,
                          rejProbCT2Two = rejProbCT2Two,
                          biasCT2Two = biasCT2Two, 
                          CICT2Two = CICT2Two,
                          rejProbCT1XZ = rejProbCT1XZ,
                          biasCT1XZ = biasCT1XZ,
                          CICT1XZ = CICT1XZ,
                          rejProbCT2YZ = rejProbCT2YZ,
                          biasCT2YZ = biasCT2YZ,
                          CICT2YZ = CICT2YZ,
                          rejProbCT1AllC = rejProbCT1AllC,
                          biasCT1AllC = biasCT1AllC,
                          CICT1AllC = CICT1AllC,
                          rejProbCT2AllC = rejProbCT2AllC,
                          biasCT2AllC = biasCT2AllC,
                          CICT2AllC = CICT2AllC,
                          ConjPower1AllData = conjPower1All,
                          DisjPower1AllData = disjPower1All,
                          ConjPower1FittedData = conjPower1Fit,
                          DisjPower1FittedData = disjPower1Fit,
                          ConjPower1AllControlData = conjPower1AllC,
                          DisjPower1AllControlData = disjPower1AllC,
                          ConjPower2AllData = conjPower2All,
                          DisjPower2AllData = disjPower2All,
                          ConjPower2FittedData = conjPower2Fit,
                          DisjPower2FittedData = disjPower2Fit,
                          ConjPower2AllControlData = conjPower2AllC,
                          DisjPower2AllControlData = disjPower2AllC
      )
    
    return(results_list)
}

#######################################################################################################
# Example for calling the function for the following setting:
# An overall ratio of 1:1:2 with block randomization and random allocation to the different types of 
# groups, so within the substudies a 1:1 ratio for treatment(s) versus control.
                                  
data <- callSimulation(t = 5,
                       n = 60,
                       deltaT1 = 0,
                       deltaT2 = 0.5,
                       meanCX = 0,
                       meanCY = 0,
                       meanCZ = 0,
                       alpha = 0.025,
                       allocProb_vec = c((1/3), (1/3), (1/3)),
                       allocProbRandXY_vec = c((1/2), (1/2)),
                       allocProbRandZ_vec = c((1/4), (1/4), (1/4), (1/4)),
                       blocksizeXY = 2,
                       blocksizeZ = 4,
                       complete = FALSE,
                       random = TRUE,
                       iterations = 10,
                       full = FALSE)
