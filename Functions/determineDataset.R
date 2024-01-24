#### Libraries
library(stats) # library for normal distribution
library(emmeans) # library for contrasts and lsmeans function
library(dplyr) # library for filter function

#### Function to generate dataset
fnSimulation <- 
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
    full # boolean for returning full output including input parameters or just the output
  ) {
    
    # set a seed for reproducibility
    set.seed(seed = 234)
    
    # calculate the overall sample size
    sampleSize <- n*t
    
    # create ID vector
    ID_vec <- c(1:sampleSize)
    
    # set counting variable for the different groups to 0
    nGroupX <- nGroupY <- nGroupZ <- nGroup <- 0
    
    # create empty data frame and vectors for future values
    time <- arm_char <- metOutcome1_vec <- tmp_vec <- vector()
    
    data <- data.frame(time,
                       arm_char,
                       tmp_vec,
                       metOutcome1_vec)
    
    # calculate the means for treatment 1 and treatment 2 in the different groups 
    meanT1X <- meanCX + deltaT1
    meanT2Y <- meanCY + deltaT2
    meanT1Z <- meanCZ + deltaT1
    meanT2Z <- meanCZ + deltaT2
    
    
    if(random == TRUE) {
      
      # determine the list of allocation for all patients when getting allocated randomly
      listAlloc_vec <- sample(x = armsNumb_vec,
                              size = sampleSize,
                              replace = TRUE,
                              prob = allocProb_vec)
    } else {
      
      # calculate how many patients should get allocated to each group for a 
      # deterministic allocation 
      # multiply by sample size because times in the next step only takes integers
      allocProb_vec1 <- allocProb_vec * sampleSize
      
      # determine the list of allocation for all patients
      listAlloc_vec <- sample(x = c(rep(armsNumb_vec,
                                        times = ceiling(allocProb_vec1))),
                              size = sampleSize,
                              replace = FALSE)
    }
    
    if(complete == TRUE){
      
      # the randomization list depends on the block size which is different for 
      # the different ratios: if the block size is smaller or equal to the 
      # number of treatments, there is no replacement within each block
      
      if(length(allocProbRandXY_vec) == 2) {
        # determine the list of complete randomization for patients in group X
        # list has length of n*t because the most extreme case is that 
        # all patients are in one group
        armX_char <- sample(x = c("T1", "C"),
                            size = sampleSize,
                            replace = TRUE,
                            prob = allocProbRandXY_vec)
        
        # determine the list of randomization for patients in group Y
        # list has length of n*t because the most extreme case is that 
        # all patients are in one group
        armY_char <- sample(x = c("T2", "C"),
                            size = sampleSize,
                            replace = TRUE,
                            prob = allocProbRandXY_vec)
      } else {
        
        # determine the list of complete randomization for patients in group X
        # list has length of n*t because the most extreme case is that 
        # all patients are in one group
        armX_char <- sample(x = c("T1", "T1", "C"),
                            size = sampleSize,
                            replace = TRUE,
                            prob = allocProbRandXY_vec)
        
        # determine the list of randomization for patients in group Y
        # list has length of n*t because the most extreme case is that 
        # all patients are in one group
        armY_char <- sample(x = c("T2", "T2", "C"),
                            size = sampleSize,
                            replace = TRUE,
                            prob = allocProbRandXY_vec)
      }  
      
      if(length(allocProbRandZ_vec) == 3) {
        # determine the list of randomization for patients in group Z
        # list has length of n*t because the most extreme case is that 
        # all patients are in one group
        armZ_char <- sample(x = c("T1", "T2", "C"),
                            size = sampleSize,
                            replace = TRUE,
                            prob = allocProbRandZ_vec)
      } else {
        # determine the list of randomization for patients in group Z
        # list has length of n*t because the most extreme case is that 
        # all patients are in one group
        armZ_char <- sample(x = c("T1", "T2", "C", "C"),
                            size = sampleSize,
                            replace = TRUE,
                            prob = allocProbRandZ_vec)
      }
      
      # use block randomization
    } else {
      
      # determine the length of the randomization list 
      nBlockXY_vec <- rep(1:ceiling(sampleSize/blocksizeXY),
                          each = blocksizeXY)
      
      # if the block size is smaller or equal to the number of treatments,
      # there is no replacement within each block
      if(blocksizeXY <= 2) {
        
        # calculate the first block for group X: each treatment occurs once and 
        # the probability to get allocated to one of the treatments is the same
        randListX_vec <- sample(x = c("T1", "C"),
                                size = blocksizeXY,
                                replace = FALSE,
                                prob = allocProbRandXY_vec)
        
        # determine the randomization list by making sure that 
        # the treatments are allocated equally
        while(length(randListX_vec) <= length(nBlockXY_vec)) {
          
          randListX_vec <- c(randListX_vec, sample(x = c("T1", "C"),
                                                   size = blocksizeXY,
                                                   replace = FALSE,
                                                   prob = allocProbRandXY_vec))
        }
        
        # calculate the first block for group Y: each treatment occurs once and 
        # the probability to get allocated to one of the treatments is the same
        randListY_vec <- sample(x = c("T2", "C"),
                                size = blocksizeXY,
                                replace = FALSE,
                                prob = allocProbRandXY_vec)
        
        # determine the randomization list by making sure that 
        # the treatments are allocated equally
        while(length(randListY_vec) <= length(nBlockXY_vec)) {
          
          randListY_vec <- c(randListY_vec, sample(x = c("T2", "C"),
                                                   size = blocksizeXY,
                                                   replace = FALSE,
                                                   prob = allocProbRandXY_vec))
        }
      } else {
        # calculate the first block for group X: each treatment can occur more than once 
        # and the probability to get allocated to one of the treatments is not the same
        randListX_vec <- sample(x = c(rep("T1", blocksizeXY/blocksizeXY * 2),
                                      "C"),
                                size = blocksizeXY,
                                replace = FALSE,
                                prob = allocProbRandXY_vec)
        
        # determine the randomization list 
        while(length(randListX_vec) <= length(nBlockXY_vec)) {
          
          randListX_vec <- c(randListX_vec, sample(x = c(rep("T1", blocksizeXY/blocksizeXY * 2),
                                                         "C"),
                                                   size = blocksizeXY,
                                                   replace = FALSE,
                                                   prob = allocProbRandXY_vec))
        } 
        
        # calculate the first block for group Y: the treatment can occur more than once and the 
        # probability to get allocated to one of the treatments is the same
        randListY_vec <- sample(x = c(rep("T2", blocksizeXY/blocksizeXY * 2), 
                                      "C"),
                                size = blocksizeXY,
                                replace = FALSE,
                                prob = allocProbRandXY_vec)
        
        # determine the randomization list 
        while(length(randListY_vec) <= length(nBlockXY_vec)) {
          
          randListY_vec <- c(randListY_vec, sample(x = c(rep("T2", blocksizeXY/blocksizeXY * 2),
                                                         "C"),
                                                   size = blocksizeXY,
                                                   replace = FALSE,
                                                   prob = allocProbRandXY_vec))
        }
      }
      
      # determine the length of the randomization list 
      nBlockZ_vec <- rep(1:ceiling(sampleSize/blocksizeZ),
                         each = blocksizeZ)
      
      # if the block size is smaller or equal to the number of treatments,
      # there is no replacement within each block
      if(blocksizeZ <= 3) {
        # calculate the first block for group Z: each treatment occurs once and 
        # the probability to get allocated to one of the treatments is the same
        randListZ_vec <- sample(x = c("T1", "T2", "C"),
                                size = blocksizeZ,
                                replace = FALSE,
                                prob = allocProbRandZ_vec)
        
        # determine the randomization list by making sure that 
        # the treatments are allocated equally
        while(length(randListZ_vec) <= length(nBlockZ_vec)) {
          
          randListZ_vec <- c(randListZ_vec, sample(x = c("T1", "T2", "C"),
                                                   size = blocksizeZ,
                                                   replace = FALSE,
                                                   prob = allocProbRandZ_vec))
        }
      } else{
        # calculate the first block for group Z: each treatment can occur more than once 
        # and the probability to get allocated to one of the treatments is not the same
        randListZ_vec <- sample(x = c("T1", "T2", "C", "C"),
                                size = blocksizeZ,
                                replace = FALSE,
                                prob = allocProbRandZ_vec)
        
        # determine the randomization list 
        while(length(randListZ_vec) <= length(nBlockZ_vec)) {
          
          randListZ_vec <- c(randListZ_vec, sample(x = c("T1", "T2", "C", "C"),
                                                   size = blocksizeZ,
                                                   replace = FALSE,
                                                   prob = allocProbRandZ_vec))
        }
      }
    }
    
    # repeat the simulation for a given number of time periods
    for(time in 1:t) {
      # repeat the simulation for a given number of patients
      for (i in 1:n) {

          # when considering different groups of patients
          # i.e. patients are not willing to get randomized to all three arms,
          # the prevalence can either be determined at random or deterministic
          
          # determine the prevalence at random
          if(random == TRUE) {
            
            # 1: group X = treatment 1, but not treatment 2
            # 2: group Y = treatment 2, but not treatment 1
            # 3: group Z = treatment 1 and 2
            
            # counting variable
            nGroup <- nGroup + 1
            
            # iterate over all patients and assign the corresponding 
            # type of group from the allocation list
            typeofgroup <- listAlloc_vec[nGroup]
            
            # determine the prevalence deterministic  
          } else {
            
            # 1: group X = treatment 1, but not treatment 2
            # 2: group Y = treatment 2, but not treatment 1
            # 3: group Z = treatment 1 and 2
            
            # counting variable
            nGroup <- nGroup + 1
            
            # iterate over all patients and assign the corresponding 
            # type of group from the allocation list
            typeofgroup <- listAlloc_vec[nGroup]
          }
          
          # determine the allocation to treatment or control arm depending on the
          # type of group by complete or block randomization
          
          if(complete == TRUE) {
            
            # patients in group 1 are either randomized to treatment 1 or control
            if(typeofgroup == 1) {
              
              # counting variable within group 1
              nGroupX <- nGroupX + 1
              
              # iterate over the number of patients in the group and assign
              # the corresponding treatment from the randomization list
              arm_char <- armX_char[nGroupX]
              
              # patients in group 2 are either randomized to treatment 2 or control
            } else if(typeofgroup == 2) { 
              
              # counting variable within group 2
              nGroupY <- nGroupY + 1
              
              # iterate over the number of patients in the group and assign
              # the corresponding treatment from the randomization list
              arm_char <- armY_char[nGroupY]
              
              # patients in group 3 are either randomized to treatment 1, 2 or control
            } else {  
              
              # counting variable in group 3
              nGroupZ <- nGroupZ + 1
              
              # iterate over the number of patients in the group and assign
              # the corresponding treatment from the randomization list
              arm_char <- armZ_char[nGroupZ]
            }
            
            # use block randomization
          } else {
            
            # for group X: patient is randomized to either treatment 1 or control 
            if(typeofgroup == 1) {
              
              # count the patients in group 1
              nGroupX <- nGroupX + 1
              
              # iterate over the number of patients in the group and assign 
              # the corresponding treatment from the randomization list
              arm_char <- randListX_vec[nGroupX]
              
              # for group Y: patient is randomized to either treatment 2 or control
            } else if(typeofgroup == 2) { 
              
              # count the patients in group 2
              nGroupY <- nGroupY + 1
              
              # iterate over the number of patients in the group and assign
              # the corresponding treatment from the randomization list
              arm_char <- randListY_vec[nGroupY]
              
              # for group Z: patient is randomized to either treatment 1, treatment 2 or control
            } else { 
              
              # count the patients in group 3
              nGroupZ <- nGroupZ + 1
              
              # iterate over the number of patients in the group and assign
              # the corresponding treatment from the randomization list
              arm_char <- randListZ_vec[nGroupZ]
            }
          } 
        
        # generate normally distributed outcome variables 
        # depending on the type of group and treatment arm
        if(arm_char == "T1" && typeofgroup == 1) {
          tmp_vec <- rnorm(1,
                           mean = meanT1X, 
                           sd = sdT1X)
        } else if (arm_char == "T2" && typeofgroup == 2) {
          tmp_vec <- rnorm(1,
                           mean = meanT2Y,
                           sd = sdT2Y)
        } else if(arm_char == "T1" && typeofgroup == 3) {
          tmp_vec <- rnorm(1,
                           mean = meanT1Z,
                           sd = sdT1Z)
        } else if(arm_char == "T2" && typeofgroup == 3) {
          tmp_vec <- rnorm(1,
                           mean = meanT2Z,
                           sd = sdT2Z)
        } else if(arm_char == "C" && typeofgroup == 1) {
          tmp_vec <- rnorm(1,
                           mean = meanCX,
                           sd = sdCX)
        } else if (arm_char == "C" && typeofgroup == 2) {
          tmp_vec <- rnorm(1,
                           mean = meanCY,
                           sd = sdCY)
        } else {
          tmp_vec <- rnorm(1,
                           mean = meanCZ,
                           sd = sdCZ)
        }
        
        # return a data frame with four variables: time when included, the type 
        # of group, the treatment arm and the outcome
        data <- rbind(data,
                      cbind(time,
                            typeofgroup, 
                            arm_char,
                            tmp_vec))
        
        # save how many patients are in each arm
        patientsPerArm <- data$arm_char
      }
    }
    # add the ID variable to the data frame and order the rows such that 
    # the ID variable is in the first row
    data$ID <- ID_vec
    data <- data[, c(5, 1:4)]
    
    # rename the rows in the data frame
    colnames(data) <- c("ID",
                        "Time",
                        "Typeofgroups", 
                        "TreatmentArm",
                        "Outcome")
    
    ##### Analysis
    # depending on whether one considers different type of groups and all 
    # control data or only the suitable control data the analysis is different
    
    # Analysis scenario 1
    # fit one-way model for all control data unadjusted for type of 
    # control
    fit1 <- lm(Outcome ~ TreatmentArm,
               data = data)
    
    # use lsmeans to compare the groups 
    fit1means <- lsmeans(fit1, 
                         "TreatmentArm")

    # calculate treatment vs control comparison
    contrastsTC <- contrast(object = fit1means, 
                            method = "trt.vs.ctrl",
                            adjust = "none")
    # calculate the contrasts one-sided
    contrastsTC <- test(contrastsTC, side = ">")
    # transform emmsgrid into data frame to extract values
    contrastsTC <- data.frame(contrastsTC)
    
    # extract the estimates
    estimateCT1All <- contrastsTC[[2]][1]
    estimateCT2All <- contrastsTC[[2]][2]
    # extract the p-values
    pValCT1All <- contrastsTC[[6]][1]
    pValCT2All <- contrastsTC[[6]][2]  
      
    # calculate confidence interval for treatment versus control contrast
    contrastsTCCI <- confint(contrast(object = fit1means, 
                                      method = "trt.vs.ctrl",
                                      adjust = "none"))    
    # extract upper and lower bounds
    CICT1All <- c(contrastsTCCI[[5]][1],
                  contrastsTCCI[[6]][1])
    CICT2All <- c(contrastsTCCI[[5]][2],
                  contrastsTCCI[[6]][2])
    
    
    
    # Analysis Scenario 3
    # edit the data frame and delete all patients in group 2 (Y) for the 
    # comparison of T1-C with only the patients who had the chance of getting
    # treated with the treatment of interest
    dataXZ <- data %>%
      filter(!Typeofgroups %in% '2' & !TreatmentArm %in% "T2")
    
    # fit one-way model 
    fit2 <- lm(Outcome ~ TreatmentArm,
               data = dataXZ)
    
    # use lsmeans to compare the groups
    fit2means <- lsmeans(fit2, 
                         "TreatmentArm")

    # calculate treatment vs control comparison
    contrasts2TC <- contrast(object = fit2means, 
                             method = "trt.vs.ctrl",
                             adjust = "none")
    # calculate the contrasts one-sided
    contrasts2TC <- test(contrasts2TC, side = ">")
    # transform emmsgrid into data frame to extract values
    contrasts2TC <- data.frame(contrasts2TC)
    
    # extract the estimate
    estimateCT1Fit <- contrasts2TC[[2]][1]
    # extract the p-value
    pValCT1Fit <- contrasts2TC[[6]][1]
    
    # calculate confidence interval for treatment versus control contrast
    contrasts2TCCI <- confint(contrast(object = fit2means, 
                                       method = "trt.vs.ctrl",
                                       adjust = "none"))    
    # extract upper and lower bounds
    CICT1Fit <- c(contrasts2TCCI[[5]][1],
                  contrasts2TCCI[[6]][1]) 
      
    
    # Analysis Scenario 3
    # edit the data frame and delete all patients in group 1 (X) for the 
    # comparison of T2-C with only with the patients who had the chance
    # of getting treated with the treatment of interest
    dataYZ <- data %>%
      filter(!Typeofgroups %in% '1' & !TreatmentArm %in% "T1")
    
    # fit one-way anova model
    fit3 <- lm(Outcome ~ TreatmentArm,
               data = dataYZ)

    # use lsmeans to compare the groups
    fit3means <- lsmeans(fit3, 
                         "TreatmentArm")
    
    # calculate treatment vs control comparison
    contrasts3TC <- contrast(object = fit3means, 
                             method = "trt.vs.ctrl",
                             adjust = "none")
    # calculate the contrasts one-sided
    contrasts3TC <- test(contrasts3TC, side = ">")
    # transform emmsgrid into data frame to extract values
    contrasts3TC <- data.frame(contrasts3TC)
    
    # extract the estimate
    estimateCT2Fit <- contrasts3TC[[2]][1]
    # extract the p-value
    pValCT2Fit <- contrasts3TC[[6]][1]  
      
    # calculate confidence interval for treatment versus control contrast
    contrasts3TCCI <- confint(contrast(object = fit3means, 
                                       method = "trt.vs.ctrl",
                                       adjust = "none"))    
    # extract upper and lower bounds of CI
    CICT2Fit<- c(contrasts3TCCI[[5]][1],
                 contrasts3TCCI[[6]][1])
    
    
    
    # Analysis Scenario 2
    # edit the data frame and delete the patients 
    # who were allocated to treatment 2
    dataT1C <- data %>%
      filter(!TreatmentArm %in% 'T2')

    # fit one-way anova model
    fit4 <- lm(Outcome ~ TreatmentArm,
               data = dataT1C)

    # use lsmeans to compare the groups
    fit4means <- lsmeans(fit4, 
                         "TreatmentArm")

    # calculate treatment vs control comparison
    contrasts4TC <- contrast(object = fit4means, 
                             method = "trt.vs.ctrl",
                             adjust = "none")
    # calculate the contrasts one-sided
    contrasts4TC <- test(contrasts4TC, side = ">")
    # transform emmsgrid into data frame to extract values
    contrasts4TC <- data.frame(contrasts4TC)  
      
    # extract the estimate
    estimateCT1FitAllC <- contrasts4TC[[2]][1]
    # extract the p-value
    pValCT1FitAllC <- contrasts4TC[[6]][1]
    
    # calculate confidence interval for treatment versus control contrast
    contrasts4TCCI <- confint(contrast(object = fit4means, 
                                       method = "trt.vs.ctrl",
                                       adjust = "none"))    
    # extract upper and lower bounds of CI
    CICT1FitAllC<- c(contrasts4TCCI[[5]][1],
                     contrasts4TCCI[[6]][1])
    
      
    # Analysis Scenario 2  
    # edit the data frame and delete the patients 
    # who were allocated to treatment 1
    dataT2C <- data %>%
      filter(!TreatmentArm %in% 'T1')
    
    # fit one-way anova model
    fit5 <- lm(Outcome ~ TreatmentArm,
               data = dataT2C)
    
    # use lsmeans to compare the groups
    fit5means <- lsmeans(fit5, 
                         "TreatmentArm")

    # calculate treatment vs control comparison
    contrasts5TC <- contrast(object = fit5means, 
                             method = "trt.vs.ctrl",
                             adjust = "none")
    # calculate the contrasts one-sided
    contrasts5TC <- test(contrasts5TC, side = ">")
    # transform emmsgrid into data frame to extract values
    contrasts5TC <- data.frame(contrasts5TC)
    
    # extract the estimate
    estimateCT2FitAllC <- contrasts5TC[[2]][1]
    # extract the p-value
    pValCT2FitAllC <- contrasts5TC[[6]][1]  
      
    # calculate confidence interval for treatment versus control contrast
    contrasts5TCCI <- confint(contrast(object = fit5means, 
                                       method = "trt.vs.ctrl",
                                       adjust = "none"))    
    # extract upper and lower bounds of CI
    CICT2FitAllC<- c(contrasts5TCCI[[5]][1],
                     contrasts5TCCI[[6]][1])
    
    
    
    
    
    # Analysis Scenario 4
    # fit two-way anova model for all control data and view summary of model
    fitTwo <- lm(Outcome ~ Typeofgroups + TreatmentArm,
                 data = data)

    # use lsmeans to compare the groups and view summary
    fitTwoMeans <- lsmeans(fitTwo,
                           "TreatmentArm")

    #calculate treatment vs control comparison
    contrastsTCTwo <- contrast(object = fitTwoMeans, 
                               method = "trt.vs.ctrl",
                               adjust = "none")
    # calculate the contrasts one-sided
    contrastsTCTwo <- test(contrastsTCTwo, side = ">")
    # transform it into a data frame to be able to extract the values
    contrastsTCTwo <- data.frame(contrastsTCTwo)
    
    # extract the estimates
    estimateCT1Two <- contrastsTCTwo[[2]][1]
    estimateCT2Two <- contrastsTCTwo[[2]][2]
    # extract the p-values
    pValCT1Two <- contrastsTCTwo[[6]][1]
    pValCT2Two <- contrastsTCTwo[[6]][2]
    
    # calculate confidence interval for treatment vs control comparison
    contrastsTCTwoCI <- confint(contrast(object = fitTwoMeans, 
                                         method = "trt.vs.ctrl",
                                         adjust = "none"))
    # transform CI contrasts into data frame to extract upper and lower bounds
    contrastsTCTwoCI <- data.frame(contrastsTCTwoCI)
    # extract upper and lower bounds
    CICT1Two <- c(contrastsTCTwoCI[[5]][1],
                  contrastsTCTwoCI[[6]][1])
    CICT2Two <- c(contrastsTCTwoCI[[5]][2],
                  contrastsTCTwoCI[[6]][2])
    
    
    # Analysis Scenario 6  
    # fit two-way anova model for only the control data of groups X and Z
    # and view summary of model
    fitTwoXZ <- lm(Outcome ~ Typeofgroups + TreatmentArm,
                   data = dataXZ)

    # use lsmeans to compare the groups and view summary
    fitTwoXZMeans <- lsmeans(fitTwoXZ,
                             "TreatmentArm")

    # calculate treatment vs control comparison
    contrastsTCTwoXZ <- contrast(object = fitTwoXZMeans, 
                                 method = "trt.vs.ctrl",
                                 adjust = "none")
    # calculate the contrasts one-sided
    contrastsTCTwoXZ <- test(contrastsTCTwoXZ, side = ">")
    # transform it into a data frame to be able to extract the values
    contrastsTCTwoXZ <- data.frame(contrastsTCTwoXZ)
    
    # extract the estimates
    estimateCT1XZ <- contrastsTCTwoXZ[[2]][1]
    # extract the p-values
    pValCT1XZ <- contrastsTCTwoXZ[[6]][1]
    
    # calculate confidence interval for treatment vs control comparison
    contrastsTCTwoXZCI <- confint(contrast(object = fitTwoXZMeans, 
                                           method = "trt.vs.ctrl",
                                           adjust = "none"))
    # transform it into a data frame to be able to extract the values
    contrastsTCTwoXZCI <- data.frame(contrastsTCTwoXZCI) 
    # get upper and lower bounds for contrasts of interest
    CICT1XZ <- c(contrastsTCTwoXZCI[[5]][1],
                 contrastsTCTwoXZCI[[6]][1])
    
    
    # Analysis Scenario 6
    # fit two-way anova model for only the control data of groups Y and Z
    # and view summary of model
    fitTwoYZ <- lm(Outcome ~ Typeofgroups + TreatmentArm,
                   data = dataYZ)

    # use lsmeans to compare the groups and view summary
    fitTwoYZMeans <- lsmeans(fitTwoYZ,
                             "TreatmentArm")

    # calculate treatment vs control comparison
    contrastsTCTwoYZ <- contrast(object = fitTwoYZMeans, 
                                 method = "trt.vs.ctrl",
                                 adjust = "none")
    # calculate the contrasts one-sided
    contrastsTCTwoYZ <- test(contrastsTCTwoYZ, side = ">")
    # transform it into a data frame to be able to extract the values
    contrastsTCTwoYZ <- data.frame(contrastsTCTwoYZ)
    
    # extract the estimates
    estimateCT2YZ <- contrastsTCTwoYZ[[2]][1]
    # extract the p-values
    pValCT2YZ <- contrastsTCTwoYZ[[6]][1]  

    # calculate confidence interval for treatment vs control comparison
    contrastsTCTwoYZCI <- confint(contrast(object = fitTwoYZMeans, 
                                           method = "trt.vs.ctrl",
                                           adjust = "none"))
    # transform it into a data frame to be able to extract the values
    contrastsTCTwoYZCI <- data.frame(contrastsTCTwoYZCI)
    # get upper and lower bounds for contrasts of interest
    CICT2YZ <- c(contrastsTCTwoYZCI[[5]][1],
                 contrastsTCTwoYZCI[[6]][1])
    
    
    # Analysis Scenario 5
    # fit two-way anova model for all control data of the data set which
    # only includes T1 and C and view summary of model
    fitTwoT1C <- lm(Outcome ~ Typeofgroups + TreatmentArm,
                    data = dataT1C)

    # use lsmeans to compare the groups and view summary
    fitTwoT1CMeans <- lsmeans(fitTwoT1C,
                              "TreatmentArm")

    # calculate treatment vs control comparison
    contrastsTCTwoT1C <- contrast(object = fitTwoT1CMeans, 
                                  method = "trt.vs.ctrl",
                                  adjust = "none")
    # calculate the contrasts one-sided
    contrastsTCTwoT1C <- test(contrastsTCTwoT1C, side = ">")
    # transform it into a data frame to be able to extract the values
    contrastsTCTwoT1C <- data.frame(contrastsTCTwoT1C)  
      
    # extract the estimates
    estimateCT1AllC <- contrastsTCTwoT1C[[2]][1]
    # extract the p-values
    pValCT1AllC <- contrastsTCTwoT1C[[6]][1]
    
    # calculate confidence interval for treatment vs control comparison
    contrastsTCTwoT1CCI <- confint(contrast(object = fitTwoT1CMeans, 
                                            method = "trt.vs.ctrl",
                                            adjust = "none"))
    # transform it into a data frame to be able to extract the values
    contrastsTCTwoT1CCI <- data.frame(contrastsTCTwoT1CCI)
    # get upper and lower bounds for contrasts of interest
    CICT1AllC <- c(contrastsTCTwoT1CCI[[5]][1],
                   contrastsTCTwoT1CCI[[6]][1])
    
    
    # Analysis Scenario 5
    # fit two-way anova model for all control data of the data set which
    # only includes T2 and C and view summary of model
    fitTwoT2C <- lm(Outcome ~ Typeofgroups + TreatmentArm,
                    data = dataT2C)
    
    # use lsmeans to compare the groups and view summary
    fitTwoT2CMeans <- lsmeans(fitTwoT2C,
                              "TreatmentArm")

    # calculate treatment vs control comparison
    contrastsTCTwoT2C <- contrast(object = fitTwoT2CMeans, 
                                  method = "trt.vs.ctrl",
                                  adjust = "none")
    # calculate the contrasts one-sided
    contrastsTCTwoT2C <- test(contrastsTCTwoT2C, side = ">")
    # transform it into a data frame to be able to extract the values
    contrastsTCTwoT2C <- data.frame(contrastsTCTwoT2C)
    
    # extract the estimates
    estimateCT2AllC <- contrastsTCTwoT2C[[2]][1]
    # extract the p-values
    pValCT2AllC <- contrastsTCTwoT2C[[6]][1]
    
    # calculate confidence interval for treatment vs control comparison
    
    contrastsTCTwoT2CCI <- confint(contrast(object = fitTwoT2CMeans, 
                                            method = "trt.vs.ctrl",
                                            adjust = "none"))
    # transform it into a data frame to be able to extract the values
    contrastsTCTwoT2CCI <- data.frame(contrastsTCTwoT2CCI)
    # get upper and lower bounds for contrasts of interest
    CICT2AllC <- c(contrastsTCTwoT2CCI[[5]][1],
                   contrastsTCTwoT2CCI[[6]][1])
    
    # optionally give back the dataset and the input parameters
    # otherwise just return the output parameters which are necessary 
    # for the next function
    if(full == TRUE) { 
      
    # return a list containing the data frame and parameters and the operating
    # characteristics and analysis in a sublist
    results_list <- list(Dataset = data,
                         Timepoints = t,
                         SampleSize = sampleSize,
                         NumbArms = armsNumb_vec,
                         Prevalence = allocProb_vec, 
                         MeanCX = meanCX,
                         MeanCY = meanCY,
                         MeanCZ = meanCZ,
                         MeanT1X = meanT1X, 
                         MeanT2Y = meanT2Y, 
                         MeanT1Z = meanT1Z,
                         MeanT2Z = meanT2Z,
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
                         list2 = list(PatientsperArm = patientsPerArm,
                                      TrtvsCrtlAll = contrastsTC,
                                      estimateCT1All = estimateCT1All,
                                      estimateCT2All = estimateCT2All,
                                      PValueControlT1All = pValCT1All,
                                      PValueControlT2All = pValCT2All,
                                      CITrtvsCrtlAll = contrastsTCCI,
                                      CICT1All = CICT1All,
                                      CICT2All = CICT2All,
                                      TrtvsCrtlCT1 = contrasts2TC,
                                      estimateCT1Fit = estimateCT1Fit,
                                      pValCT1Fit = pValCT1Fit,
                                      CITrtvsCrtlCT1Fit = contrasts2TCCI,
                                      CICT1Fit = CICT1Fit,
                                      TrtvsCrtlCT2 = contrasts3TC,
                                      estimateCT2Fit = estimateCT2Fit,
                                      pValCT2Fit = pValCT2Fit,
                                      CITrtvsCrtlCT2Fit = contrasts3TCCI,
                                      CICT2Fit = CICT2Fit,
                                      TrtvsCrtlCT1AllC = contrasts4TC,
                                      estimateCT1FitAllC = estimateCT1FitAllC,
                                      pValCT1FitAllC = pValCT1FitAllC,
                                      CITrtvsCrtlCT1FitAllC = contrasts4TCCI,
                                      CICT1FitAllC = CICT1FitAllC,
                                      TrtvsCrtlCT2AllC = contrasts5TC,
                                      estimateCT2FitAllC = estimateCT2FitAllC,
                                      pValCT2FitAllC = pValCT2FitAllC,
                                      CITrtvsCrtlCT2FitAllC = contrasts5TCCI,
                                      CICT2FitAllC = CICT2FitAllC,
                                      TrtvsCrtlTwo = contrastsTCTwo,
                                      estimateCT1Two = estimateCT1Two,
                                      estimateCT2Two = estimateCT2Two,  
                                      pValCT1Two = pValCT1Two,
                                      pValCT2Two = pValCT2Two,
                                      CITrtvsCrtlTwo = contrastsTCTwoCI,
                                      CICT1Two = CICT1Two,
                                      CICT2Two = CICT2Two,
                                      TrtvsCrtlTwoXZ = contrastsTCTwoXZ,
                                      estimateCT1XZ = estimateCT1XZ,
                                      pValCT1XZ = pValCT1XZ,
                                      CITrtvsCrtlTwoXZ = contrastsTCTwoXZCI,
                                      CICT1XZ = CICT1XZ,
                                      TrtvsCrtlTwoYZ = contrastsTCTwoYZ,
                                      estimateCT2YZ = estimateCT2YZ,
                                      pValCT2YZ = pValCT2YZ,
                                      CITrtvsCrtlTwoYZ = contrastsTCTwoYZCI,
                                      CICT2YZ = CICT2YZ,
                                      TrtvsCrtlTwoT1C = contrastsTCTwoT1C,
                                      estimateCT1AllC = estimateCT1AllC,
                                      pValCT1AllC = pValCT1AllC,
                                      CITrtvsCrtlTwoT1C = contrastsTCTwoT1CCI,
                                      CICT1AllC = CICT1AllC,
                                      TrtvsCrtlTwoT2C = contrastsTCTwoT2C,
                                      estimateCT2AllC = estimateCT2AllC,
                                      pValCT2AllC = pValCT2AllC,
                                      TrtvsCrtlTwoT2CCI = contrastsTCTwoT2CCI,
                                      CICT2AllC = CICT2AllC,
                                      AllocationList = listAlloc_vec
                           ))
    return(results_list) 
    
    } else {
      
      # return a sublist only containing the operating characteristics and analysis 
      results_list <- list(list2 = list(PatientsperArm = patientsPerArm,
                                        estimateCT1All = estimateCT1All,
                                        estimateCT2All = estimateCT2All,
                                        PValueControlT1All = pValCT1All,
                                        PValueControlT2All = pValCT2All,
                                        CICT1All = CICT1All,
                                        CICT2All = CICT2All,
                                        estimateCT1Fit = estimateCT1Fit,
                                        pValCT1Fit = pValCT1Fit,
                                        CICT1Fit = CICT1Fit,
                                        estimateCT2Fit = estimateCT2Fit,
                                        pValCT2Fit = pValCT2Fit,
                                        CICT2Fit = CICT2Fit,
                                        estimateCT1FitAllC = estimateCT1FitAllC,
                                        pValCT1FitAllC = pValCT1FitAllC,
                                        CICT1FitAllC = CICT1FitAllC,
                                        estimateCT2FitAllC = estimateCT2FitAllC,
                                        pValCT2FitAllC = pValCT2FitAllC,
                                        CICT2FitAllC = CICT2FitAllC,
                                        estimateCT1Two = estimateCT1Two,
                                        estimateCT2Two = estimateCT2Two,  
                                        pValCT1Two = pValCT1Two,
                                        pValCT2Two = pValCT2Two,
                                        CICT1Two = CICT1Two,
                                        CICT2Two = CICT2Two,
                                        estimateCT1XZ = estimateCT1XZ,
                                        pValCT1XZ = pValCT1XZ,
                                        CICT1XZ = CICT1XZ,
                                        estimateCT2YZ = estimateCT2YZ,
                                        pValCT2YZ = pValCT2YZ,
                                        CICT2YZ = CICT2YZ,
                                        estimateCT1AllC = estimateCT1AllC,
                                        pValCT1AllC = pValCT1AllC,
                                        CICT1AllC = CICT1AllC,
                                        estimateCT2AllC = estimateCT2AllC,
                                        pValCT2AllC = pValCT2AllC,
                                        CICT2AllC = CICT2AllC,
                                        AllocationList = listAlloc_vec
                           ))
      return(results_list)
    }
}

##########################################################################################
# Example to call the function for the following setting:
# overall ratio of 1:1:2 with block randomization, random allocation to the different 
# types of groups, so within the substudies a 1:1 ratio for treatment(s) versus control

data <- fnSimulation(t = 1,
                     n = 300,
                     deltaT1 = 0.5,
                     deltaT2 = 0.5,
                     meanCX = 0,
                     meanCY = 0,
                     meanCZ = 0,
                     allocProb_vec = c((1/3), (1/3), (1/3)),
                     allocProbRandXY_vec = c((1/2), (1/2)),
                     allocProbRandZ_vec = c((1/4), (1/4), (1/4), (1/4)),
                     blocksizeXY = 2,
                     blocksizeZ = 4,
                     complete = FALSE,
                     random = TRUE,
                     full = TRUE)
