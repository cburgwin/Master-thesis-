# **Design issues in multi-arm trials**

This repository contains the main code and R functions to reproduce the results of the simulation study presented in my [master thesis](https://zenodo.org/records/10533916) which has the title "Design issues in multi-arm trials".  

## Abstract

Multi-arm trials enhance drug development by offering increased flexibility and efficiency compared to traditional randomized clinical trials. The treatment efficacy in multi-arm trials is often assessed by comparing multiple treatment arms against a shared control arm. In traditional multi-arm studies, it is generally necessary for all enrolled patients to be eligible for all treatments in the trial. However, there are situations where this requirement may not be feasible, for example, treatments may not be available at all study centres. Selective exclusion of treatment arms can be considered as a solution in such cases, allowing clinicians and patients to exclude an unsuitable treatment arm. It is important to carefully consider the implications of the selective exclusion on the overall design and analysis of the trial. A key issue is which control data can be utilized and how to address the selective subgroups. To utilize different patient populations, it has been suggested to use randomization procedures in the trial which are capable of randomizing the patients between limited subsets of interventions according to the patient background, patients preference or treatment options at the study centres.

This thesis aims to enhance the methodology for optimally utilizing different patient subgroups in the analysis of a multi-arm trial while considering different compositions of control data. It was of further interest to evaluate the analysis strategy, the randomization strategy and the distribution of the different patient populations. The results from a simulation study are presented, where the performance of the proposed approaches in terms of the type I error rate and statistical power was evaluated under a wide range of scenarios. The results from the simulation study indicate that the analyses where the different patient subgroups were
not adjusted for can lead to a substantial power loss and type I error inflation as bias in the effect estimates is introduced. Therefore, the usage of adjusted analyses
is recommended. Furthermore, the results show that the preferred randomization strategy implies an equal ratio for treatment arm vs. control arm within each subgroup. Regarding the different composition of control data, the preferred analysis strategy is to base the comparisons on the patients who could have been directly randomized to one of the arms which are of interest for the comparison.

## Functions

### Function for generating the data set
The function for generating the data set includes the assignment to the three different subgroups X, Y and Z, the randomization to the three arms T1, T2 and C, the splitting of the data sets into all data, all control data and restricted and the adjusted and unadjusted analyses.

### Function for parallelization and calculating the operating characteristics
In the second function the estimates and p-values over the iterations are collected and the operating characteristics such as power, bias and confidence intervals are calculated. Besides, the iterations are parallelized to accelerate the simulating process.

### Function  for iterating over the simulation parameters
In the third function a grid is build which consists of all simulation parameters of interest. Then, for each row of the grid one iterates over the first two functions. In the end, the data are saved as csv file and one row in the csv file corresponds to one row in the grid, i.e. one simulation scenario.

## Funding
The work for this thesis was carried out at the Center for Medical Data Science of the Medical University of Vienna. This work was part of the EU-Pearl (EU
Patient-cEntric clinicAl tRial pLatforms) project which has received funding from the Innovative Medicines Initiative (IMI) 2 Joint Untertaking under grant agreement No 853966. This Joint Untertaking received support from the European Union´s Horizon 2020 research and innovation program and EFPIA and Children´s Tumor Foundation, Global Alliance for TB Drug Development non-profit organization, Springworks Therapeutcs Inc. This thesis reflects the author‘s point of view. Neither IMI nor the European Union, EFPIA, or any Associated Partners are responsible for any use that may be made of the information contained herein.
