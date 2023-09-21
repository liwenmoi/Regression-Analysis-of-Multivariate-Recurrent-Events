# Regression-Analysis-of-Multivariate-Recurrent-Events
Executable codes of the manuscript:

"Regression Analysis of Multivariate Recurrent Event Data Allowing Time-Varying Dependence with Application to Stoke Registry Data"

1. main.R: R file contains codes to implement the proposed method. It shows how to analyze the example datasets stored in "dat.RData".

2. HelperFunctions.R: helper functions that are used in main.R.
  
3. dat.RData: This RData includes 3 datasets. The sample size is 200.

	1) DATA.1, the dataset for the 1st recurrent event. Each subject can have multiple rows, depending on the number of recurrent events of this subject. The following columns are included:

                                 (1) id1: ID number

                                 (2) t1: observed recurrent event time or censoring time. If a subject is censored before experiencing any recurrent events, then the number is identical to that in column t1.cen.

                                 (3) t1.cen: censoring time for the subject. 

                                 (4) cen1: recurrent event indicator. 1=recurrent event; 0=o.w.

                                 (5) num1: total number of recurrent events of the subject

	2) DATA.2, the dataset for the 2nd recurrent event. Each subject can have multiple rows, depending on the number of recurrent events of this subject. The following columns are included:

                                 (1) id2: ID number

                                 (2) t2: observed recurrent event time or censoring time. If a subject is censored before experiencing any recurrent events, then the number is identical to that in column t2.cen.

                                 (3) t2.cen: censoring time for the subject. 

                                 (4) cen2: recurrent event indicator. 1=recurrent event; 0=o.w.

                                 (5) num2: total number of recurrent events of the subject

	3) DATA1, the dataset that records the covariates. Each row represents a subject. The following columns are included:

                                 (1) ID: ID number

                                 (2) CEN: censoring time

                                 (3) X1: covariate 1

                                 (4) X2: covariate 2
