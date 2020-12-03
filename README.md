# JNA_cSDH

This repository contains code used in the analysis reported in "Development of pre and postoperative risk stratification tools to identify postoperative complications in patients undergoing surgery for chronic subdural haematoma." doi:........

For patient confidentiality reasons the source data cannot be made available.

sdh_modelbuilding_JNA.R = This contains the code used in tidying analytical dataset, multiple imputation, and univariable screening across multiply imputed datasets.  It also demonstrates the stepwise regression used to form the final multivariable models.  

sdh_cv_loops_JNA.R = This contains the code used in internally cross-validating the models generated in the other file.  Here, a 'fold then impute' strategy is used to minimise the risk of bias from combining internal validation approaches with multiple imputation.  The code also generates the ROC curves included in the published supplementary material and the calculation of their SD.  

