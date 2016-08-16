#####################################################
#  2-locus SA with partial selfing
#
#  R code to generate output data for Fig.2
#
#
#  NOTES:  
#          
#	  
#	  	-- hf = hm = (1/2, 1/4)
#	  	-- C  =  (0, 0.25, 0.5, 0.75, 1.0)
rm(list=ls())


#####################
##  Dependencies
source('R/functions-analyses.R')


######################################
#  Additive effects (hf = hm = 0.5)

propPrPFast(n = 30000, C = 0.0, hf = 0.5, hm = 0.5); print("1/10 for Add. Eff.");
propPrPFast(n = 30000, C = 0.1, hf = 0.5, hm = 0.5); print("2/10 for Add. Eff.");
propPrPFast(n = 30000, C = 0.2, hf = 0.5, hm = 0.5); print("3/10 for Add. Eff.");
propPrPFast(n = 30000, C = 0.3, hf = 0.5, hm = 0.5); print("4/10 for Add. Eff.");
propPrPFast(n = 30000, C = 0.4, hf = 0.5, hm = 0.5); print("5/10 for Add. Eff.");
propPrPFast(n = 30000, C = 0.5, hf = 0.5, hm = 0.5); print("6/10 for Add. Eff.");
propPrPFast(n = 30000, C = 0.6, hf = 0.5, hm = 0.5); print("7/10 for Add. Eff.");
propPrPFast(n = 30000, C = 0.7, hf = 0.5, hm = 0.5); print("8/10 for Add. Eff.");
propPrPFast(n = 30000, C = 0.8, hf = 0.5, hm = 0.5); print("9/10 for Add. Eff.");
propPrPFast(n = 30000, C = 0.9, hf = 0.5, hm = 0.5); print("10/10 for Add. Eff.");


######################################
#  Dominance reversal (hf = hm = 0.25)

propPrPFast(n = 30000, C = 0.0, hf = 0.25, hm = 0.25); print("1/10 for Dom. Rev. Eff.");
propPrPFast(n = 30000, C = 0.1, hf = 0.25, hm = 0.25); print("2/10 for Dom. Rev. Eff.");
propPrPFast(n = 30000, C = 0.2, hf = 0.25, hm = 0.25); print("3/10 for Dom. Rev. Eff.");
propPrPFast(n = 30000, C = 0.3, hf = 0.25, hm = 0.25); print("4/10 for Dom. Rev. Eff.");
propPrPFast(n = 30000, C = 0.4, hf = 0.25, hm = 0.25); print("5/10 for Dom. Rev. Eff.");
propPrPFast(n = 30000, C = 0.5, hf = 0.25, hm = 0.25); print("6/10 for Dom. Rev. Eff.");
propPrPFast(n = 30000, C = 0.6, hf = 0.25, hm = 0.25); print("7/10 for Dom. Rev. Eff.");
propPrPFast(n = 30000, C = 0.7, hf = 0.25, hm = 0.25); print("8/10 for Dom. Rev. Eff.");
propPrPFast(n = 30000, C = 0.8, hf = 0.25, hm = 0.25); print("9/10 for Dom. Rev. Eff.");
propPrPFast(n = 30000, C = 0.9, hf = 0.25, hm = 0.25); print("10/10 for Dom. Rev. Eff.");
