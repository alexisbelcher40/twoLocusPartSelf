#####################################################
#  2-locus SA with partial selfing
#
#  R code to run deterministic simulations
#  of the genotypic frequency recursions
#
#  Appendix XXX:
#  
#  Article title goes here...
#  
#  Authors: Colin Olito, Crispin Jordan, Tim Connallon
#
#  NOTES:  
#          
#	  Deterministic simulation need to explore essentially the same 
#	  parameter conditions we will be presenting in Fig.1.  So, this
#	  basically means exploring sm x sf = [0,1] for 
#	  
#	  	-- hf = hm = (1/2, 1/4)
#	  	-- rf = rm = (0.5, 0.1, 0.2, 0.5)
#	  	-- C  =  (0, 0.25, 0.5, 0.75)



rm(list=ls())
#####################
##  Dependencies
source('R/functions-analyses.R')



######################################
#  Additive effects (hf = hm = 0.5)

##############
# 1: C = 0
recursionFwdSimLoop(n = 10000, gen = 5000, C = 0, hf = 0.5, hm = 0.5, r.vals = c(0.5, 0.2, 0.1, 0), threshold = 1e-7)

##############
# 2: C = 0.25
recursionFwdSimLoop(n = 10000, gen = 5000, C = 0.25, hf = 0.5, hm = 0.5, r.vals = c(0.5, 0.2, 0.1, 0), threshold = 1e-7)

##############
# 3: C = 0.5
recursionFwdSimLoop(n = 10000, gen = 5000, C = 0.5, hf = 0.5, hm = 0.5, r.vals = c(0.5, 0.2, 0.1, 0), threshold = 1e-7)

##############
# 4: C = 0.75
recursionFwdSimLoop(n = 10000, gen = 5000, C = 0.75, hf = 0.5, hm = 0.5, r.vals = c(0.5, 0.2, 0.1, 0), threshold = 1e-7)





######################################
#  Dominance reversal (hf = hm = 0.25)

##############
# 5: C = 0
recursionFwdSimLoop(n = 10000, gen = 5000, C = 0, hf = 0.25, hm = 0.25, r.vals = c(0.5, 0.2, 0.1, 0), threshold = 1e-7)

##############
# 6: C = 0.25
recursionFwdSimLoop(n = 10000, gen = 5000, C = 0.25, hf = 0.25, hm = 0.25, r.vals = c(0.5, 0.2, 0.1, 0), threshold = 1e-7)

##############
# 7: C = 0.5
recursionFwdSimLoop(n = 10000, gen = 5000, C = 0.5, hf = 0.25, hm = 0.25, r.vals = c(0.5, 0.2, 0.1, 0), threshold = 1e-7)

##############
# 8: C = 0.75
recursionFwdSimLoop(n = 10000, gen = 5000, C = 0.75, hf = 0.25, hm = 0.25, r.vals = c(0.5, 0.2, 0.1, 0), threshold = 1e-7)
