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
#	  
#	  I suspect we will be able to loop over C. Then we just have to 
#	  submit 2 separate jobs on NecTAR


rm(list=ls())
#####################
##  Dependencies
source('R/functions-Analyses.R')



####################################
#  Additive effects (hf = hm = 0.5)

###########
# 1: C = 0

r.vals      <-  c(0.5, 0.2, 0.1, 0)
sm.vals     <-  runif(10)
sf.vals     <-  runif(10)
Poly     <-  0
eigPoly  <-  0
agree    <-  0


Fii.init    <-  c(0.99,0,0,0,0,0,0,0,0,0.01)

	for (i in 1:length(r.vals)) {
		for (j in 1:length(sm.vals)) {
			for (k in 1:length(sf.vals)) {
				
				par.list  <-  list(
								   gen  =  5000,
								   C    =  0,
								   sm   =  sm.vals[j],
								   sf   =  sf.vals[k],
								   hm   =  0.5,
								   hf   =  0.5,
								   rm   =  r.vals[i],
								   rf   =  r.vals[i]
								  )

				res      <-  twoLocusSAPartSelfRecSim(par.list = par.list, Fii.init = Fii.init, threshold = 1e-7)
				Poly     <-  c(Poly, res$Poly)
				eigPoly  <-  c(eigPoly, res$eigPoly)
				agree    <-  c(agree, res$agree)
			
			}
		} 
		print(i)
	}

Poly     <-  Poly[-1]
eigPoly  <-  eigPoly[-1]
agree    <-  agree[-1]
sum(agree)/length(agree)

results.df  <-  data.frame("hf"      = rep(0.5, length(r.vals)*length(sm.vals)),
						   "hm"      = rep(0.5, length(r.vals)*length(sm.vals)),
						   "C"       = rep(0,   length(r.vals)*length(sm.vals)),
						   "r"       = c(rep(r.vals[1],length(sm.vals)), 
						   		  		 rep(r.vals[2],length(sm.vals)),
						   		  		 rep(r.vals[3],length(sm.vals)),
						   		  		 rep(r.vals[4],length(sm.vals))),
						   "sm"      = sm.vals
						   "Poly"    = Poly,
						   "eigPoly" = eigPoly,
						   "agree"   = agree
						   )



head(results.df)








