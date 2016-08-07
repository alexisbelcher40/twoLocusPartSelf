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
#  
#  Authors: Colin Olito, Crispin Jordan, Tim Connallon

#  NOTES:  
#          


####################################
#  Simulation parameter values

par.list  <-  list(
				   gen  =  10,
				   C    =  0,
				   sm   =  0.1,
				   sf   =  0.1,
				   hm   =  1/2,
				   hf   =  1/2,
				   rm   =  0.1,
				   rf   =  0.1
				   )


####################################
#  Fitness Expressions

Wf.mat   <-  matrix(
					c((1-par.list$sf)^2                          , (1-par.list$sf)*(1-par.list$hf*par.list$sf), (1-par.list$sf)*(1-par.list$hf*par.list$sf) , (1-par.list$hf*par.list$sf)^2,
					  (1-par.list$sf)*(1-par.list$hf*par.list$sf), (1-par.list$sf)                            , (1-par.list$hf*par.list$sf)^2               , (1-par.list$hf*par.list$sf),
					  (1-par.list$sf)*(1-par.list$hf*par.list$sf), (1-par.list$hf*par.list$sf)^2              , (1-par.list$sf)                             , (1-par.list$hf*par.list$sf),
					  (1-par.list$hf*par.list$sf)^2              , (1-par.list$hf*par.list$sf)                , (1-par.list$hf*par.list$sf)                 , 1), 
					nrow=4, byrow=TRUE
					)
Wf.mat

Wm.mat   <-  matrix(
					c( 1                           , (1-par.list$hm*par.list$sm)                , (1-par.list$hm*par.list$sm)                 , (1-par.list$hm*par.list$sm)^2,
					  (1-par.list$hm*par.list$sm)  , (1-par.list$sm)                            , (1-par.list$hm*par.list$sm)^2               , (1-par.list$hm*par.list$sm)*(1-par.list$sm),
					  (1-par.list$hm*par.list$sm)  , (1-par.list$hm*par.list$sm)^2              , (1-par.list$sm)                             , (1-par.list$sm)*(1-par.list$hm*par.list$sm),
					  (1-par.list$hm*par.list$sm)^2, (1-par.list$sm)*(1-par.list$hm*par.list$sm), (1-par.list$sm)*(1-par.list$hm*par.list$sm) , (1-par.list$sm)^2), 
					nrow=4, byrow=TRUE
					)

#  Check that fitness matrices are correct (make sure that sm == sf & hf==hm)
#  Wm.mat == Wf.mat[c(4,3,2,1),c(4,3,2,1)]

######################################
#  Initialize data storage structures



Fii.init  <-  c(0.99,0,0,0,0,0,0,0,0,0.01)



########################
########################
##  Simulation function


twoLocusSAPartSelf  <-  function(par.list, init.freq) {

	##  Warnings

	##  Initilize data storage structures
	Fii.gen  <-  matrix(0, ncol=10, nrow=par.list$gen)


	##  Generation Loop
	for (i in 1:par.list$gen) {

		if (i == 1)

			Fii.gen[i,1]  <-  F11.pr(Fii.init, Wf.mat, Wm.mat, par.list,...)

	}
}