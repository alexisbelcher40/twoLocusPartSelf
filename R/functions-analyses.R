#####################################################
#  2-locus SA with partial selfing
#
#  Necessary functions for analyses of the model
#  of the genotypic frequency recursions
#
#  Appendix XXX:
#  
#  Article title goes here...
#  
#  
#  Authors: Colin Olito, Crispin Jordan, Tim Connallon
#
#  NOTES:  
#          


##########################
##  Key to Fii subscripts
#  F11  =  Fii[1]
#  F12  =  Fii[2]
#  F13  =  Fii[3]
#  F14  =  Fii[4]
#  F22  =  Fii[5]
#  F23  =  Fii[6]
#  F24  =  Fii[7]
#  F33  =  Fii[8]
#  F34  =  Fii[9]
#  F44  =  Fii[10]

#########################################
## Average fitness through each sex role

Wf.av  <-  function(Fii, Wf.mat, ...){
   (Fii[1]*Wf.mat[1,1]) + (Fii[2]*Wf.mat[1,2]) + (Fii[3]*Wf.mat[1,3]) + (Fii[4]*Wf.mat[1,4]) + (Fii[5]*Wf.mat[2,2]) + (Fii[6]*Wf.mat[2,3]) + (Fii[7]*Wf.mat[2,4]) + (Fii[8]*Wf.mat[3,3]) + (Fii[9]*Wf.mat[3,4]) + (Fii[10]*Wf.mat[4,4])
}
Wm.av  <-  function(Fii, Wm.mat, ...){
   (Fii[1]*Wm.mat[1,1]) + (Fii[2]*Wm.mat[1,2]) + (Fii[3]*Wm.mat[1,3]) + (Fii[4]*Wm.mat[1,4]) + (Fii[5]*Wm.mat[2,2]) + (Fii[6]*Wm.mat[2,3]) + (Fii[7]*Wm.mat[2,4]) + (Fii[8]*Wm.mat[3,3]) + (Fii[9]*Wm.mat[3,4]) + (Fii[10]*Wm.mat[4,4])
}


#########################################
## Haplotype frequency change in gametes

#  Ovules
x1p  <-  function(Fii, Wf.mat, Wm.mat, par.list,...) {
	((2*Fii[1]*Wf.mat[1,1]) + (Fii[2]*Wf.mat[1,2]) + (Fii[3]*Wf.mat[1,3]) + (Fii[4]*Wf.mat[1,4])) / (2*Wf.av(Fii, Wf.mat)) - 
		par.list$rf*(((Fii[4]*Wf.mat[1,4]) - (Fii[6]*Wf.mat[2,3]))/(2*Wf.av(Fii, Wf.mat)))
}

x2p  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	((2*Fii[5]*Wf.mat[2,2]) + (Fii[2]*Wf.mat[1,2]) + (Fii[6]*Wf.mat[2,3]) + (Fii[7]*Wf.mat[2,4]))/(2*Wf.av(Fii, Wf.mat)) + 
		par.list$rf*(((Fii[4]*Wf.mat[1,4]) - (Fii[6]*Wf.mat[2,3]))/(2*Wf.av(Fii, Wf.mat)))
}

x3p  <-  function(Fii, Wf.mat, Wm.mat, par.list,...) {
	 ((2*Fii[8]*Wf.mat[3,3]) + (Fii[9]*Wf.mat[3,4]) + (Fii[3]*Wf.mat[1,3]) + (Fii[6]*Wf.mat[2,3]))/(2*Wf.av(Fii, Wf.mat)) + 
 		par.list$rf*(((Fii[4]*Wf.mat[1,4]) - (Fii[6]*Wf.mat[2,3]))/(2*Wf.av(Fii, Wf.mat)))
}

x4p  <-  function(Fii, Wf.mat, Wm.mat, par.list,...) {
	((2*Fii[10]*Wf.mat[4,4]) + (Fii[9]*Wf.mat[3,4]) + (Fii[4]*Wf.mat[1,4]) + (Fii[7]*Wf.mat[2,4]))/(2*Wf.av(Fii, Wf.mat)) - 
		par.list$rf*(((Fii[4]*Wf.mat[1,4]) - (Fii[6]*Wf.mat[2,3]))/(2*Wf.av(Fii, Wf.mat)))
}

#  Pollen/Sperm
y1p  <-  function(Fii, Wf.mat, Wm.mat, par.list,...) {
	((2*Fii[1]*Wm.mat[1,1]) + (Fii[2]*Wm.mat[1,2]) + (Fii[3]*Wm.mat[1,3]) + (Fii[4]*Wm.mat[1,4])) / (2*Wm.av(Fii, Wm.mat)) - 
		par.list$rm*(((Fii[4]*Wm.mat[1,4]) - (Fii[6]*Wm.mat[2,3]))/(2*Wm.av(Fii, Wm.mat)))
}

y2p  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	((2*Fii[5]*Wm.mat[2,2]) + (Fii[2]*Wm.mat[1,2]) + (Fii[6]*Wm.mat[2,3]) + (Fii[7]*Wm.mat[2,4]))/(2*Wm.av(Fii, Wm.mat)) + 
		par.list$rm*(((Fii[4]*Wm.mat[1,4]) - (Fii[6]*Wm.mat[2,3]))/(2*Wm.av(Fii, Wm.mat)))
}

y3p  <-  function(Fii, Wf.mat, Wm.mat, par.list,...) {
	 ((2*Fii[8]*Wm.mat[3,3]) + (Fii[9]*Wm.mat[3,4]) + (Fii[3]*Wm.mat[1,3]) + (Fii[6]*Wm.mat[2,3]))/(2*Wm.av(Fii, Wm.mat)) + 
 		par.list$rm*(((Fii[4]*Wm.mat[1,4]) - (Fii[6]*Wm.mat[2,3]))/(2*Wm.av(Fii, Wm.mat)))
}

y4p  <-  function(Fii, Wf.mat, Wm.mat, par.list,...) {
	((2*Fii[10]*Wm.mat[4,4]) + (Fii[9]*Wm.mat[3,4]) + (Fii[4]*Wm.mat[1,4]) + (Fii[7]*Wm.mat[2,4]))/(2*Wm.av(Fii, Wm.mat)) - 
		par.list$rm*(((Fii[4]*Wm.mat[1,4]) - (Fii[6]*Wm.mat[2,3]))/(2*Wm.av(Fii, Wm.mat)))
}


#########################################
## Genotypic frequency recursions


F11.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x1  <-  x1p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
	y1  <-  y1p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
	(1 - par.list$C)*x1*y1 + par.list$C*(Fii[1] + Fii[2]/4 + Fii[3]/4 + Fii[4]*((1 - par.list$rf)^2)/4 + Fii[6]*(par.list$rf^2)/4)
}

F12.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x1  <-  x1p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y2  <-  y2p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	x2  <-  x2p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y1  <-  y1p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*(x1*y2 + x2*y1) + par.list$C*(Fii[2]/2 + Fii[4]*par.list$rf*(1 - par.list$rf)/2 + Fii[6]*par.list$rf*(1 - par.list$rf)/2)
}

F13.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x1  <-  x1p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y3  <-  y3p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	x3  <-  x3p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y1  <-  y1p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*(x1*y3 + x3*y1) + par.list$C*(Fii[3]/2 + Fii[4]*par.list$rf*(1 - par.list$rf)/2 + Fii[6]*par.list$rf*(1 - par.list$rf)/2)
}

F14.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x1  <-  x1p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y4  <-  y4p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	x4  <-  x4p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y1  <-  y1p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*(x1*y4 + x4*y1) + par.list$C*(Fii[4]*((1 - par.list$rf)^2)/2 + Fii[6]*(par.list$rf^2)/2)
}

F22.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x2  <-  x2p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y2  <-  y2p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*x2*y2 + par.list$C*(Fii[5] + Fii[2]/4 + Fii[4]*(par.list$rf^2)/4 + Fii[6]*((1 - par.list$rf)^2)/4 + Fii[7]/4)
}

F23.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x2  <-  x2p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y3  <-  y3p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	x3  <-  x3p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y2  <-  y2p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*(x2*y3 + x3*y2) + par.list$C*(Fii[4]*(par.list$rf^2)/2 + Fii[6]*((1 - par.list$rf)^2)/2)
}

F24.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x2  <-  x2p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y4  <-  y4p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	x4  <-  x4p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y2  <-  y2p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*(x2*y4 + x4*y2) + par.list$C*(Fii[7]/2 + Fii[4]*par.list$rf*(1 - par.list$rf)/2 + Fii[6]*par.list$rf*(1 - par.list$rf)/2)
}

F33.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x3  <-  x3p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y3  <-  y3p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*x3*y3 + par.list$C*(Fii[8] + Fii[3]/4 + Fii[4]*(par.list$rf^2)/4 + Fii[6]*((1 - par.list$rf)^2)/4 + Fii[9]/4)
}

F34.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x3  <-  x3p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y4  <-  y4p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	x4  <-  x4p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y3  <-  y3p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*(x3*y4 + x4*y3) + par.list$C*(Fii[4]*par.list$rf*(1 - par.list$rf)/2 + Fii[6]*par.list$rf*(1 - par.list$rf)/2 + Fii[9]/2)
}

F44.pr  <-  function(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x4  <-  x4p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y4  <-  y4p(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*x4*y4 + par.list$C*(Fii[10] + Fii[4]*((1 - par.list$rf)^2)/4 + Fii[6]*(par.list$rf^2)/4 + Fii[7]/4 + Fii[9]/4)
}



###########################################
##  Eigenvalues from Mathematica solutions
##  using quasi-equibirium expressions for
##  genotypic frequencies

lambda.AB1  <-  function(par.list) {
	C   <-  par.list$C
	sm  <-  par.list$sm
	sf  <-  par.list$sf
	hm  <-  par.list$hm
	hf  <-  par.list$hf
	rm  <-  par.list$rm
	rf  <-  par.list$rf

	(4 - 2*C - 2*sf + 3*C*sf - C^2*sf - 2*hf*sf + 2*C^2*hf*sf + (-1 + C)*(-C + 2*(-1 + C)*hm)*(-1 + sf)*sm)/(2*(-2 + C)*(-1 + sf))
}

lambda.AB2  <-  function(par.list) {
	C   <-  par.list$C
	sm  <-  par.list$sm
	sf  <-  par.list$sf
	hm  <-  par.list$hm
	hf  <-  par.list$hf
	rm  <-  par.list$rm
	rf  <-  par.list$rf

	(1/(2*(C - 2)*(sf - 1)^2))*
	(2*C - 2*(C^2)*rm - 6*C*sf + 2*(C^2)*sf - 4* (C^2)*hf*sf + 8*C*rm*sf - 4*(C^2)*rm*sf - 8*C*hf*rm*sf + 
	 8*(C^2)*hf*rm*sf + 3*C*(sf^2) - (C^2)*(sf^2) + 2*(C^2)*(hf^2)*(sf^2) - 
	 4*C*rm*(sf^2) + 2*(C^2)*rm*(sf^2) + 4*C*(hf^2)*rm*(sf^2) - 
	 4*(C^2)*(hf^2)*rm*(sf^2) + 2*(rf + rm + 2*sf - 2) - 2*(C^2)*rf*((hf*sf - 1)^2) + 
	 2*sf*((-2)*(hf*(rf - 1) + rm) + ((hf^2)*(rf - 1) + rm - 1)*sf) - 
	 2*(C - 1)*(C + 2*hm - 2*C*hm + 2*(C - 1)*hm*rm)*((sf - 1)^2)*sm + 
	 (C - 1)*(C*(1 + 2*(hm^2)*(rm - 1)) - 2*(hm^2)*(rm - 1))*((sf - 1)^2)*(sm^2))
}

lambda.ab1  <-  function(par.list) {
	C   <-  par.list$C
	sm  <-  par.list$sm
	sf  <-  par.list$sf
	hm  <-  par.list$hm
	hf  <-  par.list$hf
	rm  <-  par.list$rm
	rf  <-  par.list$rf

	(2*hf*sf*(sm - 1) - 2*(sm + hm*sm - 2) + C*(sf*(sm - 1) - sm + 4*hm*sm - 2) + (C^2)*(sm - 2*hm*sm + sf*(2*hf + sm - 2*hf*sm - 1)))/(2*(C - 2)*(sm - 1))
}

lambda.ab2  <-  function(par.list) {
	C   <-  par.list$C
	sm  <-  par.list$sm
	sf  <-  par.list$sf
	hm  <-  par.list$hm
	hf  <-  par.list$hf
	rm  <-  par.list$rm
	rf  <-  par.list$rf

	(1/(2*(C - 2)*(sm - 1)^2))*
		(C*(2 - 2*(4*hf*rm - 1)*sf*((sm - 1)^2) + (4*(hf^2)*rm - 1)*(sf^2)*((sm - 1)^2) + (2 + 8*hm*(rm - 1) - 8*rm)*sm + (-1 - 4*(hm^2)*(rm - 1) + 4*rm - 1)*(sm^2)) + 
		2*(-2 + 2*hf*sf - (hf^2)*(sf^2) + rf*((hf*sf - 1)^2)*((sm - 1)^2) + 2*sm + 2*hm*sm - 4*hf*sf*sm + 2*(hf^2)*(sf^2)*sm - (sm^2) - (hm^2)*(sm^2) + 2*hf*sf*(sm^2) - (hf^2)*(sf^2)*(sm^2) + rm*((hm*sm - 1)^2)) -
	    (C^2)*(-2*sf + 4*hf*sf + (sf^2) - 2*(hf^2)*(sf^2) + 2*rf*((-1 + hf*sf)^2)*((-1 + sm)^2) + 2*sm - 4*hm*sm + 4*sf*sm - 8*hf*sf*sm - 2*(sf^2)*sm + 4*(hf^2)*(sf^2)*sm - (sm^2) + 2*(hm^2)*(sm^2) - 
	    2*sf*(sm^2) + 4*hf*sf*(sm^2) + (sf^2)*(sm^2) - 2*(hf^2)*(sf^2)*(sm^2) + 2*rm*(1 - 4*hf*sf*((-1 + sm)^2) + 2*(hf^2)*(sf^2)*((-1 + sm)^2) + 2*(-2 + hm)*sm - (-2 + (hm^2))*(sm^2))))

}



##############################################################
##  Invasion conditions based on Eigenvalues from Mathematica
##  from analytic results using quasi-equibirium 
##  expressions for genotypic frequencies

##  Kidwell et al. (1977) conditions for invasion of 
##  female 
invKidwell.lAB  <-  function(sm) {
	sm/(1 + sm)
}
invKidwell.lab  <-  function(sm) {
	sm/(1 - sm)
}

##  Invasion based on lambda.AB1 for female beneficial 
##  allele under obligate outcrossing
inv.lAB1.obOut  <-  function(hf, hm, sm) {
	(hm*sm)/(1 - hf + hm*sm)
}

##  Invasion based on lambda.ab1 for female beneficial 
##  allele under obligate outcrossing
inv.lab1.obOut  <-  function(hf, hm, sm) {
	((-1 + hm)*sm)/(hf*(-1 + sm))
}

##  Invasion based on lambda.AB2 (w/ recombination)
##  for female beneficial allele under obligate outcrossing
inv.lAB2.obOut  <-  function(hf, hm, sm, rf, rm) {
	(1 + hf*(-1 + rf) + hm*sm*(2 - hm*sm) + rm*((-1 + hm*sm)^2) - 
		sqrt(-((-1 + hf)^2)*(-1 + rf)*(1 + hm*sm*(2 - hm*sm) + rm*(-1 + hm*sm)^2)))/
			(1 + (hf^2)*(-1 + rf) + hm*sm*(2 - hm*sm) + rm*((-1 + hm*sm)^2))	
}

##  Invasion based on lambda.ab2 (w/ recombination)
##  for male beneficial allele under obligate outcrossing
inv.lab2.obOut  <-  function(hf, hm, sm, rf, rm) {
	(hf*(-1 + rf)*((-1 + sm)^2) + 
			sqrt(-(hf^2)*(-1 + rf)*((-1 + sm)^2)*(1 + rm + 2*(-2 + hm - hm*rm)*sm + (2 + (hm^2)*(-1 + rm))*sm^2))) / 
		((hf^2)*(-1 + rf)*((-1 + sm)^2))
}



##  Invasion based on lambda.AB1 for female beneficial 
##  allele (w/ selfing)
inv.lAB1  <-  function(hf, hm, sm, C) {
	((-1 + C)*(-C - 2*hm + 2*C*hm)*sm)/(2 + C - (C^2) - 2*hf + 2*(C^2)*hf + C*sm - (C^2)*sm + 2*hm*sm - 4*C*hm*sm + 2*(C^2)*hm*sm)
}

##  Invasion based on lambda.ab1 for female beneficial 
##  allele (w/ selfing)
inv.lab1  <-  function(hf, hm, sm, C) {
	-(((-1 + C)*(2 - C + 2*(-1 + C)*hm)*sm)/((1 + C)*(-C + 2*(-1 + C)*hf)*(-1 + sm)))
}

##  Invasion based on lambda.AB2 (w/ recombination)
##  for female beneficial allele (w/ selfing)
inv.lAB2  <-  function(hf, hm, sm, rf, rm, C) {
(-2 - C + C^2 + 2*hf - 2*C^2*hf - 2*hf*rf + 2*C^2*hf*rf - 2*rm + 
   4*C*rm - 2*C^2*rm - 4*C*hf*rm + 4*C^2*hf*rm - 2*C*sm + 2*C^2*sm - 
   4*hm*sm + 8*C*hm*sm - 4*C^2*hm*sm + 4*hm*rm*sm - 8*C*hm*rm*sm + 
   4*C^2*hm*rm*sm + C*sm^2 - C^2*sm^2 + 2*hm^2*sm^2 - 4*C*hm^2*sm^2 + 
   2*C^2*hm^2*sm^2 - 2*hm^2*rm*sm^2 + 4*C*hm^2*rm*sm^2 - 
   2*C^2*hm^2*rm*sm^2 + sqrt((C^2*(1 + 2*hf*(-1 + rf + 2*rm) + 
           2*sm - 4*hm*sm - sm^2 + 2*hm^2*sm^2 - 
           2*rm*(-1 + hm*sm)^2) - 
        2*(1 + hf*(-1 + rf) + 2*hm*sm - hm^2*sm^2 + 
           rm*(-1 + hm*sm)^2) + 
        C*(-1 + (-2 + 8*hm)*sm + (1 - 4*hm^2)*sm^2 - 
           4*rm*(hf - (-1 + hm*sm)^2)))^2 + (-1 + C)*(2*(1 + C)*rf + 
         sm*(-2*hm*(-2 + hm*sm) + C*(2 - 4*hm - sm + 2*hm^2*sm)) + 
         2*rm*((-1 + hm*sm)^2 + C*(1 + 2*hm*sm - hm^2*sm^2)))*(2*(1 +
             hf^2*(-1 + rf) + 2*hm*sm - hm^2*sm^2 + 
            rm*(-1 + hm*sm)^2) + 
         C^2*(-1 - 2*hf^2*(-1 + rf + 2*rm) - 2*sm + 4*hm*sm + sm^2 - 
            2*hm^2*sm^2 + 2*rm*(-1 + hm*sm)^2) + 
         C*(1 + (2 - 8*hm)*sm + (-1 + 4*hm^2)*sm^2 + 
            4*rm*(hf^2 - (-1 + hm*sm)^2)))))/(C^2*(1 + 
      2*hf^2*(-1 + rf + 2*rm) + 2*sm - 4*hm*sm - sm^2 + 2*hm^2*sm^2 - 
      2*rm*(-1 + hm*sm)^2) - 
   2*(1 + hf^2*(-1 + rf) + 2*hm*sm - hm^2*sm^2 + rm*(-1 + hm*sm)^2) + 
   C*(-1 + (-2 + 8*hm)*sm + (1 - 4*hm^2)*sm^2 - 
      4*rm*(hf^2 - (-1 + hm*sm)^2)))
}


##  Invasion based on lambda.ab2 (w/ recombination)
##  for male beneficial allele (w/ selfing)
inv.lab2  <-  function(hf, hm, sm, rf, rm, C) {
	(C + (C^2) + 2*hf - 2*(C^2)*hf - 2*hf*rf + 2*(C^2)*hf*rf - 4*C*hf*rm + 
	 4*(C^2)*hf*rm - 2*C*sm - 2*(C^2)*sm - 4*hf*sm + 4*(C^2)*hf*sm + 
	 4*hf*rf*sm - 4*(C^2)*hf*rf*sm + 8*C*hf*rm*sm - 8*(C^2)*hf*rm*sm + 
	 C*(sm^2) + (C^2)*(sm^2) + 2*hf*(sm^2) - 2*(C^2)*hf*(sm^2) - 2*hf*rf*(sm^2) + 
	 2*(C^2)*hf*rf*(sm^2) - 4*C*hf*rm*(sm^2) + 4*(C^2)*hf*rm*(sm^2) - 
	 ((1/2)* sqrt(4*((C - 2*hf*(-1 + rf) - 4*C*hf*rm + 
	   (C^2)*(1 + 2*hf*(-1 + rf + 2*rm)))^2)*((-1 + sm)^4) - 
	  4*(-1 + C)*(C - 2*(hf^2)*(-1 + rf) - 4*C*(hf^2)*rm + 
	  (C^2)*(1 + 2*(hf^2)*(-1 + rf + 2*rm)))*((-1 + sm)^2)*(2*(1 + C)*rf*((-1 + sm)^2) + 
	  sm*(-2*(-1 + hm)*(-2 + sm + hm*sm) + C*(2 - 4*hm - sm + 2*(hm^2)*sm)) + 
	  2*rm*(((-1 + hm*sm)^2) + C*(1 + 2*(-2 + hm)*sm - (-2 + (hm^2))*(sm^2))))))) / 
		((C - 2*(hf^2)*(-1 + rf) - 4*C*(hf^2)*rm + (C^2)*(1 + 2*(hf^2)*(-1 + rf + 2*rm)))*((-1 + sm)^2))
}



######################################################
##  Analytic solutions (results based on Eigenvalues) 

#' 2-locus SA w/ partial selfing Eigenvalue invasion analysis
#'
#' @title 2-locus SA w/ partial selfing Eigenvalue invasion analysis
#' @param par.list A list with desired parameter values for the simulation with structure:
#' par.list  <-  list(
#'				   gen  =  5000,
#'				   C    =  0,
#'				   sm   =  0.7,
#'				   sf   =  0.7,
#'				   hm   =  0.5,
#'				   hf   =  0.5,
#'				   rm   =  0.5,
#'				   rf   =  0.5
#'				   )
#' @return Returns a list with each of the candidate leading eigenvalues calculated at the boundaries p,q = 0 and p,q = 1, 
#' as well as a coded categorical results indicating the outcome of the invasion analysis where
#' 0 = indeterminate
#' 1 = positive selection for the female-beneficial allele
#' 2 = positive selection for the male-beneficial allele
#' 3 = unstable internal equilibrium
#' 4 = protected polymorphism
#' @seealso 
#' @export
#' @author Colin Olito.
#' @examples
#' eigenInvAnalysis(par.list) 
eigenInvAnalysis  <-  function(par.list) {

	##  Warnings
	if(any(par.list[2:8] < 0) | any(par.list[2:8] > 1) | any(par.list[7:8] > 0.5))
		stop('The chosen parameter values fall outside of the reasonable bounds')

	##  Calculate Eigenvalues from analytic solutions 
	##  using quasi-equibirium genotypic frequencies

	l.AB1  <- lambda.AB1(par.list)
	l.AB2  <- lambda.AB2(par.list)
	l.ab1  <- lambda.ab1(par.list)
	l.ab2  <- lambda.ab2(par.list)

	
	## Categorize outcome of invasion analysis
	invKey  <-  list(
					 "1"  =  "Positive selection for female-beneficial allele",
					 "2"  =  "Positive selection for male-beneficial allele",
					 "3"  =  "Unstable internal equilibrium",
					 "4"  =  "Protected polymorphism"
					 )

	Inv   <-  0
	PrP   <-  0
	rPrP  <-  0
	# Positive selection for female-beneficial allele
		if (any(c(l.AB1, l.AB2) > 1) & l.ab1 < 1 & l.ab2 < 1 )
			Inv  <-  1

	# Positive selection for male-beneficial allele
		if (l.AB1 < 1 & l.AB2 < 1 & any(c(l.ab1, l.ab2) > 1))
			Inv  <-  2

	# Unstable internal equilibrium
		if (l.AB1 < 1 & l.AB2 <  1 & l.ab1 < 1 & l.ab2 < 1 )
			Inv  <-  3

	# Protected polymorphism
		if (any(c(l.AB1, l.AB2) > 1) & any(c(l.ab1, l.ab2) > 1 )) {
			Inv  <-  4
			PrP  <-  1			
		}
	
	##  Is polymorphism facilitated by recombination? 
	if (any(c(l.AB1,l.ab1) < 1) & l.AB2 > 1 &  l.ab2 > 1)
		rPrP  <-  1

	##  Output list
	res  <-  list(
				  "par.list"  =  par.list,
				  "l.AB1"     =  l.AB1,
				  "l.AB2"     =  l.AB2,
				  "l.ab1"     =  l.ab1,
				  "l.ab2"     =  l.ab2,
				  "InvKey"    =  invKey,
				  "Inv"       =  Inv,
				  "PrP"       =  PrP,
				  "rPrP"      =  rPrP
 				 )
	return(res)
}







########################
##  Simulation function
########################

#' Forward deterministic simulation of genotypic recursions for
#' 2-locus SA w/ partial selfing model
#'
#' @title Forward deterministic simulation of genotypic recursions
#' @param par.list A list with desired parameter values for the simulation with structure:
#' par.list  <-  list(
#'				   gen  =  5000,
#'				   C    =  0,
#'				   sm   =  0.7,
#'				   sf   =  0.7,
#'				   hm   =  0.5,
#'				   hf   =  0.5,
#'				   rm   =  0.5,
#'				   rf   =  0.5
#'				   )
#' @param Fii.init A vector of initial genotypic frequencies (must have length = 10).
#' c(0.99,0,0,0,0,0,0,0,0,0.01) for invasion of aabb into population 'fixed' for AABB.
#' c(0.01,0,0,0,0,0,0,0,0,0.99) for invasion of AABB into population 'fixed' for aabb.
#' @return Returns a list with timeseries for each genotype, equilibrium frequencies, and a numeric (0,1) for whether the 
#' equilibrium was polymorphic (with tolerance 1E-6).
#' @seealso 
#' @export
#' @author Colin Olito.
#' @examples
#' recursionFwdSim(par.list, Fii.init, threshold = 1e-6) 
recursionFwdSim  <-  function(par.list, Fii.init, threshold = 1e-6) {

	##  Warnings
	if(any(par.list[2:8] < 0) | any(par.list[2:8] > 1) | any(par.list[7:8] > 0.5))
		stop('The chosen parameter values fall outside of the reasonable bounds')

	if(any(Fii.init < 0) | any(Fii.init > 1))
		stop('The chosen initial genotypic frequencies are not valid')


	##  Fitness Matrices
	Wf.mat   <-  matrix(
						c((1-par.list$sf)^2                          , (1-par.list$sf)*(1-par.list$hf*par.list$sf), (1-par.list$sf)*(1-par.list$hf*par.list$sf) , (1-par.list$hf*par.list$sf)^2,
						  (1-par.list$sf)*(1-par.list$hf*par.list$sf), (1-par.list$sf)                            , (1-par.list$hf*par.list$sf)^2               , (1-par.list$hf*par.list$sf),
						  (1-par.list$sf)*(1-par.list$hf*par.list$sf), (1-par.list$hf*par.list$sf)^2              , (1-par.list$sf)                             , (1-par.list$hf*par.list$sf),
						  (1-par.list$hf*par.list$sf)^2              , (1-par.list$hf*par.list$sf)                , (1-par.list$hf*par.list$sf)                 , 1), 
						nrow=4, byrow=TRUE
						)

	Wm.mat   <-  matrix(
						c( 1                           , (1-par.list$hm*par.list$sm)                , (1-par.list$hm*par.list$sm)                 , (1-par.list$hm*par.list$sm)^2,
						  (1-par.list$hm*par.list$sm)  , (1-par.list$sm)                            , (1-par.list$hm*par.list$sm)^2               , (1-par.list$hm*par.list$sm)*(1-par.list$sm),
						  (1-par.list$hm*par.list$sm)  , (1-par.list$hm*par.list$sm)^2              , (1-par.list$sm)                             , (1-par.list$sm)*(1-par.list$hm*par.list$sm),
						  (1-par.list$hm*par.list$sm)^2, (1-par.list$sm)*(1-par.list$hm*par.list$sm), (1-par.list$sm)*(1-par.list$hm*par.list$sm) , (1-par.list$sm)^2), 
						nrow=4, byrow=TRUE
						)	

	##  Initilize data storage structures
	Fii.gen  <-  matrix(0, ncol=10, nrow=par.list$gen)

	##  Generation Loop
		# initialize
		Fii.gen[1,1]   <-  F11.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
		Fii.gen[1,2]   <-  F12.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
		Fii.gen[1,3]   <-  F13.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
		Fii.gen[1,4]   <-  F14.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
		Fii.gen[1,5]   <-  F22.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
		Fii.gen[1,6]   <-  F23.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
		Fii.gen[1,7]   <-  F24.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
		Fii.gen[1,8]   <-  F33.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
		Fii.gen[1,9]   <-  F34.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
		Fii.gen[1,10]  <-  F44.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)


	# Start simulation
	i      <-  2
	diffs  <-  rep(1,10)

	while (i < par.list$gen & any(diffs[diffs != 0] > threshold)) {
		Fii.gen[i,1]   <-  F11.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
		Fii.gen[i,2]   <-  F12.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
		Fii.gen[i,3]   <-  F13.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
		Fii.gen[i,4]   <-  F14.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
		Fii.gen[i,5]   <-  F22.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
		Fii.gen[i,6]   <-  F23.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
		Fii.gen[i,7]   <-  F24.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
		Fii.gen[i,8]   <-  F33.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
		Fii.gen[i,9]   <-  F34.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)
		Fii.gen[i,10]  <-  F44.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list)	
		
		diffs  <-  Fii.gen[i,] - Fii.gen[i-1,]
		i      <-  i+1
	}

	##  Is equilibrium polymorphic?
	if (any(Fii.gen[i-1,2:9] > 1e-5))
		 Poly  <- 1
	else Poly  <- 0

	##  Calculate Eigenvalues from analytic solutions 
	##  using quasi-equibirium genotypic frequencies

	l.AB1  <- lambda.AB1(par.list)
	l.AB2  <- lambda.AB2(par.list)
	l.ab1  <- lambda.ab1(par.list)
	l.ab2  <- lambda.ab2(par.list)

	if (any(c(l.AB1, l.AB2) > 1) & any(c(l.ab1, l.ab2) > 1 ))
		 eigPoly  <-  1
	else eigPoly  <-  0

	##  Does simulation result agree with Eigenvalues?

	if (Poly == eigPoly)
		 agree  <-  1
	else agree  <-  0


	##  Output list
	res  <-  list(
				  "par.list" =  par.list,
				  "Fii.gen"  =  Fii.gen[1:i-1,],
				  "EQ.freq"  =  Fii.gen[i-1,],
				  "l.AB1"    =  l.AB1,
				  "l.AB2"    =  l.AB2,
				  "l.ab1"    =  l.ab1,
				  "l.ab2"    =  l.ab2,
				  "Poly"     =  Poly,
				  "eigPoly"  =  eigPoly,
				  "agree"    =  agree
 				 )
	return(res)
}









#' Simulation loop wrapping forward deterministic simulations 
#' of genotypic recursions 
#' #'
#' @title Forward deterministic simulation of genotypic recursions.
#' @param n number of randomly generated values for sf & sm. 
#' Determines resolution with which parameter space is explored.
#' @param gen Maximum number of generations for each simulation (as in par.list).
#' @param C The fixed selfing rate (as in par.list).
#' @param hf Dominance through female expression (as in par.list).
#' @param hm Dominance through male expression (as in par.list).
#' @param r.vals Values of recombination rate to explore(as in par.list).
#' @param threshold Threshold difference between genotypic frequencies before simulation cuts off.
#' @return Returns a data frame with parameter values, a variable describing whether 
#' the final state of the simulation was polymorphic polymorphism, whether evaluating the eigenvalues
#' predicts polymorphism, and whether these two methods agree with one another.
#' @seealso `recursionFwdSim`
#' @export
#' @author Colin Olito.
#' @examples
#' recursionFwdSimLoop(n = 10000, gen = 5000, C = 0, hf = 0.5, hm = 0.5, r.vals = c(0.5, 0.2, 0.1, 0), threshold = 1e-7)

recursionFwdSimLoop  <-  function(n = 10000, gen = 5000, C = 0, hf = 0.5, hm = 0.5, r.vals = c(0.5, 0.2, 0.1, 0), threshold = 1e-7) {

	## Warnings
	if(any(c(C,hf,hm,r.vals) < 0) | any(c(C,hf,hm) > 1) | any(r.vals > 0.5))
		stop('At least one of the chosen parameter values fall outside of the reasonable bounds')

	if(threshold > 1e-7)
		stop('Carefully consider whether you want to change this threshold, 
			  as it will effect whether the simulations agree with the analytic results')

	#  initialize selection coeficients and storage structures
	sm.vals     <-  runif(n)
	sf.vals     <-  runif(n)
	Poly     <-  0
	eigPoly  <-  0
	agree    <-  0

	#  initial genotypic frequencies (can change to c(0.01,0,0,0,0,0,0,0,0,0.99))
	Fii.init    <-  c(0.99,0,0,0,0,0,0,0,0,0.01)

	##  Simulation Loop over values of r, sm, sf for fixed selfing rate (C)
	for (i in 1:length(r.vals)) {
		for (j in 1:length(sm.vals)) {
			for (k in 1:length(sf.vals)) {
				
				par.list  <-  list(
								   gen  =  gen,
								   C    =  C,
								   sm   =  sm.vals[j],
								   sf   =  sf.vals[k],
								   hm   =  hm,
								   hf   =  hf,
								   rm   =  r.vals[i],
								   rf   =  r.vals[i]
								  )

				res      <-  recursionFwdSim(par.list = par.list, Fii.init = Fii.init, threshold = threshold)
				Poly     <-  c(Poly, res$Poly)
				eigPoly  <-  c(eigPoly, res$eigPoly)
				agree    <-  c(agree, res$agree)
			
			}
		} 
	}

	#  trim leading zero from storage vectors
	Poly     <-  Poly[-1]
	eigPoly  <-  eigPoly[-1]
	agree    <-  agree[-1]

	#  Compile results as data.frame
	results.df  <-  data.frame("hf"      = rep(0.5, length(r.vals)*length(sm.vals)),
							   "hm"      = rep(0.5, length(r.vals)*length(sm.vals)),
							   "C"       = rep(0,   length(r.vals)*length(sm.vals)),
							   "r"       = c(rep(r.vals[1],length(sm.vals)), 
							   		  		 rep(r.vals[2],length(sm.vals)),
							   		  		 rep(r.vals[3],length(sm.vals)),
							   		  		 rep(r.vals[4],length(sm.vals))),
							   "sm"      = sm.vals,
							   "sf"      = sf.vals,
							   "Poly"    = Poly,
							   "eigPoly" = eigPoly,
							   "agree"   = agree
							   )

	#  Write results.df to .txt file
	filename  <-  paste("./output/simResults/recFwdSimLoop.out", "_C", C, "_hf", hf, "_hm", hm, ".txt", sep="")
	write.table(results.df, file=filename, col.names = TRUE, row.names = FALSE)

	#  Return results.df in case user wants it
	return(results.df)
}



































#' Proportion of parameter space resulting in protected polymorphism
#' determined by evaluating eigenvalues.
#'
#' @titleProportion of parameter space resulting in protected polymorphism
#' determined by evaluating eigenvalues.
#' @param n number of randomly generated values for sf & sm. 
#' Determines resolution with which parameter space is explored.
#' @param C A vector of selfing rates to explore.
#' @param hf Dominance through female expression (as in par.list).
#' @param hm Dominance through male expression (as in par.list).
#' @return Returns a data frame with ...
#' @seealso `recursionFwdSim`
#' @export
#' @author Colin Olito.
#' @examples
#' recursionFwdSimLoop(n = 10000, gen = 5000, C = 0, hf = 0.5, hm = 0.5, r.vals = c(0.5, 0.2, 0.1, 0), threshold = 1e-7)

propPrP  <-  function(n = 10000, C = 0, hf = 0.5, hm = 0.5) {

	## Warnings
	if(any(c(C,hf,hm) < 0) | any(c(C,hf,hm) > 1))
		stop('At least one of the chosen parameter values fall outside of the reasonable bounds')

	#  initialize selection coeficients and storage structures
	r.vals      <-  seq(0, 0.5, by=0.005)
	sm.vals     <-  runif(n)
	sf.vals     <-  runif(n)
	PrP     <-  c()
	rPrP    <-  c()


	##  Loop over values of r, sm, sf for fixed selfing rate (C)
	##  calculating proportion of parameter space where PrP is 
	##  predicted each time.

	for (i in 1:length(r.vals)) {
		eigPoly     <-  0
		rPoly       <-  0

		for (j in 1:length(sm.vals)) {
			for (k in 1:length(sf.vals)) {
				
				par.list  <-  list(
								   gen  =  NA,
								   C    =  C,
								   sm   =  sm.vals[j],
								   sf   =  sf.vals[k],
								   hm   =  hm,
								   hf   =  hf,
								   rm   =  r.vals[i],
								   rf   =  r.vals[i]
								  )

				res      <-  eigenInvAnalysis(par.list = par.list)
				eigPoly  <-  c(eigPoly, res$PrP)
				rPoly    <-  c(rPoly,    res$rPrP)
			
			}
		}
		#  trim leading zero from storage vectors
		eigPoly  <-  eigPoly[-1]
		rPoly    <-  rPoly[-1]
		
		#  Calculate proportion of parameter space resulting in PrP
		PrP[i]   <-  sum(eigPoly)/length(eigPoly)
		rPrP[i]  <-  sum(rPoly)/length(rPoly)

	}


	#  Compile results as data.frame
	results.df  <-  data.frame("hf"      = rep(hf, length(r.vals)),
							   "hm"      = rep(hm, length(r.vals)),
							   "C"       = rep(C,  length(r.vals)),
							   "r"       = r.vals,
							   "PrP"     = PrP,
							   "rPrP"    = rPrP
							   )

	#  Write results.df to .txt file
	filename  <-  paste("./output/data/propPrp.out", "_C", C, "_hf", hf, "_hm", hm, "_n", n, ".txt", sep="")
	write.table(results.df, file=filename, col.names = TRUE, row.names = FALSE)

	#  Return results.df in case user wants it
#	return(results.df)
}




