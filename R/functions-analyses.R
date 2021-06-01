#####################################################
#  2-locus SA with partial selfing
#
#  Necessary functions for analyses of the model
#  of the genotypic frequency recursions
#
#  Consequences of genetic linkage for the maintenance 
#  of sexually antagonistic polymorphism in hermaphrodites
#  
#  Author: Colin Olito
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
x1  <-  function(Fii, Wf.mat, Wm.mat, par.list,...) {
	((2*Fii[1]*Wf.mat[1,1]) + (Fii[2]*Wf.mat[1,2]) + (Fii[3]*Wf.mat[1,3]) + (Fii[4]*Wf.mat[1,4])) / (2*Wf.av(Fii, Wf.mat)) - 
		par.list$r*(((Fii[4]*Wf.mat[1,4]) - (Fii[6]*Wf.mat[2,3]))/(2*Wf.av(Fii, Wf.mat)))
}

x2  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	((2*Fii[5]*Wf.mat[2,2]) + (Fii[2]*Wf.mat[1,2]) + (Fii[6]*Wf.mat[2,3]) + (Fii[7]*Wf.mat[2,4]))/(2*Wf.av(Fii, Wf.mat)) + 
		par.list$r*(((Fii[4]*Wf.mat[1,4]) - (Fii[6]*Wf.mat[2,3]))/(2*Wf.av(Fii, Wf.mat)))
}

x3  <-  function(Fii, Wf.mat, Wm.mat, par.list,...) {
	 ((2*Fii[8]*Wf.mat[3,3]) + (Fii[9]*Wf.mat[3,4]) + (Fii[3]*Wf.mat[1,3]) + (Fii[6]*Wf.mat[2,3]))/(2*Wf.av(Fii, Wf.mat)) + 
 		par.list$r*(((Fii[4]*Wf.mat[1,4]) - (Fii[6]*Wf.mat[2,3]))/(2*Wf.av(Fii, Wf.mat)))
}

x4  <-  function(Fii, Wf.mat, Wm.mat, par.list,...) {
	((2*Fii[10]*Wf.mat[4,4]) + (Fii[9]*Wf.mat[3,4]) + (Fii[4]*Wf.mat[1,4]) + (Fii[7]*Wf.mat[2,4]))/(2*Wf.av(Fii, Wf.mat)) - 
		par.list$r*(((Fii[4]*Wf.mat[1,4]) - (Fii[6]*Wf.mat[2,3]))/(2*Wf.av(Fii, Wf.mat)))
}

#  Pollen/Sperm
y1  <-  function(Fii, Wf.mat, Wm.mat, par.list,...) {
	((2*Fii[1]*Wm.mat[1,1]) + (Fii[2]*Wm.mat[1,2]) + (Fii[3]*Wm.mat[1,3]) + (Fii[4]*Wm.mat[1,4])) / (2*Wm.av(Fii, Wm.mat)) - 
		par.list$r*(((Fii[4]*Wm.mat[1,4]) - (Fii[6]*Wm.mat[2,3]))/(2*Wm.av(Fii, Wm.mat)))
}

y2  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	((2*Fii[5]*Wm.mat[2,2]) + (Fii[2]*Wm.mat[1,2]) + (Fii[6]*Wm.mat[2,3]) + (Fii[7]*Wm.mat[2,4]))/(2*Wm.av(Fii, Wm.mat)) + 
		par.list$r*(((Fii[4]*Wm.mat[1,4]) - (Fii[6]*Wm.mat[2,3]))/(2*Wm.av(Fii, Wm.mat)))
}

y3  <-  function(Fii, Wf.mat, Wm.mat, par.list,...) {
	 ((2*Fii[8]*Wm.mat[3,3]) + (Fii[9]*Wm.mat[3,4]) + (Fii[3]*Wm.mat[1,3]) + (Fii[6]*Wm.mat[2,3]))/(2*Wm.av(Fii, Wm.mat)) + 
 		par.list$r*(((Fii[4]*Wm.mat[1,4]) - (Fii[6]*Wm.mat[2,3]))/(2*Wm.av(Fii, Wm.mat)))
}

y4  <-  function(Fii, Wf.mat, Wm.mat, par.list,...) {
	((2*Fii[10]*Wm.mat[4,4]) + (Fii[9]*Wm.mat[3,4]) + (Fii[4]*Wm.mat[1,4]) + (Fii[7]*Wm.mat[2,4]))/(2*Wm.av(Fii, Wm.mat)) - 
		par.list$r*(((Fii[4]*Wm.mat[1,4]) - (Fii[6]*Wm.mat[2,3]))/(2*Wm.av(Fii, Wm.mat)))
}


#########################################
## Genotypic frequency recursions


F11.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x1  <-  x1(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y1  <-  y1(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*x1*y1 + par.list$C*((Fii[1]*Wf.mat[1,1] + Fii[2]*Wf.mat[1,2]/4 + Fii[3]*Wf.mat[1,3]/4 + Fii[4]*Wf.mat[1,4]*((1 - par.list$r)^2)/4 + Fii[6]*Wf.mat[2,3]*(par.list$r^2)/4)/(Wf.av(Fii=Fii, Wf.mat=Wf.mat)))
}

F12.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x1  <-  x1(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y2  <-  y2(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	x2  <-  x2(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y1  <-  y1(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*(x1*y2 + x2*y1) + par.list$C*((Fii[2]*Wf.mat[1,2]/2 + Fii[4]*Wf.mat[1,4]*par.list$r*(1 - par.list$r)/2 + Fii[6]*Wf.mat[2,3]*par.list$r*(1 - par.list$r)/2)/(Wf.av(Fii=Fii, Wf.mat=Wf.mat)))
}

F13.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x1  <-  x1(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y3  <-  y3(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	x3  <-  x3(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y1  <-  y1(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*(x1*y3 + x3*y1) + par.list$C*((Fii[3]*Wf.mat[1,3]/2 + Fii[4]*Wf.mat[1,4]*par.list$r*(1 - par.list$r)/2 + Fii[6]*Wf.mat[2,3]*par.list$r*(1 - par.list$r)/2)/(Wf.av(Fii=Fii, Wf.mat=Wf.mat)))
}

F14.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x1  <-  x1(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y4  <-  y4(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	x4  <-  x4(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y1  <-  y1(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*(x1*y4 + x4*y1) + par.list$C*((Fii[4]*Wf.mat[1,4]*((1 - par.list$r)^2)/2 + Fii[6]*Wf.mat[2,3]*(par.list$r^2)/2)/(Wf.av(Fii=Fii, Wf.mat=Wf.mat)))
}

F22.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x2  <-  x2(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y2  <-  y2(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*x2*y2 + par.list$C*((Fii[5]*Wf.mat[2,2] + Fii[2]*Wf.mat[1,2]/4 + Fii[4]*Wf.mat[1,4]*(par.list$r^2)/4 + Fii[6]*Wf.mat[2,3]*((1 - par.list$r)^2)/4 + Fii[7]*Wf.mat[2,4]/4)/(Wf.av(Fii=Fii, Wf.mat=Wf.mat)))
}

F23.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x2  <-  x2(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y3  <-  y3(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	x3  <-  x3(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y2  <-  y2(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*(x2*y3 + x3*y2) + par.list$C*((Fii[4]*Wf.mat[1,4]*(par.list$r^2)/2 + Fii[6]*Wf.mat[2,3]*((1 - par.list$r)^2)/2)/(Wf.av(Fii=Fii, Wf.mat=Wf.mat)))
}

F24.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x2  <-  x2(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y4  <-  y4(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	x4  <-  x4(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y2  <-  y2(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*(x2*y4 + x4*y2) + par.list$C*((Fii[7]*Wf.mat[2,4]/2 + Fii[4]*Wf.mat[1,4]*par.list$r*(1 - par.list$r)/2 + Fii[6]*Wf.mat[2,3]*par.list$r*(1 - par.list$r)/2)/(Wf.av(Fii=Fii, Wf.mat=Wf.mat)))
}

F33.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x3  <-  x3(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y3  <-  y3(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*x3*y3 + par.list$C*((Fii[8]*Wf.mat[3,3] + Fii[3]*Wf.mat[1,3]/4 + Fii[4]*Wf.mat[1,4]*(par.list$r^2)/4 + Fii[6]*Wf.mat[2,3]*((1 - par.list$r)^2)/4 + Fii[9]*Wf.mat[3,4]/4)/(Wf.av(Fii=Fii, Wf.mat=Wf.mat)))
}

F34.pr  <-  function(Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x3  <-  x3(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y4  <-  y4(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	x4  <-  x4(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y3  <-  y3(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*(x3*y4 + x4*y3) + par.list$C*((Fii[4]*Wf.mat[1,4]*par.list$r*(1 - par.list$r)/2 + Fii[6]*Wf.mat[2,3]*par.list$r*(1 - par.list$r)/2 + Fii[9]*Wf.mat[3,4]/2)/(Wf.av(Fii=Fii, Wf.mat=Wf.mat)))
}

F44.pr  <-  function(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...) {
	x4  <-  x4(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	y4  <-  y4(Fii = Fii, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list,...)
	(1 - par.list$C)*x4*y4 + par.list$C*((Fii[10]*Wf.mat[4,4] + Fii[4]*Wf.mat[1,4]*((1 - par.list$r)^2)/4 + Fii[6]*Wf.mat[2,3]*(par.list$r^2)/4 + Fii[7]*Wf.mat[2,4]/4 + Fii[9]*Wf.mat[3,4]/4)/(Wf.av(Fii=Fii, Wf.mat=Wf.mat)))
}



###########################################
##  Eigenvalues from Mathematica solutions
##  using quasi-equibirium expressions for
##  genotypic frequencies.

##  Additive fitness effects

lambda.AB1.add  <-  function(par.list) {
	C   <-  par.list$C
	sm  <-  par.list$sm
	sf  <-  par.list$sf

	(4 + sf*(-3 + sm) - sm + C*(-2 - sf*(-3 + sm) + sm)) / (2*(-2 + C)*(-1 + sf))
}

lambda.AB2.add  <-  function(par.list) {
	C   <-  par.list$C
	sm  <-  par.list$sm
	sf  <-  par.list$sf
	r   <-  par.list$r

(-8 + 4*C + 8*r - 8*C^2*r + 12*sf - 12*C*sf - 12*r*sf + 8*C*r*sf + 
 4*C^2*r*sf - 5*sf^2 + 6*C*sf^2 - C^2*sf^2 + 5*r*sf^2 - 6*C*r*sf^2 + 
 C^2*r*sf^2 - 4*(-1 + C)*(1 + (-1 + C)*r)*(-1 + sf)^2*sm + (-1 + C)*
 (1 + C + (-1 + C)*r)*(-1 + sf)^2*sm^2) / (4*(-2 + C)*(-1 + sf)^2)
}

lambda.ab1.add  <-  function(par.list) {
	C   <-  par.list$C
	sm  <-  par.list$sm
	sf  <-  par.list$sf

(4 + sf*(-1 + sm) - 3*sm + C*(-2 + sf*(-1 + sm) + sm)) / (2*(-2 + C)*(-1 + sm))
}

lambda.ab2.add  <-  function(par.list) {
	C   <-  par.list$C
	sm  <-  par.list$sm
	sf  <-  par.list$sf
	r   <-  par.list$r

1 / (4*(-2 + C)*(-1 + sm)^2)*((-1 + r)*(8 - 4*sf*(-1 + sm)^2 + sf^2*(-1 + sm)^2 + 
      sm*(-12 + 5*sm)) + 2*C*(2 - 2*(-1 + 2*r)*sf*(-1 + sm)^2 + (-1 + r)*sf^2*(-1 + sm)^2 + 
      sm*(-2 + r*(-4 + 3*sm))) - C^2*(sf^2*(-1 + sm)^2 - sm^2 + r*(8 - 12*sf*(-1 + sm)^2 + 
      3*sf^2*(-1 + sm)^2 + sm*(-20 + 11*sm))))
}





##  Dominance Reversal fitness effects

lambda.AB1.domRev  <-  function(par.list) {
	C   <-  par.list$C
	sm  <-  par.list$sm
	sf  <-  par.list$sf

	(8 + C*(-4 + 6*sf) + sf*(-5 + sm) - sm - C^2*(sf + (-1 + sf)*sm)) / (4*(-2 + C)*(-1 + sf))
}

lambda.AB2.domRev  <-  function(par.list) {
	C   <-  par.list$C
	sm  <-  par.list$sm
	sf  <-  par.list$sf
	r   <-  par.list$r

	(1 / (16*(-2 + C)*(-1 + sf)^2))*
	((-1 + r)*(32 - 40*sf + 17*sf^2 - 
    8*(-1 + sf)^2*sm + (-1 + sf)^2*sm^2) + 
	2*C*(8 - 3*sf*(8 - 4*sf + r*(-8 + 5*sf)) + 
	8*r*(-1 + sf)^2*sm - (3 + r)*(-1 + sf)^2*sm^2) + 
	C^2*(-32*r + 8*sf - 8*r*sf - 7*sf^2 + 13*r*sf^2 - 
	8*(1 + r)*(-1 + sf)^2*sm + (7 + r)*(-1 + sf)^2*sm^2))
}

lambda.ab1.domRev  <-  function(par.list) {
	C   <-  par.list$C
	sm  <-  par.list$sm
	sf  <-  par.list$sf

	-((-8 + sf + C*(4 + (2 + C)*sf) + 5*sm - (C^2 + (1 + C)^2*sf)*sm)/(4*(-2 + C)*(-1 + sm)))
}

lambda.ab2.domRev  <-  function(par.list) {
	C   <-  par.list$C
	sm  <-  par.list$sm
	sf  <-  par.list$sf
	r   <-  par.list$r

	(1/(16*(-2 + C)*(-1 + sm)^2))*
	(2*C*(8 - 8*(-1 + r)*sf*(-1 + sm)^2 + (-4 + r)*sf^2*(-1 + sm)^2 - 
	3*sm^2 + 3*r*sm*(-8 + 5*sm)) + (-1 + r)*(32 - 8*sf*(-1 + sm)^2 +
	sf^2*(-1 + sm)^2 + sm*(-40 + 17*sm)) + 
	C^2*(8*sf*(-1 + sm)^2 - 7*sf^2*(-1 + sm)^2 + sm*(-8 + 7*sm) + 
	r*(-32 + 24*sf*(-1 + sm)^2 - 3*sf^2*(-1 + sm)^2 + (88 - 47*sm)*sm)))
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

#############################
##  Additive fitness effects

########
##  Invasion based on lambda.AB1 for female beneficial 
##  allele under obligate outcrossing
inv.lAB1.add.obOut  <-  function(hf, hm, sm) {
	sm/(1 + sm)
}

##  Invasion based on lambda.ab1 for female beneficial 
##  allele under obligate outcrossing
inv.lab1.add.obOut  <-  function(hf, hm, sm) {
	sm/(1 - sm)
}


##  Invasion based on lambda.AB2 (w/ recombination)
##  for female beneficial allele under obligate outcrossing
inv.lAB2.add.obOut  <-  function(hf, hm, sm, r) {

	(2 - (-4 + sm)*sm - sqrt(-(-1 + r)*(4 +r*(-2 + sm)^2 - (-4 + sm)*sm)) + 
		r*(6 + (-4 + sm)*sm)) / (3 - (-4 + sm)*sm + r*(5 + (-4 + sm)*sm))
}

##  Invasion based on lambda.ab2 (w/ recombination)
##  for male beneficial allele under obligate outcrossing
inv.lab2.add.obOut  <-  function(hf, hm, sm, r) {

	2 + (sqrt(-(((-1 + r)*(4 + r*(-2 + sm)^2 + sm*(-12 + 7*sm))) / (-1 + sm)^2)) / (-1 + r))
}


########
##  Invasion based on lambda.AB1 for female beneficial 
##  allele (w/ selfing)
inv.lAB1.add  <-  function(hf, hm, sm, C) {
	(sm - C*sm)/(1 + C + sm - C*sm)
}


##  Invasion based on lambda.ab1 for female beneficial 
##  allele (w/ selfing)
inv.lab1.add  <-  function(hf, hm, sm, C) {
	((-1 + C)*sm) / ((1 + C)*(-1 + sm))
}

##  Invasion based on lambda.AB2 (w/ recombination)
##  for female beneficial allele (w/ selfing)
inv.lAB2.add  <-  function(hf, hm, sm, r, C) {

	(2 + 2*C + 6*r - 2*C*(2 + C)*r + 4*sm - 
   4*(C + (-1 + C)^2*r)*sm + (-1 + C)*(1 + C + (-1 + C)*r)*sm^2 - 
   sqrt(4*((1 + C)^2 - 
       2*(-1 + C)^2*C*(1 + C)*r + (-1 + C)^3*(1 + 3*C)*r^2) - 
    4*(-1 + C)*(1 + (-1 + C)*r)*((1 + 
         C)^2 + (-1 + C)*(1 + 3*C)*r)*sm + (-1 + C)*(1 + 
       C + (-1 + C)*r)*((1 + C)^2 + (-1 + C)*(1 + 3*C)*r)*sm^2)) / 
	(3 + C*(2 + C*(-1 + r) - 6*r) + 5*r + 4*sm - 
   4*(C + (-1 + C)^2*r)*sm + (-1 + C)*(1 + C + (-1 + C)*r)*sm^2)
}

##  Invasion based on lambda.ab2 (w/ recombination)
##  for male beneficial allele (w/ selfing)
inv.lab2.add  <-  function(hf, hm, sm, r, C) {

	1 / (1 - 2*C*(-1 + r) - r + 
  C^2*(1 + 3*r))*(2 + 2*C - 2*r - 4*C*r + 6*C^2*r - 
   2*sqrt(1 / ((-2 + C)^2*(-1 + sm)^2)*(4 - r^2*(-2 + sm)^2 - 
         12*sm + 2*r*(4 - 3*sm)*sm + 7*sm^2 - 
         2*C*(-1 + r)*(4 - 10*sm + 5*sm^2) + 
         C^4*(3*r^2*(-2 + sm)^2 + sm^2 - 4*r*(2 - 5*sm + 2*sm^2)) - 
         2*C^3*(4*r^2*(-2 + sm)^2 + (-2 + sm)*sm + 
            r*(-4 + 2*sm + 3*sm^2)) + 
         C^2*(4 + 6*r^2*(-2 + sm)^2 - 4*sm + 
            r*(8 - 44*sm + 30*sm^2)))) + 
   C*sqrt(1 / ((-2 + C)^2*(-1 + sm)^2)*(4 - r^2*(-2 + sm)^2 - 
         12*sm + 2*r*(4 - 3*sm)*sm + 7*sm^2 - 
         2*C*(-1 + r)*(4 - 10*sm + 5*sm^2) + 
         C^4*(3*r^2*(-2 + sm)^2 + sm^2 - 4*r*(2 - 5*sm + 2*sm^2)) - 
         2*C^3*(4*r^2*(-2 + sm)^2 + (-2 + sm)*sm + 
            r*(-4 + 2*sm + 3*sm^2)) + 
         C^2*(4 + 6*r^2*(-2 + sm)^2 - 4*sm + 
            r*(8 - 44*sm + 30*sm^2)))))
}




#######################################
##  Dominance Reversal fitness effects

########
##  Invasion based on lambda.AB1 for female beneficial 
##  allele under obligate outcrossing
inv.lAB1.domRev.obOut  <-  function(hf, hm, sm) {
	sm/(3 + sm)
}

##  Invasion based on lambda.ab1 for female beneficial 
##  allele under obligate outcrossing
inv.lab1.domRev.obOut  <-  function(hf, hm, sm) {
	-((3*sm)/(-1 + sm))
}

##  Invasion based on lambda.AB2 (w/ recombination)
##  for female beneficial allele under obligate outcrossing
inv.lAB2.domRev.obOut  <-  function(hf, hm, sm, r) {

	(12 + 8*sm - sm^2 - 3*sqrt(-(-1 + r)*(16 + r*(-4 + sm)^2 - (-8 + sm)*sm)) + r*(20 + (-8 + sm)*sm)) / 
		(15 - (-8 + sm)*sm + r*(17 + (-8 + sm)*sm))
}

##  Invasion based on lambda.ab2 (w/ recombination)
##  for male beneficial allele under obligate outcrossing
inv.lab2.domRev.obOut  <-  function(hf, hm, sm, r) {

	4 + sqrt(-(((-1 + r)*(16 + r*(-4 + sm)^2 + sm*(-56 + 31*sm))) / (-1 + sm)^2)) / (-1 + r)
}

########
##  Invasion based on lambda.AB1 for female beneficial 
##  allele (w/ selfing)
inv.lAB1.domRev  <-  function(hf, hm, sm, C) {

	((-1 + C)*sm) / (-3 + C + (-1 + C)*sm)
}


##  Invasion based on lambda.ab1 for female beneficial 
##  allele (w/ selfing)
inv.lab1.domRev  <-  function(hf, hm, sm, C) {

	-(((-3 + C)*(-1 + C)*sm) / ((1 + C)^2*(-1 + sm)))
}

##  Invasion based on lambda.AB2 (w/ recombination)
##  for female beneficial allele (w/ selfing)
inv.lAB2.domRev  <-  function(hf, hm, sm, r, C) {
	
	(12 + 20*r + 8*sm - 8*r*sm - sm^2 + r*sm^2 + 
		C^2*(-4 - 8*sm + 7*sm^2 + r*(4 - 8*sm + sm^2)) - 
		2*C*(-4 + 3*sm^2 + 
		r*(12 - 8*sm + sm^2)) - sqrt(-9*(-16 + r^2*(-4 + sm)^2 - 
			8*sm - 2*r*(-8 + sm)*sm + sm^2) + 
		2*C^2*(27*r^2*(-4 + sm)^2 + 8*(-2 - 5*sm + sm^2) + 
			r*(128 + 176*sm + 5*sm^2)) + 
		2*C*(96 + 32*sm - 31*sm^2 + r*(-256 - 32*sm + 31*sm^2)) + 
		C^4*(16 + 27*r^2*(-4 + sm)^2 + 8*sm - 7*sm^2 + 
			4*r*(-64 - 52*sm + 47*sm^2)) - 
		2*C^3*(32 + 36*r^2*(-4 + sm)^2 + 32*sm - 31*sm^2 + 
			r*(-256 - 32*sm + 139*sm^2)))) / 
	(15 + 8*sm - sm^2 + r*(17 - 8*sm + sm^2) + 
		C^2*(-7 - 8*sm + 7*sm^2 + r*(13 - 8*sm + sm^2)) - 
		2*C*(-4 + 3*sm^2 + r*(15 - 8*sm + sm^2)))
}

##  Invasion based on lambda.ab2 (w/ recombination)
##  for male beneficial allele (w/ selfing)

inv.lab2.domRev  <-  function(hf, hm, sm, r, C) {

	1 / (1 - 2*C*(-4 + r) - r + 
		C^2*(7 + 3*r))*(4 + 8*C + 4*C^2 - 4*r - 8*C*r + 12*C^2*r - 
		2*sqrt(1 / ((-2 + C)^2*(-1 + sm)^2)*(16 - r^2*(-4 + sm)^2 - 
			56*sm + 31*sm^2 - 6*r*sm*(-8 + 5*sm) + 
			2*C*(32 - 144*sm + 81*sm^2 + r*(64 - 48*sm + 15*sm^2)) + 
			C^4*(16 + 3*r^2*(-4 + sm)^2 - 88*sm + 65*sm^2 - 
				4*r*(32 - 100*sm + 53*sm^2)) - 
			2*C^3*(-32 + 4*r^2*(-4 + sm)^2 - 16*sm + 17*sm^2 + 
				r*(64 - 112*sm + 59*sm^2)) + 
			2*C^2*(3*r^2*(-4 + sm)^2 + 8*(6 - 7*sm + 2*sm^2) + 
				r*(64 - 288*sm + 165*sm^2)))) + 
		C*sqrt(1 / ((-2 + C)^2*(-1 + sm)^2)*(16 - r^2*(-4 + sm)^2 - 
			56*sm + 31*sm^2 - 6*r*sm*(-8 + 5*sm) + 
			2*C*(32 - 144*sm + 81*sm^2 + r*(64 - 48*sm + 15*sm^2)) + 
			C^4*(16 + 3*r^2*(-4 + sm)^2 - 88*sm + 65*sm^2 - 
				4*r*(32 - 100*sm + 53*sm^2)) - 
			2*C^3*(-32 + 4*r^2*(-4 + sm)^2 - 16*sm + 17*sm^2 + 
				r*(64 - 112*sm + 59*sm^2)) + 
			2*C^2*(3*r^2*(-4 + sm)^2 + 8*(6 - 7*sm + 2*sm^2) + 
				r*(64 - 288*sm + 165*sm^2)))))
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
#'				   r    =  0.5
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
	if(any(par.list[2:7] < 0) | any(par.list[2:7] > 1) | par.list[7] > 0.5)
		stop('The chosen parameter values fall outside of the reasonable bounds')

	if(par.list$hf  !=  par.list$hm)
		stop('please set male and female dominance values to be equal, and either 0.5 or 0.25')

	if(par.list$hf != 0.5 & par.list$hf != 0.25)
		stop('please set male and female dominance values to be equal, and either 0.5 or 0.25')

	##  Calculate Eigenvalues from analytic solutions 
	##  using quasi-equibirium genotypic frequencies
	if (par.list$hf == 0.5) {
		l.AB1  <- lambda.AB1.add(par.list)
		l.AB2  <- lambda.AB2.add(par.list)
		l.ab1  <- lambda.ab1.add(par.list)
		l.ab2  <- lambda.ab2.add(par.list)
	}
	if (par.list$hf == 0.25) {
		l.AB1  <- lambda.AB1.domRev(par.list)
		l.AB2  <- lambda.AB2.domRev(par.list)
		l.ab1  <- lambda.ab1.domRev(par.list)
		l.ab2  <- lambda.ab2.domRev(par.list)
	}

	
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
#'				   r    =  0.5
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
recursionFwdSim  <-  function(par.list, threshold = 1e-6) {

	##  Warnings
	if(any(par.list[2:7] < 0) | any(par.list[2:7] > 1) | par.list[7] > 0.5)
		stop('The chosen parameter values fall outside of the reasonable bounds')

	if(par.list$hf  !=  par.list$hm)
		stop('please set male and female dominance values to be equal, and either 0.5 or 0.25')

	if(par.list$hf != 0.5 & par.list$hf != 0.25)
		stop('please set male and female dominance values to be equal, and either 0.5 or 0.25')



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

	##  Initial frequencies
	if(par.list$sf > par.list$sm)
		Fii.init    <-  c(0.01,0,0,0,0,0,0,0,0,0.99)
	if(par.list$sf < par.list$sm)
		Fii.init    <-  c(0.99,0,0,0,0,0,0,0,0,0.01)


	##  Generation Loop
		# initialize
		Fii.gen[1,1]   <-  round(F11.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		Fii.gen[1,2]   <-  round(F12.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		Fii.gen[1,3]   <-  round(F13.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		Fii.gen[1,4]   <-  round(F14.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		Fii.gen[1,5]   <-  round(F22.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		Fii.gen[1,6]   <-  round(F23.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		Fii.gen[1,7]   <-  round(F24.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		Fii.gen[1,8]   <-  round(F33.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		Fii.gen[1,9]   <-  round(F34.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		Fii.gen[1,10]  <-  round(F44.pr(Fii = Fii.init, Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)


	# Start simulation
	i      <-  2
	diffs  <-  rep(1,10)

	while (i < par.list$gen & any(diffs[diffs != 0] > threshold)) {
		Fii.gen[i,1]   <-  round(F11.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		Fii.gen[i,2]   <-  round(F12.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		Fii.gen[i,3]   <-  round(F13.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		Fii.gen[i,4]   <-  round(F14.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		Fii.gen[i,5]   <-  round(F22.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		Fii.gen[i,6]   <-  round(F23.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		Fii.gen[i,7]   <-  round(F24.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		Fii.gen[i,8]   <-  round(F33.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		Fii.gen[i,9]   <-  round(F34.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		Fii.gen[i,10]  <-  round(F44.pr(Fii = Fii.gen[i-1,], Wf.mat = Wf.mat, Wm.mat = Wm.mat, par.list = par.list), digits=6)
		
		diffs  <-  Fii.gen[i,] - Fii.gen[i-1,]
		i      <-  i+1
	}

	##  Is equilibrium polymorphic?
	if (any(Fii.gen[i-1,c(1,10)] > 0.999)) # & all(Fii.gen[i-1,2:9] < 1e-4))
		 Poly  <-  0
	else Poly  <-  1

	##  Calculate Eigenvalues from analytic solutions 
	##  using quasi-equibirium genotypic frequencies
	if (par.list$hf == 0.5) {
		l.AB1  <- lambda.AB1.add(par.list)
		l.AB2  <- lambda.AB2.add(par.list)
		l.ab1  <- lambda.ab1.add(par.list)
		l.ab2  <- lambda.ab2.add(par.list)
	}
	if (par.list$hf == 0.25) {
		l.AB1  <- lambda.AB1.domRev(par.list)
		l.AB2  <- lambda.AB2.domRev(par.list)
		l.ab1  <- lambda.ab1.domRev(par.list)
		l.ab2  <- lambda.ab2.domRev(par.list)
	}

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

recursionFwdSimLoop  <-  function(n = 10000, gen = 5000, sRange = c(0,1), C = 0, hf = 0.5, hm = 0.5, r.vals = c(0.5, 0.2, 0.1, 0), threshold = 1e-7) {

	## Warnings
	if(any(c(C,hf,hm,r.vals) < 0) | any(c(C,hf,hm) > 1) | any(r.vals > 0.5))
		stop('At least one of the chosen parameter values fall outside of the reasonable bounds')

	if(threshold > 1e-7)
		stop('Carefully consider whether you want to change this threshold, 
			  as it will effect whether the simulations agree with the analytic results')

	#  initialize selection coeficients and storage structures
	s.vals   <-  matrix(runif(2*n, min=sRange[1], max=sRange[2]), ncol=2)
	Poly     <-  c()
	eigPoly  <-  c()
	agree    <-  c()


	##  Simulation Loop over values of r, sm, sf for fixed selfing rate (C)
	for (i in 1:length(r.vals)) {
			for (j in 1:nrow(s.vals)) {
				
				par.list  <-  list(
								   gen  =  gen,
								   C    =  C,
								   sf   =  s.vals[j,1],
								   sm   =  s.vals[j,2],
								   hm   =  hm,
								   hf   =  hf,
								   r    =  r.vals[i]
								  )
				res      <-  recursionFwdSim(par.list = par.list, threshold = threshold)
				Poly[(i-1)*nrow(s.vals) + j]     <-  res$Poly
				eigPoly[(i-1)*nrow(s.vals) + j]  <-  res$eigPoly
				agree[(i-1)*nrow(s.vals) + j]    <-  res$agree
		}
	}

	#  Compile results as data.frame
	results.df  <-  data.frame("hf"      = rep(0.5, length(r.vals)*length(s.vals)),
							   "hm"      = rep(0.5, length(r.vals)*length(s.vals)),
							   "C"       = rep(0,   length(r.vals)*length(s.vals)),
							   "r"       = c(rep(r.vals[1],length(s.vals)), 
							   		  		 rep(r.vals[2],length(s.vals)),
							   		  		 rep(r.vals[3],length(s.vals)),
							   		  		 rep(r.vals[4],length(s.vals))),
							   "sf"      = s.vals[,1],
							   "sm"      = s.vals[,2],
							   "Poly"    = Poly,
							   "eigPoly" = eigPoly,
							   "agree"   = agree
							   )

	#  Write results.df to .txt file
	filename  <-  paste("./output/data/simResults/recFwdSimLoop.out", "_C", C, "_hf", hf, "_hm", hm, "_sMax",sRange[2], ".txt", sep="")
	write.table(results.df, file=filename, col.names = TRUE, row.names = FALSE)

	#  Return results.df in case user wants it
	return(results.df)
}










#' Quick and dirty invasion analysis, modified to work with apply()
#'
#' @titleProportion Quick and dirty invasion analysis, modified to work with apply()
#' @param x vector of sm values over which to perform Invasion Analysis
#' @param par.list List of parameters.  as in all other functions
#' @seealso `propPrPFast`, `eigenInvAnalysis`
#' @export
#' @author Colin Olito.
#' @examples
#' propPrPFast(n = 10000, C = 0, hf = 0.5, hm = 0.5)

fastInv  <-  function(x, par.list, ...) {

	if(par.list$hf  !=  par.list$hm)
		stop('please set male and female dominance values to be equal, and either 0.5 or 0.25')

	if(par.list$hf != 0.5 & par.list$hf != 0.25)
		stop('please set male and female dominance values to be equal, and either 0.5 or 0.25')


	par.list$sf  <-  x[1]
	par.list$sm  <-  x[2]

	##  Calculate Eigenvalues from analytic solutions 
	##  using quasi-equibirium genotypic frequencies
	if (par.list$hf == 0.5) {
		l.AB1  <- lambda.AB1.add(par.list)
		l.AB2  <- lambda.AB2.add(par.list)
		l.ab1  <- lambda.ab1.add(par.list)
		l.ab2  <- lambda.ab2.add(par.list)
	}
	if (par.list$hf == 0.25) {
		l.AB1  <- lambda.AB1.domRev(par.list)
		l.AB2  <- lambda.AB2.domRev(par.list)
		l.ab1  <- lambda.ab1.domRev(par.list)
		l.ab2  <- lambda.ab2.domRev(par.list)
	}

	PrP  <-  0

	# Protected polymorphism
	if (any(c(l.AB1, l.AB2) > 1) & any(c(l.ab1, l.ab2) > 1 )) {
		PrP  <-  1			
		}
	
	##  Is polymorphism facilitated by recombination? 
	if (any(c(l.AB1,l.ab1) < 1) & l.AB2 > 1 &  l.ab2 > 1)
		PrP  <-  2
	PrP
}




#' Proportion of parameter space resulting in protected polymorphism
#' determined by evaluating eigenvalues. - Modified to implement apply()
#' for innermost loop.  about 1/3 times faster than regular propPrP()
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
#' propPrPFast(n = 10000, C = 0, hf = 0.5, hm = 0.5)

propPrPFast  <-  function(n = 1000, C = 0, hf = 0.5, hm = 0.5, sRange = c(0,1), weakSel = FALSE) {

	## Warnings
	if(any(c(C,hf,hm) < 0) | any(c(C,hf,hm) > 1))
		stop('At least one of the chosen parameter values fall outside of the reasonable bounds')

	if(weakSel == TRUE & max(sRange) > 0.1)
		stop('There appears to be a disagreement between the chosen range of selection coefficients and the weak selection option. I recommend using sRange = c(0,0.1)')

	#  initialize selection coeficients and storage structures
	if(weakSel == TRUE)
		r.vals      <-  seq(0, 0.01, by=0.0005)
	else 
		r.vals      <-  seq(0, 0.5, by=0.01)
	PrP     <-  c()
	rPrP    <-  c()


	##  Loop over values of r and sm for fixed selfing rate (C)
	##  calculating proportion of parameter space where PrP is 
	##  predicted each time.

	s.vals    <-  matrix(runif(2*n, min=sRange[1], max=sRange[2]), ncol=2)
	for (i in 1:length(r.vals)) {
		poly  <-  rep(0, times=nrow(s.vals))
		par.list  <-  list(
						   gen  =  NA,
						   C    =  C,
						   sm   =  NA,
						   sf   =  NA,
						   hm   =  hm,
						   hf   =  hf,
						   r    =  r.vals[i]
						  )

		poly  <-  apply(s.vals, 1,  fastInv, par.list=par.list)

		#  Calculate proportion of parameter space resulting in PrP
		PrP[i]   <-  sum(poly == 1 | poly == 2)/length(poly)
		rPrP[i]  <-  sum(poly == 2)/length(poly)
		print(r.vals[i])
		rm(poly)
	}

	#  Compile results as data.frame
	results.df  <-  data.frame("hf"      =  rep(hf, length(r.vals)),
							   "hm"      =  rep(hm, length(r.vals)),
							   "C"       =  rep(C,  length(r.vals)),
							   "r"       =  r.vals,
							   "PrP"     =  PrP,
							   "rPrP"    =  rPrP,
							   "slPrP"   =  PrP - rPrP
							   )

	#  Write results.df to .txt file
	if(weakSel == TRUE)
		sel  <-  "_weak"
	else 
		sel  <-  ""
	filename  <-  paste("./output/data/Fig2Data/propPrp.out", sel, "_C", C, "_hf", hf, "_hm", hm, "_n", n, ".txt", sep="")
	write.table(results.df, file=filename, col.names = TRUE, row.names = FALSE)

	#  Return results.df in case user wants it
	return(results.df)
}






Bottom line they really are a hermaphrodite and only the pussy works facts and true now  give me my money back fbi can see all this 
give my money back
