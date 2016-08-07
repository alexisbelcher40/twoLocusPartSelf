#####################################################
#  2-locus SA with partial selfing
#
#  Necessary functions for deterministic simulations
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
 		par.list$rf*(((Fii[4]*Wf.mat[1,4]) - (Fii[6]*Wf.mat[2,3]))/(2*Wfav))
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

y2p  <-  function(Fii, Wm.mat, Wm.mat, par.list,...){
	((2*Fii[5]*Wm.mat[2,2]) + (Fii[2]*Wm.mat[1,2]) + (Fii[6]*Wm.mat[2,3]) + (Fii[7]*Wm.mat[2,4]))/(2*Wm.av(Fii, Wm.mat)) + 
		par.list$rm*(((Fii[4]*Wm.mat[1,4]) - (Fii[6]*Wm.mat[2,3]))/(2*Wm.av(Fii, Wm.mat)))
}

y3p  <-  function(Fii, Wm.mat, Wm.mat, par.list,...) {
	 ((2*Fii[8]*Wm.mat[3,3]) + (Fii[9]*Wm.mat[3,4]) + (Fii[3]*Wm.mat[1,3]) + (Fii[6]*Wm.mat[2,3]))/(2*Wm.av(Fii, Wm.mat)) + 
 		par.list$rm*(((Fii[4]*Wm.mat[1,4]) - (Fii[6]*Wm.mat[2,3]))/(2*Wmav))
}

y4p  <-  function(Fii, Wm.mat, Wm.mat, par.list,...) {
	((2*Fii[10]*Wm.mat[4,4]) + (Fii[9]*Wm.mat[3,4]) + (Fii[4]*Wm.mat[1,4]) + (Fii[7]*Wm.mat[2,4]))/(2*Wm.av(Fii, Wm.mat)) - 
		par.list$rm*(((Fii[4]*Wm.mat[1,4]) - (Fii[6]*Wm.mat[2,3]))/(2*Wm.av(Fii, Wm.mat)))
}


#########################################
## Genotypic frequency recursions


Fii[1].pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	(1 – C)*x1*y1 + C*(Fii[1] + Fii[2]/4 + Fii[3]/4 + Fii[4]*((1 – r)^2)/4 + Fii[6]*(r^2)/4)	
} 

Fii[2].pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	(1 – C)*(x1*y2 + x2*y1) + C*(Fii[2]/2 + Fii[4]*r*(1 – r)/2 + Fii[6]*r*(1 – r)/2
}

Fii[3].pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	(1 – C)*(x1*y3 + x3*y1) + C*(Fii[3]/2 + Fii[4]*r*(1 – r)/2 + Fii[6]*r*(1 – r)/2
}

Fii[4].pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	(1 – C)*(x1*y4 + x4*y1) + C*(Fii[4]*((1 – r)^2)/2 + Fii[6]*(r^2)/2
}

Fii[5].pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	(1 – C)*x2*y2 + C*(Fii[5] + Fii[2]/4 + Fii[4]*(r^2)/4 + Fii[6]*((1 – r)^2)/4 + Fii[7]/4)
}

Fii[6].pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){(
	1 – C)*(x2*y3 + x3*y2) + C*(Fii[4]*(r^2)/2 + Fii[6]*((1 – r)^2)/2)
}

Fii[7].pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	(1 – C)*(x2*y4 + x4*y2) + C*(Fii[7]/2 + Fii[4]*r*(1 – r)/2 + Fii[6]*r*(1 – r)/2)
}

Fii[9].pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	(1 – C)*(x3*y4 + x4*y3) + C*(Fii[4]*r*(1 – r)/2 + Fii[6]*r*(1 – r)/2 + Fii[9]/2)
}

Fii[8].pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	(1 – C)*x3*y3 + C(Fii[8] + Fii[3]/4 + Fii[4]*(r^2)/4 + Fii[6]*(1 – r)2/4 + Fii[9]/4)
}

Fii[10].pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	(1 – C)*x4*y4 + C(Fii[10] + Fii[4]*((1 – r)^2)/4 + Fii[6]*(r^2)/4 + Fii[7]/4 + Fii[9]/4)
}