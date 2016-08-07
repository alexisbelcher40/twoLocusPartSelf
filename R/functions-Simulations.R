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

y2p  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	((2*Fii[5]*Wm.mat[2,2]) + (Fii[2]*Wm.mat[1,2]) + (Fii[6]*Wm.mat[2,3]) + (Fii[7]*Wm.mat[2,4]))/(2*Wm.av(Fii, Wm.mat)) + 
		par.list$rm*(((Fii[4]*Wm.mat[1,4]) - (Fii[6]*Wm.mat[2,3]))/(2*Wm.av(Fii, Wm.mat)))
}

y3p  <-  function(Fii, Wf.mat, Wm.mat, par.list,...) {
	 ((2*Fii[8]*Wm.mat[3,3]) + (Fii[9]*Wm.mat[3,4]) + (Fii[3]*Wm.mat[1,3]) + (Fii[6]*Wm.mat[2,3]))/(2*Wm.av(Fii, Wm.mat)) + 
 		par.list$rm*(((Fii[4]*Wm.mat[1,4]) - (Fii[6]*Wm.mat[2,3]))/(2*Wmav))
}

y4p  <-  function(Fii, Wf.mat, Wm.mat, par.list,...) {
	((2*Fii[10]*Wm.mat[4,4]) + (Fii[9]*Wm.mat[3,4]) + (Fii[4]*Wm.mat[1,4]) + (Fii[7]*Wm.mat[2,4]))/(2*Wm.av(Fii, Wm.mat)) - 
		par.list$rm*(((Fii[4]*Wm.mat[1,4]) - (Fii[6]*Wm.mat[2,3]))/(2*Wm.av(Fii, Wm.mat)))
}


#########################################
## Genotypic frequency recursions


F11.pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	(1 – par.list$C)*x1*y1 + par.list$C*(Fii[1] + Fii[2]/4 + Fii[3]/4 + Fii[4]*((1 – par.list$rf)^2)/4 + Fii[6]*(par.list$rf^2)/4)	
} 

F12.pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	(1 – par.list$C)*(x1*y2 + x2*y1) + par.list$C*(Fii[2]/2 + Fii[4]*par.list$rf*(1 – par.list$rf)/2 + Fii[6]*par.list$rf*(1 – par.list$rf)/2
}

F13.pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	(1 – par.list$C)*(x1*y3 + x3*y1) + par.list$C*(Fii[3]/2 + Fii[4]*par.list$rf*(1 – par.list$rf)/2 + Fii[6]*par.list$rf*(1 – par.list$rf)/2
}

F14.pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	(1 – par.list$C)*(x1*y4 + x4*y1) + par.list$C*(Fii[4]*((1 – par.list$rf)^2)/2 + Fii[6]*(par.list$rf^2)/2
}

F22.pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	(1 – par.list$C)*x2*y2 + par.list$C*(Fii[5] + Fii[2]/4 + Fii[4]*(par.list$rf^2)/4 + Fii[6]*((1 – par.list$rf)^2)/4 + Fii[7]/4)
}

F23.pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	(1 – par.list$C)*(x2*y3 + x3*y2) + par.list$C*(Fii[4]*(par.list$rf^2)/2 + Fii[6]*((1 – par.list$rf)^2)/2)
}

F24.pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	(1 – par.list$C)*(x2*y4 + x4*y2) + par.list$C*(Fii[7]/2 + Fii[4]*par.list$rf*(1 – par.list$rf)/2 + Fii[6]*par.list$rf*(1 – par.list$rf)/2)
}

F33.pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	(1 – par.list$C)*(x3*y4 + x4*y3) + par.list$C*(Fii[4]*par.list$rf*(1 – par.list$rf)/2 + Fii[6]*par.list$rf*(1 – par.list$rf)/2 + Fii[9]/2)
}

F34.pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	(1 – par.list$C)*x3*y3 + par.list$C(Fii[8] + Fii[3]/4 + Fii[4]*(par.list$rf^2)/4 + Fii[6]*(1 – par.list$rf)2/4 + Fii[9]/4)
}

F44.pr  <-  function(Fii, Wf.mat, Wm.mat, par.list,...){
	(1 – par.list$C)*x4*y4 + par.list$C(Fii[10] + Fii[4]*((1 – par.list$rf)^2)/4 + Fii[6]*(par.list$rf^2)/4 + Fii[7]/4 + Fii[9]/4)
}