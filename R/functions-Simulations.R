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

Wf.av  <-  function(Fii.list, Wf.list, ...){
   (F11*wf11) + (F12*wf12) + (F13*wf13) + (F14*wf14) + (F22*wf22) + (F23*wf23) + (F24*wf24) + (F33*wf33) + (F34*wf34) + (F44*wf44)
}
Wm.av  <-  function(Fii.list, Wm.list, ...){
   (F11*wm11) + (F12*wm12) + (F13*wm13) + (F14*wm14) + (F22*wm22) + (F23*wm23) + (F24*wm24) + (F33*wm33) + (F34*wm34) + (F44*wm44)
}


#########################################
## Haplotype frequency change in gametes

#  Ovules
x1p  <-  function(Fii.list, Wf.list, Wm.list, par.list,...) {
	((2*F11*wf11) + (F12*wf12) + (F13*wf13) + (F14*wf14)) / (2*Wf.av(Fii.list, Wf.list)) - 
		rf*(((F14*wf14) - (F23*wf23))/(2*Wf.av(Fii.list, Wf.list)))
}

x2p  <-  function(Fii.list, Wf.list, Wm.list, par.list,...){
	((2*F22*wf22) + (F12*wf12) + (F23*wf23) + (F24*wf24))/(2*Wf.av(Fii.list, Wf.list)) + 
		rf*(((F14*wf14) - (F23*wf23))/(2*Wf.av()))
}

x3p  <-  function(Fii.list, Wf.list, Wm.list, par.list,...) {
	 ((2*F33*wf33) + (F34*wf34) + (F13*wf13) + (F23*wf23))/(2*Wf.av()) + 
 		rf*(((F14*wf14) - (F23*wf23))/(2*Wfav))
}

x4p  <-  function(Fii.list, Wf.list, Wm.list, par.list,...) {
	((2*F44*wf44) + (F34*wf34) + (F14*wf14) + (F24*wf24))/(2*Wf.av()) - 
		rf*(((F14*wf14) - (F23*wf23))/(2*Wf.av()))
}

#  Pollen/Sperm
y1p  <-  function(Fii.list, Wf.list, Wm.list, par.list,...) {
	((2*F11*wm11) + (F12*wm12) + (F13*wm13) + (F14*wm14)) / (2*Wm.av(Fii.list, Wm.list)) - 
		rf*(((F14*wm14) - (F23*wm23))/(2*Wm.av(Fii.list, Wm.list)))
}

y2p  <-  function(Fii.list, Wm.list, Wm.list, par.list,...){
	((2*F22*wm22) + (F12*wm12) + (F23*wm23) + (F24*wm24))/(2*Wm.av(Fii.list, Wm.list)) + 
		rf*(((F14*wm14) - (F23*wm23))/(2*Wm.av()))
}

y3p  <-  function(Fii.list, Wm.list, Wm.list, par.list,...) {
	 ((2*F33*wm33) + (F34*wm34) + (F13*wm13) + (F23*wm23))/(2*Wm.av()) + 
 		rf*(((F14*wm14) - (F23*wm23))/(2*Wmav))
}

y4p  <-  function(Fii.list, Wm.list, Wm.list, par.list,...) {
	((2*F44*wm44) + (F34*wm34) + (F14*wm14) + (F24*wm24))/(2*Wm.av()) - 
		rf*(((F14*wm14) - (F23*wm23))/(2*Wm.av()))
}


#########################################
## Genotypic frequency recursions


F11.pr = (1 – C)*x1*y1 + C*(F11 + F12/4 + F13/4 + F14*((1 – r)^2)/4 + F23*(r^2)/4)
F12.pr = (1 – C)*(x1*y2 + x2*y1) + C*(F12/2 + F14*r*(1 – r)/2 + F23*r*(1 – r)/2
F13.pr = (1 – C)*(x1*y3 + x3*y1) + C*(F13/2 + F14*r*(1 – r)/2 + F23*r*(1 – r)/2
F14.pr = (1 – C)*(x1*y4 + x4*y1) + C*(F14*((1 – r)^2)/2 + F23*(r^2)/2
F22.pr = (1 – C)*x2*y2 + C*(F22 + F12/4 + F14*(r^2)/4 + F23*((1 – r)^2)/4 + F24/4)
F23.pr = (1 – C)*(x2*y3 + x3*y2) + C*(F14*(r^2)/2 + F23*((1 – r)^2)/2)
F24.pr = (1 – C)*(x2*y4 + x4*y2) + C*(F24/2 + F14*r*(1 – r)/2 + F23*r*(1 – r)/2)
F34.pr = (1 – C)*(x3*y4 + x4*y3) + C*(F14*r*(1 – r)/2 + F23*r*(1 – r)/2 + F34/2)
F33.pr = (1 – C)*x3*y3 + C(F33 + F13/4 + F14*(r^2)/4 + F23*(1 – r)2/4 + F34/4)
F44.pr = (1 – C)*x4*y4 + C(F44 + F14*((1 – r)^2)/4 + F23*(r^2)/4 + F24/4 + F34/4)