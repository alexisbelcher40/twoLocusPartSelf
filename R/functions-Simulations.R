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
