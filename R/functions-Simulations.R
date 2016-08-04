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

Wf.av  <-  function(Fii.list, Wf.list, par.list, ...){
   (F11*wf11) + (F12*wf12) + (F13*wf13) + (F14*wf14) + (F22*wf22) + (F23*wf23) + (F24*wf24) + (F33*wf33) + (F34*wf34) + (F44*wf44)
}
Wm.av  <-  function(Fii.list, Wm.list, par.list, ...){
   (F11*wm11) + (F12*wm12) + (F13*wm13) + (F14*wm14) + (F22*wm22) + (F23*wm23) + (F24*wm24) + (F33*wm33) + (F34*wm34) + (F44*wm44)
}


#########################################
## Haplotype frequency change in gametes


x1p  <-  function(Fii, Wf, Wm, ...) {

}


x1p[wf11_, wf12_, wf13_, wf14_, wf22_, wf23_, wf24_, wf33_, wf34_, wf44_][rf_][C_] := 

((2*\[Phi]11[q1, C]*wf11) + (\[Phi]12[q1, q2, C]*
      wf12) + (\[Phi]13[q1, q3, C]*wf13) + (\[Phi]14[q1, q4, C]*
      wf14)) / (2*Wfav[q1, q2, q3, q4][wf11, wf12, wf13, wf14, wf22, wf23, wf24, wf33, wf34, wf44][C]) - rf*(((\[Phi]14[q1, q4, C]*wf14) - (\[Phi]23[q2, q3, C]*wf23))/(
     2*Wfav[q1, q2, q3, q4][wf11, wf12, wf13, wf14, wf22, wf23, wf24, 
        wf33, wf34, wf44][C]));