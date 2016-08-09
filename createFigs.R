#####################################################
#  2-locus SA with partial selfing
#
#  Functions to generate Figures for: 
#    Article title goes here...
#  
#  Authors: Colin Olito, Crispin Jordan, Tim Connallon
#
#
#  NOTES:  
#          

rm(list=ls())

library(extrafont)
library(fontcm)
loadfonts()

#source('paths.R')
source('R/functions-analyses.R')
source('R/functions-figures.R')


###############
# PAPER FIGURES
###############

toPdf(Fig.1(), figPath(name='Fig1.pdf'), width=7, height=7)
embed_fonts(figPath(name='Fig1.pdf'))

#toPdf(fig2(), figPath(name='fig2.pdf'), width=7, height=7)
#embed_fonts(figPath(name='fig2.pdf'))

#toPdf(fig3(), figPath(name='fig3.pdf'), width=7, height=7)
#embed_fonts(figPath(name='fig3.pdf'))

#toPdf(fig4(), figPath(name='fig4.pdf'), width=7, height=7)
#embed_fonts(figPath(name='fig4.pdf'))

#toPdf(fig5(), figPath(name='fig5.pdf'), width=9, height=3)
#embed_fonts(figPath(name='fig5.pdf'))

##################
# APPENDIX FIGURES
##################
#toPdf(figA1(), figPath(name='figA1.pdf'), width=8, height=8)
#embed_fonts(figPath(name='figA1.pdf'))

#toPdf(figA2(), figPath(name='figA2.pdf'), width=9, height=3)
#embed_fonts(figPath(name='figA2.pdf'))

#toPdf(figA3(), figPath(name='figA3.pdf'), width=9, height=3.5)
#embed_fonts(figPath(name='figA3.pdf'))
