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

toPdf(Fig.1wk(), figPath(name='Fig1wk.pdf'), width=7, height=7)
embed_fonts(figPath(name='Fig1wk.pdf'))

toPdf(Fig.2(), figPath(name='Fig2.pdf'), width=7, height=7)
embed_fonts(figPath(name='Fig2.pdf'))

toPdf(recSimFig_add(), figPath(name='recSimFig_add.pdf'), width=7, height=7)
embed_fonts(figPath(name='recSimFig_add.pdf'))

toPdf(recSimFig_domRev(), figPath(name='recSimFig_domRev.pdf'), width=7, height=7)
embed_fonts(figPath(name='recSimFig_domRev.pdf'))

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
