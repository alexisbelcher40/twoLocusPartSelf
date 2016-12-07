#####################################################
#  2-locus SA with partial selfing
#
#  Functions to generate Figures for: 
#    Article title goes here...
#  
#  Author: Colin Olito
#
#
#  NOTES: Run this file, either from terminal using Rscript,
#		  or interactively in R. This should create all the 
#		  figures needed to correctly compile the mansucript
#		  LaTeX file.  
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

toPdf(Fig.2wk(), figPath(name='Fig2wk.pdf'), width=7, height=7)
embed_fonts(figPath(name='Fig2wk.pdf'))

toPdf(recSimFig_add(), figPath(name='recSimFig_add.pdf'), width=7, height=7)
embed_fonts(figPath(name='recSimFig_add.pdf'))

toPdf(recSimFig_domRev(), figPath(name='recSimFig_domRev.pdf'), width=7, height=7)
embed_fonts(figPath(name='recSimFig_domRev.pdf'))

