#####################################################
#  2-locus SA with partial selfing
#
#  R code to run deterministic simulations
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

rm(list=ls())
#####################
##  Dependencies
source('R/functions-Simulations.R')

####################################
#  Simulation parameter values

par.list  <-  list(
				   gen  =  5000,
				   C    =  0,
				   sm   =  0.7,
				   sf   =  0.7,
				   hm   =  0.5,
				   hf   =  0.5,
				   rm   =  0.5,
				   rf   =  0.5
				   )



######################################
#  Initial genotypic frequencies

Fii.init  <-  c(0.99,0,0,0,0,0,0,0,0,0.01)
Fii.init  <-  c(0.01,0,0,0,0,0,0,0,0,0.99)


res  <-  twoLocusSAPartSelf(par.list = par.list, Fii.init = Fii.init)
res$Poly
res$EQ.freq

test  <-  res$Fii.gen
plot(test[,1] ~ c(1:nrow(test)), lwd=4, type='l', ylim=c(0,1))
for (i in 2:10) {
	lines(test[,i] ~ c(1:nrow(test)), lwd=4, col=i)
}

eigenInvAnalysis(par.list)


sf.seq  <-  seq(0.1,1,len=10)
sm.seq  <-  seq(0.1,1,len=10)

testAgree  <-  matrix(0,nrow=length(sf.seq),ncol=length(sm.seq))

par.list  <-  list(
				   gen  =  5000,
				   C    =  0,
				   sm   =  0,
				   sf   =  0,
				   hm   =  0.5,
				   hf   =  0.5,
				   rm   =  0.5,
				   rf   =  0.5
				   )



for (i in 1:length(sf.seq)) {
	for (j in 1:length(sm.seq)) {

		par.list$sf  <-  sf.seq[i]
		par.list$sm  <-  sm.seq[j]
		testAgree[i,j]  <-  twoLocusSAPartSelf(par.list = par.list, Fii.init = Fii.init)$agree
	}
print(i)
}

testAgree


library(rasterVis)
library(grid)
library(gridExtra)
library(extrafont)
library(fontcm)
loadfonts()

Greys  <- rev(gray.colors(25, start=0.1, end=1.0, gamma=1))
my.col.at <- seq(0,1, length=25)
my.lab.at <- seq(0,1, by=0.2)
myColorkey <- list(at=my.col.at, ## where the colors change
                   labels=list(
                       at=my.lab.at ## where to print labels
                     ))
length(my.col.at)

r <- raster(testAgree[11:1,])

testplot  <-  levelplot(r, margin=FALSE, contour=FALSE,
			main=expression(paste(testAgree)), 
			xlab=expression(paste(italic(s[m]))), 
			ylab=expression(paste(italic(s[f]))))
testplot



invKidwell.lAB(0.2)
invKidwell.lab(0.2)

inv.lAB1.obOut(hf=0.5, hm=0.5, sm=0.2)
inv.lab1.obOut(hf=0.5, hm=0.5, sm=0.2)
inv.lAB2.obOut(hf = 0.5, hm = 0.5, sm = 0.5, rf = 0.1, rm = 0.1)	
inv.lab2.obOut(hf = 0.5, hm = 0.5, sm = 0.2, rf = 0.1, rm = 0.1)	




inv.lAB1(hf=0.5, hm=0.5, sm=0.2, C=0.5)
inv.lab1(hf=0.5, hm=0.5, sm=0.2, C=0.5)
inv.lAB2(hf = 0.5, hm = 0.5, sf = 0.5, sm = 0.5, rf = 0.1, rm = 0.1, C = 0.5)	
inv.lab2(hf = 0.5, hm = 0.5, sf = 0.5, sm = 0.5, rf = 0.1, rm = 0.1, C = 0.5)	
