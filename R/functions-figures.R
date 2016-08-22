###############
# DEPENDENCIES
###############
# library(...)


#######################
# AUXILLIARY FUNCTIONS
#######################

toPdf <- function(expr, filename, ...) {
  toDev(expr, pdf, filename, ...)
}

figPath  <-  function(name) {
  file.path('output/figures', name)
}

toDev <- function(expr, dev, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf('Creating %s\n', filename))
  dev(filename, family='CM Roman', ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}



####################
# PLOTTING FUNCTIONS
####################

#' Plot text or points according to relative axis position.
#'
#' @title Plot text or points according to relative axis position
#' @param px Relative x-axis position (in proportion) where character is to be plotted.
#' @param py Relative y-axis position (in proportion) where character is to be plotted.
#' @param lab Plotted text. Works if argument \code{\link[graphics]{text}} is \code{TRUE}.
#' @param adj See argument of same name in R base function \code{\link[graphics]{par}}.
#' @param text Logical. Should text or points be plotted?
#' @param log Used if the original plot uses the argument log, e.g. \code{log='x'}, \code{log='y'} or \code{log='xy'}.
#' @param ... Additional arguments to R base function \code{\link[graphics]{text}}.
#' @export
proportionalLabel <- function(px, py, lab, adj=c(0, 1), text=TRUE, log=FALSE, ...) {
    usr  <-  par('usr')
    x.p  <-  usr[1] + px*(usr[2] - usr[1])
    y.p  <-  usr[3] + py*(usr[4] - usr[3])
    if(log=='x') {
        x.p<-10^(x.p)
    }
    if(log=='y') {
        y.p<-10^(y.p)
    }
    if(log=='xy') {
        x.p<-10^(x.p)
        y.p<-10^(y.p)
    }
    if(text){
        text(x.p, y.p, lab, adj=adj, ...)
    } else {
        points(x.p, y.p, ...)
    }
}

#' Draw equally-spaced white lines on plot window.
#'
#' @title Equally-spaced white lines on plot window
#' @param ... Additional arguments to internal function \code{\link{proportionalLabel}}.
#' @author Diego Barneche
#' @export
plotGrid  <-  function(lineCol='white',...) {
    proportionalLabel(rep(0.2, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.4, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.6, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.8, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.2, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.4, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.6, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.8, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
}


#' Internal. Create nice rounded numbers for plotting.
#'
#' @title Rounded numbers for plotting
#' @param value A numeric vector.
#' @param precision Number of rounding digits.
#' @return A character vector.
#' @author Diego Barneche.
rounded  <-  function(value, precision=1) {
  sprintf(paste0('%.', precision, 'f'), round(value, precision))
}


#' Creates transparent colours
#'
#' @title Creates transparent colours
#' @param col Colour.
#' @param opacity equivalent to alpha transparency parameter
#' @export
transparentColor <- function(col, opacity=0.5) {
    if (length(opacity) > 1 && any(is.na(opacity))) {
        n        <-  max(length(col), length(opacity))
        opacity  <-  rep(opacity, length.out=n)
        col      <-  rep(col, length.out=n)
        ok       <-  !is.na(opacity)
        ret      <-  rep(NA, length(col))
        ret[ok]  <-  Recall(col[ok], opacity[ok])
        ret
    } else {
        tmp  <-  col2rgb(col)/255
        rgb(tmp[1,], tmp[2,], tmp[3,], alpha=opacity)
    }
}






##############################################################
##############################################################







#' Fig.1: Analytic results Kidwell Plots
#' 
#'
#' @title Create a 2x2 panel plot with our analytic results showing Kidwell plots
#' 
#' @export

Fig.1  <-  function() {

# Color scheme
    colfunc <- colorRampPalette(c("#252525", "grey70"))
    COLS  <-  colfunc(6)
#    COLS  <-  c("black", "#525252", "#737373", "#bdbdbd")

#  Create vector of male selection coefficiets for invasion functions
sm  <-  seq(0,1,by=0.0001)

# Set plot layout
layout.mat <- matrix(c(1,2,3,4,5,6), nrow=2, ncol=3, byrow=TRUE)
layout <- layout(layout.mat,respect=TRUE)

##  Row 1: Additive allelic effects
    ##  Panel One: C = 0
        
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi.obOut   <-  inv.lab1.add.obOut(hf=0.5, hm=0.5, sm=sm)
        twoLoc.Hi.obOut[twoLoc.Hi.obOut > 1]  <-  1.00000001
        twoLoc.Hi.obOut[10001]  <-  1.00000001
        twoLoc.Lo.obOut  <-  inv.lAB1.add.obOut(hf=0.5, hm=0.5, sm=sm)
        
        r0.5.Hi  <-  inv.lab2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0.5)
        r0.5.Hi[r0.5.Hi > 1]  <-  1.00000001
        r0.5.Hi[r0.5.Hi == 'NaN']  <-  1.00000001
        r0.5.Lo  <-  inv.lAB2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0.5)

        r0.4.Hi  <-  inv.lab2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0.4)
        r0.4.Hi[r0.4.Hi > 1]  <-  1.00000001
        r0.4.Hi[r0.4.Hi == 'NaN']  <-  1.00000001
        r0.4.Lo  <-  inv.lAB2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0.4)

        r0.3.Hi  <-  inv.lab2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0.3)
        r0.3.Hi[r0.3.Hi > 1]  <-  1.00000001
        r0.3.Hi[r0.3.Hi == 'NaN']  <-  1.00000001
        r0.3.Lo  <-  inv.lAB2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0.3)

        r0.2.Hi  <-  inv.lab2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0.2)
        r0.2.Hi[r0.2.Hi > 1]  <-  1.00000001
        r0.2.Hi[r0.2.Hi == 'NaN']  <-  1.00000001
        r0.2.Lo  <-  inv.lAB2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0.2)

        r0.1.Hi  <-  inv.lab2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0.1)
        r0.1.Hi[r0.1.Hi > 1]  <-  1.00000001
        r0.1.Hi[r0.1.Hi == 'NaN']  <-  1.00000001
        r0.1.Lo  <-  inv.lAB2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0.1)

        r0.Hi  <-  inv.lab2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0)
        r0.Hi[r0.Hi > 1]  <-  1.00000001
        r0.Hi[r0.Hi == 'NaN']  <-  1.00000001
        r0.Lo  <-  inv.lAB2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0)

        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        #  w/ recombination
        lines(r0.5.Hi[r0.5.Hi > twoLoc.Hi.obOut & r0.5.Hi <= 1] ~ sm[r0.5.Hi > twoLoc.Hi.obOut & r0.5.Hi <= 1], lwd=2, col=COLS[1], lty=1)
        lines(r0.5.Lo[r0.5.Lo < twoLoc.Lo.obOut] ~ sm[r0.5.Lo < twoLoc.Lo.obOut], lwd=2, col=COLS[1], lty=1)
        lines(r0.4.Hi[r0.4.Hi > twoLoc.Hi.obOut & r0.4.Hi <= 1] ~ sm[r0.4.Hi > twoLoc.Hi.obOut & r0.4.Hi <= 1], lwd=2, col=COLS[2], lty=1)
        lines(r0.4.Lo[r0.4.Lo < twoLoc.Lo.obOut] ~ sm[r0.4.Lo < twoLoc.Lo.obOut], lwd=2, col=COLS[2], lty=1)
        lines(r0.3.Hi[r0.3.Hi > twoLoc.Hi.obOut & r0.3.Hi <= 1] ~ sm[r0.3.Hi > twoLoc.Hi.obOut & r0.3.Hi <= 1], lwd=2, col=COLS[3], lty=1)
        lines(r0.3.Lo[r0.3.Lo < twoLoc.Lo.obOut] ~ sm[r0.3.Lo < twoLoc.Lo.obOut], lwd=2, col=COLS[3], lty=1)
        lines(r0.2.Hi[r0.2.Hi > twoLoc.Hi.obOut & r0.2.Hi <= 1] ~ sm[r0.2.Hi > twoLoc.Hi.obOut & r0.2.Hi <= 1], lwd=2, col=COLS[4], lty=1)
        lines(r0.2.Lo[r0.2.Lo < twoLoc.Lo.obOut] ~ sm[r0.2.Lo < twoLoc.Lo.obOut], lwd=2, col=COLS[4], lty=1)
        lines(r0.1.Hi[r0.1.Hi > twoLoc.Hi.obOut & r0.1.Hi <= 1] ~ sm[r0.1.Hi > twoLoc.Hi.obOut & r0.1.Hi <= 1], lwd=2, col=COLS[5], lty=1)
        lines(r0.1.Lo[r0.1.Lo < twoLoc.Lo.obOut] ~ sm[r0.1.Lo < twoLoc.Lo.obOut], lwd=2, col=COLS[5], lty=1)
        lines(r0.Hi[r0.Hi > twoLoc.Hi.obOut & r0.Hi <= 1] ~ sm[r0.Hi > twoLoc.Hi.obOut & r0.Hi <= 1], lwd=2, col=COLS[6], lty=1)
        lines(r0.Lo[r0.Lo < twoLoc.Lo.obOut] ~ sm[r0.Lo < twoLoc.Lo.obOut], lwd=2, col=COLS[6], lty=1)
        # Using only first eigenvalue (ignoring recombination)
        polygon(c(sm,rev(sm)), c(twoLoc.Hi.obOut, rev(twoLoc.Lo.obOut)), col=transparentColor('grey80', 0.6), border='grey70')
        lines(twoLoc.Hi.obOut[twoLoc.Hi.obOut <= 1] ~ sm[twoLoc.Hi.obOut <= 1], lwd=2, col='black')
        lines(twoLoc.Lo.obOut ~ sm, lwd=2, col='black')
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(h), " = 1/2")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.25, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = ', 0)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi.obOut)
        rm(twoLoc.Lo.obOut)
        rm(r0.2.Hi)
        rm(r0.2.Lo)
        rm(r0.1.Hi)
        rm(r0.1.Lo)
        rm(r0.Hi)
        rm(r0.Lo)



    ##  Panel Two: C = 0.25
        
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi   <-  inv.lab1.add(hf=0.5, hm=0.5, sm=sm, C=0.25)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.add(hf=0.5, hm=0.5, sm=sm, C=0.25)
        
        r0.5.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.5, C=0.25)
        r0.5.Hi[r0.5.Hi > 1]  <-  1.00000001
        r0.5.Hi[r0.5.Hi == 'NaN']  <-  1.00000001
        r0.5.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.5, C=0.25)

        r0.4.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.4, C=0.25)
        r0.4.Hi[r0.4.Hi > 1]  <-  1.00000001
        r0.4.Hi[r0.4.Hi == 'NaN']  <-  1.00000001
        r0.4.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.4, C=0.25)

        r0.3.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.3, C=0.25)
        r0.3.Hi[r0.3.Hi > 1]  <-  1.00000001
        r0.3.Hi[r0.3.Hi == 'NaN']  <-  1.00000001
        r0.3.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.3, C=0.25)

        r0.2.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.2, C=0.25)
        r0.2.Hi[r0.2.Hi > 1]  <-  1.00000001
        r0.2.Hi[r0.2.Hi == 'NaN']  <-  1.00000001
        r0.2.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.2, C=0.25)

        r0.1.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.1, C=0.25)
        r0.1.Hi[r0.1.Hi > 1]  <-  1.00000001
        r0.1.Hi[r0.1.Hi == 'NaN']  <-  1.00000001
        r0.1.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.1, C=0.25)

        r0.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0, C=0.25)
        r0.Hi[r0.Hi > 1]  <-  1.00000001
        r0.Hi[r0.Hi == 'NaN']  <-  1.00000001
        r0.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0, C=0.25)
        
          # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        #  w/ recombination
        lines(r0.5.Hi[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1] ~ sm[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1], lwd=2, col=COLS[1], lty=1)
        lines(r0.5.Lo[r0.5.Lo < twoLoc.Lo] ~ sm[r0.5.Lo < twoLoc.Lo], lwd=2, col=COLS[1], lty=1)
        lines(r0.4.Hi[r0.4.Hi > twoLoc.Hi & r0.4.Hi <= 1] ~ sm[r0.4.Hi > twoLoc.Hi & r0.4.Hi <= 1], lwd=2, col=COLS[2], lty=1)
        lines(r0.4.Lo[r0.4.Lo < twoLoc.Lo] ~ sm[r0.4.Lo < twoLoc.Lo], lwd=2, col=COLS[2], lty=1)
        lines(r0.3.Hi[r0.3.Hi > twoLoc.Hi & r0.3.Hi <= 1] ~ sm[r0.3.Hi > twoLoc.Hi & r0.3.Hi <= 1], lwd=2, col=COLS[3], lty=1)
        lines(r0.3.Lo[r0.3.Lo < twoLoc.Lo] ~ sm[r0.3.Lo < twoLoc.Lo], lwd=2, col=COLS[3], lty=1)
        lines(r0.2.Hi[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1] ~ sm[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1], lwd=2, col=COLS[4], lty=1)
        lines(r0.2.Lo[r0.2.Lo < twoLoc.Lo] ~ sm[r0.2.Lo < twoLoc.Lo], lwd=2, col=COLS[4], lty=1)
        lines(r0.1.Hi[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1] ~ sm[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1], lwd=2, col=COLS[5], lty=1)
        lines(r0.1.Lo[r0.1.Lo < twoLoc.Lo] ~ sm[r0.1.Lo < twoLoc.Lo], lwd=2, col=COLS[5], lty=1)
        lines(r0.Hi[r0.Hi > twoLoc.Hi & r0.Hi <= 1] ~ sm[r0.Hi > twoLoc.Hi & r0.Hi <= 1], lwd=2, col=COLS[6], lty=1)
        lines(r0.Lo[r0.Lo < twoLoc.Lo] ~ sm[r0.Lo < twoLoc.Lo], lwd=2, col=COLS[6], lty=1)
        # Using only first eigenvalue (ignoring recombination)
        polygon(c(sm,rev(sm)),c(twoLoc.Hi,rev(twoLoc.Lo)), col=transparentColor('grey80', 0.6), border='grey70')
        lines(twoLoc.Hi[twoLoc.Hi <= 1] ~ sm[twoLoc.Hi <= 1], lwd=2, col='black')
        lines(twoLoc.Lo ~ sm, lwd=2, col='black')
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = ',0.25)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.2.Hi)
        rm(r0.2.Lo)
        rm(r0.1.Hi)
        rm(r0.1.Lo)
        rm(r0.Hi)
        rm(r0.Lo)




    ##  Panel Three: C = 0.5
        
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi   <-  inv.lab1.add(hf=0.5, hm=0.5, sm=sm, C=0.5)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.add(hf=0.5, hm=0.5, sm=sm, C=0.5)
        
        r0.5.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.5, C=0.5)
        r0.5.Hi[r0.5.Hi > 1]  <-  1.00000001
        r0.5.Hi[r0.5.Hi == 'NaN']  <-  1.00000001
        r0.5.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.5, C=0.5)

        r0.4.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.4, C=0.5)
        r0.4.Hi[r0.4.Hi > 1]  <-  1.00000001
        r0.4.Hi[r0.4.Hi == 'NaN']  <-  1.00000001
        r0.4.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.4, C=0.5)

        r0.3.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.3, C=0.5)
        r0.3.Hi[r0.3.Hi > 1]  <-  1.00000001
        r0.3.Hi[r0.3.Hi == 'NaN']  <-  1.00000001
        r0.3.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.3, C=0.5)

        r0.2.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.2, C=0.5)
        r0.2.Hi[r0.2.Hi > 1]  <-  1.00000001
        r0.2.Hi[r0.2.Hi == 'NaN']  <-  1.00000001
        r0.2.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.2, C=0.5)

        r0.1.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.1, C=0.5)
        r0.1.Hi[r0.1.Hi > 1]  <-  1.00000001
        r0.1.Hi[r0.1.Hi == 'NaN']  <-  1.00000001
        r0.1.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.1, C=0.5)

        r0.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0, C=0.5)
        r0.Hi[r0.Hi > 1]  <-  1.00000001
        r0.Hi[r0.Hi == 'NaN']  <-  1.00000001
        r0.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0, C=0.5)

        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        lines(r0.5.Hi[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1] ~ sm[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1], lwd=2, col=COLS[1], lty=1)
        lines(r0.5.Lo[r0.5.Lo < twoLoc.Lo] ~ sm[r0.5.Lo < twoLoc.Lo], lwd=2, col=COLS[1], lty=1)
        lines(r0.4.Hi[r0.4.Hi > twoLoc.Hi & r0.4.Hi <= 1] ~ sm[r0.4.Hi > twoLoc.Hi & r0.4.Hi <= 1], lwd=2, col=COLS[2], lty=1)
        lines(r0.4.Lo[r0.4.Lo < twoLoc.Lo] ~ sm[r0.4.Lo < twoLoc.Lo], lwd=2, col=COLS[2], lty=1)
        lines(r0.3.Hi[r0.3.Hi > twoLoc.Hi & r0.3.Hi <= 1] ~ sm[r0.3.Hi > twoLoc.Hi & r0.3.Hi <= 1], lwd=2, col=COLS[3], lty=1)
        lines(r0.3.Lo[r0.3.Lo < twoLoc.Lo] ~ sm[r0.3.Lo < twoLoc.Lo], lwd=2, col=COLS[3], lty=1)
        lines(r0.2.Hi[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1] ~ sm[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1], lwd=2, col=COLS[4], lty=1)
        lines(r0.2.Lo[r0.2.Lo < twoLoc.Lo] ~ sm[r0.2.Lo < twoLoc.Lo], lwd=2, col=COLS[4], lty=1)
        lines(r0.1.Hi[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1] ~ sm[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1], lwd=2, col=COLS[5], lty=1)
        lines(r0.1.Lo[r0.1.Lo < twoLoc.Lo] ~ sm[r0.1.Lo < twoLoc.Lo], lwd=2, col=COLS[5], lty=1)
        lines(r0.Hi[r0.Hi > twoLoc.Hi & r0.Hi <= 1] ~ sm[r0.Hi > twoLoc.Hi & r0.Hi <= 1], lwd=2, col=COLS[6], lty=1)
        lines(r0.Lo[r0.Lo < twoLoc.Lo] ~ sm[r0.Lo < twoLoc.Lo], lwd=2, col=COLS[6], lty=1)

        # Using only first eigenvalue (ignoring recombination)
        polygon(c(sm,rev(sm)),c(twoLoc.Hi,rev(twoLoc.Lo)), col=transparentColor('grey80', 0.6), border='grey70')
        lines(twoLoc.Hi[twoLoc.Hi <= 1] ~ sm[twoLoc.Hi <= 1], lwd=2, col='black')
        lines(twoLoc.Lo ~ sm, lwd=2, col='black')
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = ',0.5)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        legend(
            x       =  usr[2]*0.39,
            y       =  usr[4],
        #    title   =  expression(paste(Outcome~of~invasion~analysis)),
            legend  =  c(
                        expression(paste(italic(r), " = ", 0)),
                        expression(paste(italic(r), " = ", 0.1)),
                        expression(paste(italic(r), " = ", 0.2)),
                        expression(paste(italic(r), " = ", 0.3)),
                        expression(paste(italic(r), " = ", 0.4)),
                        expression(paste(italic(r), " = ", 0.5)),
                        expression(paste(1~locus))),
            lty     =  1,
            lwd     =  3,
            col     =  c(rev(COLS),'black'),
            cex     =  0.75,
            xjust   =  1,
            yjust   =  1,
            bty     =  'n',
            border  =  NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.2.Hi)
        rm(r0.2.Lo)
        rm(r0.1.Hi)
        rm(r0.1.Lo)
        rm(r0.Hi)
        rm(r0.Lo)


############


##  Row 2: Dominance Reversal
    ##  Panel Four: C = 0
        
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi.obOut   <-  inv.lab1.domRev.obOut(hf=0.25, hm=0.25, sm=sm)
        twoLoc.Hi.obOut[twoLoc.Hi.obOut > 1]  <-  1.00000001
        twoLoc.Hi.obOut[10001]  <-  1.00000001
        twoLoc.Lo.obOut  <-  inv.lAB1.domRev.obOut(hf=0.25, hm=0.25, sm=sm)
        
        r0.5.Hi  <-  inv.lab2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0.5)
        r0.5.Hi[r0.5.Hi > 1]  <-  1.00000001
        r0.5.Hi[r0.5.Hi == 'NaN']  <-  1.00000001
        r0.5.Lo  <-  inv.lAB2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0.5)

        r0.4.Hi  <-  inv.lab2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0.4)
        r0.4.Hi[r0.4.Hi > 1]  <-  1.00000001
        r0.4.Hi[r0.4.Hi == 'NaN']  <-  1.00000001
        r0.4.Lo  <-  inv.lAB2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0.4)

        r0.3.Hi  <-  inv.lab2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0.3)
        r0.3.Hi[r0.3.Hi > 1]  <-  1.00000001
        r0.3.Hi[r0.3.Hi == 'NaN']  <-  1.00000001
        r0.3.Lo  <-  inv.lAB2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0.3)

        r0.2.Hi  <-  inv.lab2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0.2)
        r0.2.Hi[r0.2.Hi > 1]  <-  1.00000001
        r0.2.Hi[r0.2.Hi == 'NaN']  <-  1.00000001
        r0.2.Lo  <-  inv.lAB2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0.2)

        r0.1.Hi  <-  inv.lab2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0.1)
        r0.1.Hi[r0.1.Hi > 1]  <-  1.00000001
        r0.1.Hi[r0.1.Hi == 'NaN']  <-  1.00000001
        r0.1.Lo  <-  inv.lAB2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0.1)

        r0.Hi  <-  inv.lab2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0)
        r0.Hi[r0.Hi > 1]  <-  1.00000001
        r0.Hi[r0.Hi == 'NaN']  <-  1.00000001
        r0.Lo  <-  inv.lAB2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0)

        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        #  w/ recombination
        lines(r0.5.Hi[r0.5.Hi > twoLoc.Hi.obOut & r0.5.Hi <= 1] ~ sm[r0.5.Hi > twoLoc.Hi.obOut & r0.5.Hi <= 1], lwd=2, col=COLS[1], lty=1)
        lines(r0.5.Lo[r0.5.Lo < twoLoc.Lo.obOut] ~ sm[r0.5.Lo < twoLoc.Lo.obOut], lwd=2, col=COLS[1], lty=1)
        lines(r0.4.Hi[r0.4.Hi > twoLoc.Hi.obOut & r0.4.Hi <= 1] ~ sm[r0.4.Hi > twoLoc.Hi.obOut & r0.4.Hi <= 1], lwd=2, col=COLS[2], lty=1)
        lines(r0.4.Lo[r0.4.Lo < twoLoc.Lo.obOut] ~ sm[r0.4.Lo < twoLoc.Lo.obOut], lwd=2, col=COLS[2], lty=1)
        lines(r0.3.Hi[r0.3.Hi > twoLoc.Hi.obOut & r0.3.Hi <= 1] ~ sm[r0.3.Hi > twoLoc.Hi.obOut & r0.3.Hi <= 1], lwd=2, col=COLS[3], lty=1)
        lines(r0.3.Lo[r0.3.Lo < twoLoc.Lo.obOut] ~ sm[r0.3.Lo < twoLoc.Lo.obOut], lwd=2, col=COLS[3], lty=1)
        lines(r0.2.Hi[r0.2.Hi > twoLoc.Hi.obOut & r0.2.Hi <= 1] ~ sm[r0.2.Hi > twoLoc.Hi.obOut & r0.2.Hi <= 1], lwd=2, col=COLS[4], lty=1)
        lines(r0.2.Lo[r0.2.Lo < twoLoc.Lo.obOut] ~ sm[r0.2.Lo < twoLoc.Lo.obOut], lwd=2, col=COLS[4], lty=1)
        lines(r0.1.Hi[r0.1.Hi > twoLoc.Hi.obOut & r0.1.Hi <= 1] ~ sm[r0.1.Hi > twoLoc.Hi.obOut & r0.1.Hi <= 1], lwd=2, col=COLS[5], lty=1)
        lines(r0.1.Lo[r0.1.Lo < twoLoc.Lo.obOut] ~ sm[r0.1.Lo < twoLoc.Lo.obOut], lwd=2, col=COLS[5], lty=1)
        lines(r0.Hi[r0.Hi > twoLoc.Hi.obOut & r0.Hi <= 1] ~ sm[r0.Hi > twoLoc.Hi.obOut & r0.Hi <= 1], lwd=2, col=COLS[6], lty=1)
        lines(r0.Lo[r0.Lo < twoLoc.Lo.obOut] ~ sm[r0.Lo < twoLoc.Lo.obOut], lwd=2, col=COLS[6], lty=1)
        # Using only first eigenvalue (ignoring recombination)
        polygon(c(sm,rev(sm)), c(twoLoc.Hi.obOut, rev(twoLoc.Lo.obOut)), col=transparentColor('grey80', 0.6), border='grey70')
        lines(twoLoc.Hi.obOut[twoLoc.Hi.obOut <= 1] ~ sm[twoLoc.Hi.obOut <= 1], lwd=2, col='black')
        lines(twoLoc.Lo.obOut ~ sm, lwd=2, col='black')
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(h), " = 1/4")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.25, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.25, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi.obOut)
        rm(twoLoc.Lo.obOut)
        rm(r0.2.Hi)
        rm(r0.2.Lo)
        rm(r0.1.Hi)
        rm(r0.1.Lo)
        rm(r0.Hi)
        rm(r0.Lo)

    ##  Panel Five: C = 0.25
        
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi   <-  inv.lab1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.25)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.25)
        
        r0.5.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.5, C=0.25)
        r0.5.Hi[r0.5.Hi > 1]  <-  1.00000001
        r0.5.Hi[r0.5.Hi == 'NaN']  <-  1.00000001
        r0.5.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.5, C=0.25)

        r0.4.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.4, C=0.25)
        r0.4.Hi[r0.4.Hi > 1]  <-  1.00000001
        r0.4.Hi[r0.4.Hi == 'NaN']  <-  1.00000001
        r0.4.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.4, C=0.25)

        r0.3.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.3, C=0.25)
        r0.3.Hi[r0.3.Hi > 1]  <-  1.00000001
        r0.3.Hi[r0.3.Hi == 'NaN']  <-  1.00000001
        r0.3.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.3, C=0.25)

        r0.2.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.2, C=0.25)
        r0.2.Hi[r0.2.Hi > 1]  <-  1.00000001
        r0.2.Hi[r0.2.Hi == 'NaN']  <-  1.00000001
        r0.2.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.2, C=0.25)

        r0.1.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.1, C=0.25)
        r0.1.Hi[r0.1.Hi > 1]  <-  1.00000001
        r0.1.Hi[r0.1.Hi == 'NaN']  <-  1.00000001
        r0.1.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.1, C=0.25)

        r0.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0, C=0.25)
        r0.Hi[r0.Hi > 1]  <-  1.00000001
        r0.Hi[r0.Hi == 'NaN']  <-  1.00000001
        r0.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0, C=0.25)
        
          # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        lines(r0.5.Hi[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1] ~ sm[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1], lwd=2, col=COLS[1], lty=1)
        lines(r0.5.Lo[r0.5.Lo < twoLoc.Lo] ~ sm[r0.5.Lo < twoLoc.Lo], lwd=2, col=COLS[1], lty=1)
        lines(r0.4.Hi[r0.4.Hi > twoLoc.Hi & r0.4.Hi <= 1] ~ sm[r0.4.Hi > twoLoc.Hi & r0.4.Hi <= 1], lwd=2, col=COLS[2], lty=1)
        lines(r0.4.Lo[r0.4.Lo < twoLoc.Lo] ~ sm[r0.4.Lo < twoLoc.Lo], lwd=2, col=COLS[2], lty=1)
        lines(r0.3.Hi[r0.3.Hi > twoLoc.Hi & r0.3.Hi <= 1] ~ sm[r0.3.Hi > twoLoc.Hi & r0.3.Hi <= 1], lwd=2, col=COLS[3], lty=1)
        lines(r0.3.Lo[r0.3.Lo < twoLoc.Lo] ~ sm[r0.3.Lo < twoLoc.Lo], lwd=2, col=COLS[3], lty=1)
        lines(r0.2.Hi[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1] ~ sm[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1], lwd=2, col=COLS[4], lty=1)
        lines(r0.2.Lo[r0.2.Lo < twoLoc.Lo] ~ sm[r0.2.Lo < twoLoc.Lo], lwd=2, col=COLS[4], lty=1)
        lines(r0.1.Hi[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1] ~ sm[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1], lwd=2, col=COLS[5], lty=1)
        lines(r0.1.Lo[r0.1.Lo < twoLoc.Lo] ~ sm[r0.1.Lo < twoLoc.Lo], lwd=2, col=COLS[5], lty=1)
        lines(r0.Hi[r0.Hi > twoLoc.Hi & r0.Hi <= 1] ~ sm[r0.Hi > twoLoc.Hi & r0.Hi <= 1], lwd=2, col=COLS[6], lty=1)
        lines(r0.Lo[r0.Lo < twoLoc.Lo] ~ sm[r0.Lo < twoLoc.Lo], lwd=2, col=COLS[6], lty=1)
        # Using only first eigenvalue (ignoring recombination)
        polygon(c(sm,rev(sm)),c(twoLoc.Hi,rev(twoLoc.Lo)), col=transparentColor('grey80', 0.6), border='grey70')
        lines(twoLoc.Hi[twoLoc.Hi <= 1] ~ sm[twoLoc.Hi <= 1], lwd=2, col='black')
        lines(twoLoc.Lo ~ sm, lwd=2, col='black')
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, -0.25, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.05, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.2.Hi)
        rm(r0.2.Lo)
        rm(r0.1.Hi)
        rm(r0.1.Lo)
        rm(r0.Hi)
        rm(r0.Lo)




    ##  Panel Six: C = 0.5
        
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi   <-  inv.lab1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.5)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.5)
        
        r0.5.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.5, C=0.5)
        r0.5.Hi[r0.5.Hi > 1]  <-  1.00000001
        r0.5.Hi[r0.5.Hi == 'NaN']  <-  1.00000001
        r0.5.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.5, C=0.5)

        r0.4.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.4, C=0.5)
        r0.4.Hi[r0.4.Hi > 1]  <-  1.00000001
        r0.4.Hi[r0.4.Hi == 'NaN']  <-  1.00000001
        r0.4.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.4, C=0.5)

        r0.3.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.3, C=0.5)
        r0.3.Hi[r0.3.Hi > 1]  <-  1.00000001
        r0.3.Hi[r0.3.Hi == 'NaN']  <-  1.00000001
        r0.3.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.3, C=0.5)

        r0.2.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.2, C=0.5)
        r0.2.Hi[r0.2.Hi > 1]  <-  1.00000001
        r0.2.Hi[r0.2.Hi == 'NaN']  <-  1.00000001
        r0.2.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.2, C=0.5)

        r0.1.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.1, C=0.5)
        r0.1.Hi[r0.1.Hi > 1]  <-  1.00000001
        r0.1.Hi[r0.1.Hi == 'NaN']  <-  1.00000001
        r0.1.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.1, C=0.5)

        r0.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0, C=0.5)
        r0.Hi[r0.Hi > 1]  <-  1.00000001
        r0.Hi[r0.Hi == 'NaN']  <-  1.00000001
        r0.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0, C=0.5)
        
          # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        lines(r0.5.Hi[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1] ~ sm[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1], lwd=2, col=COLS[1], lty=1)
        lines(r0.5.Lo[r0.5.Lo < twoLoc.Lo] ~ sm[r0.5.Lo < twoLoc.Lo], lwd=2, col=COLS[1], lty=1)
        lines(r0.4.Hi[r0.4.Hi > twoLoc.Hi & r0.4.Hi <= 1] ~ sm[r0.4.Hi > twoLoc.Hi & r0.4.Hi <= 1], lwd=2, col=COLS[2], lty=1)
        lines(r0.4.Lo[r0.4.Lo < twoLoc.Lo] ~ sm[r0.4.Lo < twoLoc.Lo], lwd=2, col=COLS[2], lty=1)
        lines(r0.3.Hi[r0.3.Hi > twoLoc.Hi & r0.3.Hi <= 1] ~ sm[r0.3.Hi > twoLoc.Hi & r0.3.Hi <= 1], lwd=2, col=COLS[3], lty=1)
        lines(r0.3.Lo[r0.3.Lo < twoLoc.Lo] ~ sm[r0.3.Lo < twoLoc.Lo], lwd=2, col=COLS[3], lty=1)
        lines(r0.2.Hi[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1] ~ sm[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1], lwd=2, col=COLS[4], lty=1)
        lines(r0.2.Lo[r0.2.Lo < twoLoc.Lo] ~ sm[r0.2.Lo < twoLoc.Lo], lwd=2, col=COLS[4], lty=1)
        lines(r0.1.Hi[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1] ~ sm[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1], lwd=2, col=COLS[5], lty=1)
        lines(r0.1.Lo[r0.1.Lo < twoLoc.Lo] ~ sm[r0.1.Lo < twoLoc.Lo], lwd=2, col=COLS[5], lty=1)
        lines(r0.Hi[r0.Hi > twoLoc.Hi & r0.Hi <= 1] ~ sm[r0.Hi > twoLoc.Hi & r0.Hi <= 1], lwd=2, col=COLS[6], lty=1)
        lines(r0.Lo[r0.Lo < twoLoc.Lo] ~ sm[r0.Lo < twoLoc.Lo], lwd=2, col=COLS[6], lty=1)
        # Using only first eigenvalue (ignoring recombination)
        polygon(c(sm,rev(sm)),c(twoLoc.Hi,rev(twoLoc.Lo)), col=transparentColor('grey80', 0.6), border='grey70')
        lines(twoLoc.Hi[twoLoc.Hi <= 1] ~ sm[twoLoc.Hi <= 1], lwd=2, col='black')
        lines(twoLoc.Lo ~ sm, lwd=2, col='black')
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, -0.25, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.05, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.2.Hi)
        rm(r0.2.Lo)
        rm(r0.1.Hi)
        rm(r0.1.Lo)
        rm(r0.Hi)
        rm(r0.Lo)

}







#' Fig.1wk: Analytic results Kidwell Plots -- weak selection!
#' 
#'
#' @title Create a 2x2 panel plot with our analytic results showing Kidwell plots under weak selection
#' 
#' @export

Fig.1wk  <-  function() {

# Color scheme
    colfunc <- colorRampPalette(c("#252525", "grey70"))
    COLS  <-  colfunc(6)
#    COLS  <-  c("black", "#525252", "#737373", "#bdbdbd")

#  Create vector of male selection coefficiets for invasion functions
sm  <-  seq(0,0.1,by=0.00001)

# Set plot layout
layout.mat <- matrix(c(1,2,3,4,5,6), nrow=2, ncol=3, byrow=TRUE)
layout <- layout(layout.mat,respect=TRUE)

##  Row 1: Additive allelic effects

    ##  Panel One: C = 0
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi.obOut   <-  inv.lab1.add.obOut(hf=0.5, hm=0.5, sm=sm)
        twoLoc.Hi.obOut[twoLoc.Hi.obOut > 0.1]  <-  0.100000001
        twoLoc.Hi.obOut[10001]  <-  0.100000001
        twoLoc.Lo.obOut  <-  inv.lAB1.add.obOut(hf=0.5, hm=0.5, sm=sm)
        
        r0.Hi  <-  inv.lab2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0)
        r0.Hi[r0.Hi > 0.1]  <-  0.100000001
        r0.Hi[r0.Hi == 'NaN']  <-  0.100000001
        r0.Lo  <-  inv.lAB2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0)

        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.1), ylim = c(0,0.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        #  w/ recombination
        lines(r0.Hi[r0.Hi > twoLoc.Hi.obOut & r0.Hi <= 0.1] ~ sm[r0.Hi > twoLoc.Hi.obOut & r0.Hi <= 0.1], lwd=2, col=COLS[6], lty=1)
        lines(r0.Lo[r0.Lo < twoLoc.Lo.obOut] ~ sm[r0.Lo < twoLoc.Lo.obOut], lwd=2, col=COLS[6], lty=1)
        # Using only first eigenvalue (ignoring recombination)
        polygon(c(sm,rev(sm)), c(twoLoc.Hi.obOut, rev(twoLoc.Lo.obOut)), col=transparentColor('grey80', 0.6), border='grey70')
        lines(twoLoc.Hi.obOut[twoLoc.Hi.obOut <= 0.1] ~ sm[twoLoc.Hi.obOut <= 0.1], lwd=2, col='black')
        lines(twoLoc.Lo.obOut ~ sm, lwd=2, col='black')
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(h), " = 1/2")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.25, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = ', 0)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi.obOut)
        rm(twoLoc.Lo.obOut)
        rm(r0.Hi)
        rm(r0.Lo)



    ##  Panel Two: C = 0.25
        
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi   <-  inv.lab1.add(hf=0.5, hm=0.5, sm=sm, C=0.25)
        twoLoc.Hi[twoLoc.Hi > 0.1]  <-  0.100000001
        twoLoc.Lo  <-  inv.lAB1.add(hf=0.5, hm=0.5, sm=sm, C=0.25)
        
        r0.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0, C=0.25)
        r0.Hi[r0.Hi > 0.1]  <-  0.100000001
        r0.Hi[r0.Hi == 'NaN']  <-  0.100000001
        r0.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0, C=0.25)
        
          # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.1), ylim = c(0,0.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        #  w/ recombination
        lines(r0.Hi[r0.Hi > twoLoc.Hi & r0.Hi <= 0.1] ~ sm[r0.Hi > twoLoc.Hi & r0.Hi <= 0.1], lwd=2, col=COLS[6], lty=1)
        lines(r0.Lo[r0.Lo < twoLoc.Lo] ~ sm[r0.Lo < twoLoc.Lo], lwd=2, col=COLS[6], lty=1)
        # Using only first eigenvalue (ignoring recombination)
        polygon(c(sm,rev(sm)),c(twoLoc.Hi,rev(twoLoc.Lo)), col=transparentColor('grey80', 0.6), border='grey70')
        lines(twoLoc.Hi[twoLoc.Hi <= 0.1] ~ sm[twoLoc.Hi <= 0.1], lwd=2, col='black')
        lines(twoLoc.Lo ~ sm, lwd=2, col='black')
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = ',0.25)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.Hi)
        rm(r0.Lo)




    ##  Panel Three: C = 0.5
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi   <-  inv.lab1.add(hf=0.5, hm=0.5, sm=sm, C=0.5)
        twoLoc.Hi[twoLoc.Hi > 0.1]  <-  0.100000001
        twoLoc.Lo  <-  inv.lAB1.add(hf=0.5, hm=0.5, sm=sm, C=0.5)
        
        r0.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0, C=0.5)
        r0.Hi[r0.Hi > 0.1]  <-  0.100000001
        r0.Hi[r0.Hi == 'NaN']  <-  0.100000001
        r0.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0, C=0.5)

          # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.1), ylim = c(0,0.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        lines(r0.Hi[r0.Hi > twoLoc.Hi & r0.Hi <= 0.1] ~ sm[r0.Hi > twoLoc.Hi & r0.Hi <= 0.1], lwd=2, col=COLS[6], lty=1)
        lines(r0.Lo[r0.Lo < twoLoc.Lo] ~ sm[r0.Lo < twoLoc.Lo], lwd=2, col=COLS[6], lty=1)
        # Using only first eigenvalue (ignoring recombination)
        polygon(c(sm,rev(sm)),c(twoLoc.Hi,rev(twoLoc.Lo)), col=transparentColor('grey80', 0.6), border='grey70')
        lines(twoLoc.Hi[twoLoc.Hi <= 0.1] ~ sm[twoLoc.Hi <= 0.1], lwd=2, col='black')
        lines(twoLoc.Lo ~ sm, lwd=2, col='black')
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(C), ' = ',0.5)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        legend(
            x       =  usr[2]*0.39,
            y       =  usr[4],
        #    title   =  expression(paste(Outcome~of~invasion~analysis)),
            legend  =  c(
                        expression(paste(italic(r), " = ", 0)),
                        expression(paste(1~locus))),
            lty     =  1,
            lwd     =  3,
            col     =  c(rev(COLS)[1],'black'),
            cex     =  0.75,
            xjust   =  1,
            yjust   =  1,
            bty     =  'n',
            border  =  NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.Hi)
        rm(r0.Lo)


############


##  Row 2: Dominance Reversal
    ##  Panel Four: C = 0
        
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi.obOut   <-  inv.lab1.domRev.obOut(hf=0.25, hm=0.25, sm=sm)
        twoLoc.Hi.obOut[twoLoc.Hi.obOut > 0.1]  <-  0.100000001
        twoLoc.Hi.obOut[10001]  <-  0.100000001
        twoLoc.Lo.obOut  <-  inv.lAB1.domRev.obOut(hf=0.25, hm=0.25, sm=sm)
        
        r0.Hi  <-  inv.lab2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0)
        r0.Hi[r0.Hi > 0.1]  <-  0.100000001
        r0.Hi[r0.Hi == 'NaN']  <-  0.100000001
        r0.Lo  <-  inv.lAB2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0)

        # Make the plot
#        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.1), ylim = c(0,0.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        #  w/ recombination
        lines(r0.Hi[r0.Hi > twoLoc.Hi.obOut & r0.Hi <= 0.1] ~ sm[r0.Hi > twoLoc.Hi.obOut & r0.Hi <= 0.1], lwd=2, col=COLS[6], lty=1)
        lines(r0.Lo[r0.Lo < twoLoc.Lo.obOut] ~ sm[r0.Lo < twoLoc.Lo.obOut], lwd=2, col=COLS[6], lty=1)
        # Using only first eigenvalue (ignoring recombination)
        polygon(c(sm,rev(sm)), c(twoLoc.Hi.obOut, rev(twoLoc.Lo.obOut)), col=transparentColor('grey80', 0.6), border='grey70')
        lines(twoLoc.Hi.obOut[twoLoc.Hi.obOut <= 0.1] ~ sm[twoLoc.Hi.obOut <= 0.1], lwd=2, col='black')
        lines(twoLoc.Lo.obOut ~ sm, lwd=2, col='black')
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(h), " = 1/4")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.25, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.25, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi.obOut)
        rm(twoLoc.Lo.obOut)
        rm(r0.Hi)
        rm(r0.Lo)

    ##  Panel Five: C = 0.25
        
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi   <-  inv.lab1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.25)
        twoLoc.Hi[twoLoc.Hi > 0.1]  <-  0.100000001
        twoLoc.Hi[10001]  <-  0.100000001
        twoLoc.Lo  <-  inv.lAB1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.25)
        
        r0.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0, C=0.25)
        r0.Hi[r0.Hi > 0.1]  <-  0.100000001
        r0.Hi[r0.Hi == 'NaN']  <-  0.100000001
        r0.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0, C=0.25)
        
          # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.1), ylim = c(0,0.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        lines(r0.Hi[r0.Hi > twoLoc.Hi & r0.Hi <= 0.1] ~ sm[r0.Hi > twoLoc.Hi & r0.Hi <= 0.1], lwd=2, col=COLS[6], lty=1)
        lines(r0.Lo[r0.Lo < twoLoc.Lo] ~ sm[r0.Lo < twoLoc.Lo], lwd=2, col=COLS[6], lty=1)
        # Using only first eigenvalue (ignoring recombination)
        polygon(c(sm,rev(sm)),c(twoLoc.Hi,rev(twoLoc.Lo)), col=transparentColor('grey80', 0.6), border='grey70')
        lines(twoLoc.Hi[twoLoc.Hi <= 0.1] ~ sm[twoLoc.Hi <= 0.1], lwd=2, col='black')
        lines(twoLoc.Lo ~ sm, lwd=2, col='black')
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, -0.25, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.05, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.Hi)
        rm(r0.Lo)




    ##  Panel Six: C = 0.5
        
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi   <-  inv.lab1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.5)
        twoLoc.Hi[twoLoc.Hi > 0.1]  <-  0.100000001
        twoLoc.Lo  <-  inv.lAB1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.5)
        
        r0.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0, C=0.5)
        r0.Hi[r0.Hi > 0.1]  <-  0.100000001
        r0.Hi[r0.Hi == 'NaN']  <-  0.100000001
        r0.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0, C=0.5)
        
          # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.1), ylim = c(0,0.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        lines(r0.Hi[r0.Hi > twoLoc.Hi & r0.Hi <= 0.1] ~ sm[r0.Hi > twoLoc.Hi & r0.Hi <= 0.1], lwd=2, col=COLS[6], lty=1)
        lines(r0.Lo[r0.Lo < twoLoc.Lo] ~ sm[r0.Lo < twoLoc.Lo], lwd=2, col=COLS[6], lty=1)
        # Using only first eigenvalue (ignoring recombination)
        polygon(c(sm,rev(sm)),c(twoLoc.Hi,rev(twoLoc.Lo)), col=transparentColor('grey80', 0.6), border='grey70')
        lines(twoLoc.Hi[twoLoc.Hi <= 0.1] ~ sm[twoLoc.Hi <= 0.1], lwd=2, col='black')
        lines(twoLoc.Lo ~ sm, lwd=2, col='black')
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, -0.25, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.05, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.Hi)
        rm(r0.Lo)

}






#' Fig.2: polymorphism ~ recombination rate plots.
#' 
#'
#' @title Create a 2x2 panel plot with our analytic results showing Kidwell plots
#' 
#' @export

Fig.2  <- function() {

    ## Read data files for plotting

        # Additive Effects
        C0.0.h.5  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0_hf0.5_hm0.5_n30000.txt', head=TRUE)
        C0.1.h.5  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0.1_hf0.5_hm0.5_n30000.txt', head=TRUE)
        C0.2.h.5  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0.2_hf0.5_hm0.5_n30000.txt', head=TRUE)
        C0.3.h.5  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0.3_hf0.5_hm0.5_n30000.txt', head=TRUE)
        C0.4.h.5  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0.4_hf0.5_hm0.5_n30000.txt', head=TRUE)
        C0.5.h.5  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0.5_hf0.5_hm0.5_n30000.txt', head=TRUE)
        C0.6.h.5  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0.6_hf0.5_hm0.5_n30000.txt', head=TRUE)
        C0.7.h.5  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0.7_hf0.5_hm0.5_n30000.txt', head=TRUE)
        C0.8.h.5  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0.8_hf0.5_hm0.5_n30000.txt', head=TRUE)
        C0.9.h.5  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0.9_hf0.5_hm0.5_n30000.txt', head=TRUE)

        # Dominance Reversal
        C0.0.h.25  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0_hf0.25_hm0.25_n30000.txt', head=TRUE)
        C0.1.h.25  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0.1_hf0.25_hm0.25_n30000.txt', head=TRUE)
        C0.2.h.25  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0.2_hf0.25_hm0.25_n30000.txt', head=TRUE)
        C0.3.h.25  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0.3_hf0.25_hm0.25_n30000.txt', head=TRUE)
        C0.4.h.25  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0.4_hf0.25_hm0.25_n30000.txt', head=TRUE)
        C0.5.h.25  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0.5_hf0.25_hm0.25_n30000.txt', head=TRUE)
        C0.6.h.25  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0.6_hf0.25_hm0.25_n30000.txt', head=TRUE)
        C0.7.h.25  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0.7_hf0.25_hm0.25_n30000.txt', head=TRUE)
        C0.8.h.25  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0.8_hf0.25_hm0.25_n30000.txt', head=TRUE)
        C0.9.h.25  <-  read.table('./output/data/Fig2Data/TESTpropPrp.out_C0.9_hf0.25_hm0.25_n30000.txt', head=TRUE)


    # Color scheme
    colfunc <- colorRampPalette(c("grey60", "black"))
    COLS  <-  colfunc(10)
#    COLS  <-  c("#000000", "#252525", "#525252", "#737373", "#969696")

    # Set plot layout
    layout.mat <- matrix(c(1,2,3,4), nrow=2, ncol=2, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

    ####################
    ##  ADDITIVE EFFECTS
    ##  Panel 1:  Proportion of parameter space resulting in Protected Polymorphism
        par(omi=rep(0.75, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0.2,0.55), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
        plotGrid(lineCol='white')
        box()
        lines(C0.0.h.5$PrP ~ C0.0.h.5$r, lwd=2, lty=1, col=COLS[1])
        lines(C0.1.h.5$PrP ~ C0.1.h.5$r, lwd=2, lty=1, col=COLS[2])
        lines(C0.2.h.5$PrP ~ C0.2.h.5$r, lwd=2, lty=1, col=COLS[3])
        lines(C0.3.h.5$PrP ~ C0.3.h.5$r, lwd=2, lty=1, col=COLS[4])
        lines(C0.4.h.5$PrP ~ C0.4.h.5$r, lwd=2, lty=1, col=COLS[5])
        lines(C0.5.h.5$PrP ~ C0.5.h.5$r, lwd=2, lty=1, col=COLS[6])
        lines(C0.6.h.5$PrP ~ C0.6.h.5$r, lwd=2, lty=1, col=COLS[7])
        lines(C0.7.h.5$PrP ~ C0.7.h.5$r, lwd=2, lty=1, col=COLS[8])
        lines(C0.8.h.5$PrP ~ C0.8.h.5$r, lwd=2, lty=1, col=COLS[9])
        lines(C0.9.h.5$PrP ~ C0.9.h.5$r, lwd=2, lty=1, col=COLS[10])
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.5, 1.15, 'Protected polymorphism', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(h), " = 1/2")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.25, 0.5, expression(paste("Proportion parameters space")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        legend(
            x       =  usr[2]*0.6,
            y       =  usr[4],
        #    title   =  expression(paste(Outcome~of~invasion~analysis)),
            legend  =  c(
                        expression(paste(italic(C), " = ", 0.0)),
                        expression(paste(italic(C), " = ", 0.1)),
                        expression(paste(italic(C), " = ", 0.2)),
                        expression(paste(italic(C), " = ", 0.3)),
                        expression(paste(italic(C), " = ", 0.4))),
            lty     =  1,
            lwd     =  3,
            col     =  COLS[1:5],
            cex     =  0.75,
            xjust   =  1,
            yjust   =  1,
            bty     =  'n',
            border  =  NA)
        legend(
            x       =  usr[2]*0.98,
            y       =  usr[4],
        #    title   =  expression(paste(Outcome~of~invasion~analysis)),
            legend  =  c(
                        expression(paste(italic(C), " = ", 0.5)),
                        expression(paste(italic(C), " = ", 0.6)),
                        expression(paste(italic(C), " = ", 0.7)),
                        expression(paste(italic(C), " = ", 0.8)),
                        expression(paste(italic(C), " = ", 0.9))),
            lty     =  1,
            lwd     =  3,
            col     =  COLS[6:10],
            cex     =  0.75,
            xjust   =  1,
            yjust   =  1,
            bty     =  'n',
            border  =  NA)

#    ##  Panel 2: Increase in parameter space resulting in polymorphism due to recombination
#        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,0.2), ylab='', xlab='', cex.lab=1.2)
#        usr  <-  par('usr')
#        rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
#        plotGrid(lineCol='white')
#        box()
#        lines(C0.0.h.5$rPrP ~ C0.0.h.5$r, lwd=2, lty=1, col=COLS[1])
#        lines(C0.1.h.5$rPrP ~ C0.1.h.5$r, lwd=2, lty=1, col=COLS[2])
#        lines(C0.2.h.5$rPrP ~ C0.2.h.5$r, lwd=2, lty=1, col=COLS[3])
#        lines(C0.3.h.5$rPrP ~ C0.3.h.5$r, lwd=2, lty=1, col=COLS[4])
#        lines(C0.4.h.5$rPrP ~ C0.4.h.5$r, lwd=2, lty=1, col=COLS[5])
#        lines(C0.5.h.5$rPrP ~ C0.5.h.5$r, lwd=2, lty=1, col=COLS[6])
#        lines(C0.6.h.5$rPrP ~ C0.6.h.5$r, lwd=2, lty=1, col=COLS[7])
#        lines(C0.7.h.5$rPrP ~ C0.7.h.5$r, lwd=2, lty=1, col=COLS[8])
#        lines(C0.8.h.5$rPrP ~ C0.8.h.5$r, lwd=2, lty=1, col=COLS[9])
#        lines(C0.9.h.5$rPrP ~ C0.9.h.5$r, lwd=2, lty=1, col=COLS[10])
#        axis(1, las=1, labels=NA)
#        axis(2, las=1)
#        proportionalLabel(0.5, 1.15, 'Poly. due to recombination', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.03, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)


    ##  Panel 3: Proportional increase in parameter space relative to single locus case
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(1,2.25), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
        plotGrid(lineCol='white')
        box()
        lines(C0.0.h.5$PrP / C0.0.h.5$slPrP ~ C0.0.h.5$r, lwd=2, lty=1, col=COLS[1])
        lines(C0.1.h.5$PrP / C0.1.h.5$slPrP ~ C0.1.h.5$r, lwd=2, lty=1, col=COLS[2])
        lines(C0.2.h.5$PrP / C0.2.h.5$slPrP ~ C0.2.h.5$r, lwd=2, lty=1, col=COLS[3])
        lines(C0.3.h.5$PrP / C0.3.h.5$slPrP ~ C0.3.h.5$r, lwd=2, lty=1, col=COLS[4])
        lines(C0.4.h.5$PrP / C0.4.h.5$slPrP ~ C0.4.h.5$r, lwd=2, lty=1, col=COLS[5])
        lines(C0.5.h.5$PrP / C0.5.h.5$slPrP ~ C0.5.h.5$r, lwd=2, lty=1, col=COLS[6])
        lines(C0.6.h.5$PrP / C0.6.h.5$slPrP ~ C0.6.h.5$r, lwd=2, lty=1, col=COLS[7])
        lines(C0.7.h.5$PrP / C0.7.h.5$slPrP ~ C0.7.h.5$r, lwd=2, lty=1, col=COLS[8])
        lines(C0.8.h.5$PrP / C0.8.h.5$slPrP ~ C0.8.h.5$r, lwd=2, lty=1, col=COLS[9])
        lines(C0.9.h.5$PrP / C0.9.h.5$slPrP ~ C0.9.h.5$r, lwd=2, lty=1, col=COLS[10])
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.5, 1.15, 'Prop. increase in Poly.', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.03, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)


    ######################
    ##  DOMINANCE REVERSAL
    ##  Panel 4:  Proportion of parameter space resulting in Protected Polymorphism
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0.2,0.8), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
        plotGrid(lineCol='white')
        box()
        lines(C0.0.h.25$PrP ~ C0.0.h.25$r, lwd=2, lty=1, col=COLS[1])
        lines(C0.1.h.25$PrP ~ C0.1.h.25$r, lwd=2, lty=1, col=COLS[2])
        lines(C0.2.h.25$PrP ~ C0.2.h.25$r, lwd=2, lty=1, col=COLS[3])
        lines(C0.3.h.25$PrP ~ C0.3.h.25$r, lwd=2, lty=1, col=COLS[4])
        lines(C0.4.h.25$PrP ~ C0.4.h.25$r, lwd=2, lty=1, col=COLS[5])
        lines(C0.5.h.25$PrP ~ C0.5.h.25$r, lwd=2, lty=1, col=COLS[6])
        lines(C0.6.h.25$PrP ~ C0.6.h.25$r, lwd=2, lty=1, col=COLS[7])
        lines(C0.7.h.25$PrP ~ C0.7.h.25$r, lwd=2, lty=1, col=COLS[8])
        lines(C0.8.h.25$PrP ~ C0.8.h.25$r, lwd=2, lty=1, col=COLS[9])
        lines(C0.9.h.25$PrP ~ C0.9.h.25$r, lwd=2, lty=1, col=COLS[10])
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(0.03, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(h), " = 1/4")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.25, 0.5, expression(paste("Proportion parameters space")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.25, expression(paste(italic(r))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)

#    ##  Panel 5: Increase in parameter space resulting in polymorphism due to recombination
#        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(0,0.2), ylab='', xlab='', cex.lab=1.2)
#        usr  <-  par('usr')
#        rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
#        plotGrid(lineCol='white')
#        box()
#        lines(C0.0.h.25$rPrP ~ C0.0.h.25$r, lwd=2, lty=1, col=COLS[1])
#        lines(C0.1.h.25$rPrP ~ C0.1.h.25$r, lwd=2, lty=1, col=COLS[2])
#        lines(C0.2.h.25$rPrP ~ C0.2.h.25$r, lwd=2, lty=1, col=COLS[3])
#        lines(C0.3.h.25$rPrP ~ C0.3.h.25$r, lwd=2, lty=1, col=COLS[4])
#        lines(C0.4.h.25$rPrP ~ C0.4.h.25$r, lwd=2, lty=1, col=COLS[5])
#        lines(C0.5.h.25$rPrP ~ C0.5.h.25$r, lwd=2, lty=1, col=COLS[6])
#        lines(C0.6.h.25$rPrP ~ C0.6.h.25$r, lwd=2, lty=1, col=COLS[7])
#        lines(C0.7.h.25$rPrP ~ C0.7.h.25$r, lwd=2, lty=1, col=COLS[8])
#        lines(C0.8.h.25$rPrP ~ C0.8.h.25$r, lwd=2, lty=1, col=COLS[9])
#        lines(C0.9.h.25$rPrP ~ C0.9.h.25$r, lwd=2, lty=1, col=COLS[10])
#        axis(1, las=1)
#        axis(2, las=1)
#        proportionalLabel(0.03, 1.05, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, -0.25, expression(paste(italic(r))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)


    ##  Panel 6: Proportional increase in parameter space relative to single locus case
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,0.5), ylim = c(1,2.25), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
        plotGrid(lineCol='white')
        box()
        lines(C0.0.h.25$PrP / C0.0.h.25$slPrP ~ C0.0.h.25$r, lwd=2, lty=1, col=COLS[1])
        lines(C0.1.h.25$PrP / C0.1.h.25$slPrP ~ C0.1.h.25$r, lwd=2, lty=1, col=COLS[2])
        lines(C0.2.h.25$PrP / C0.2.h.25$slPrP ~ C0.2.h.25$r, lwd=2, lty=1, col=COLS[3])
        lines(C0.3.h.25$PrP / C0.3.h.25$slPrP ~ C0.3.h.25$r, lwd=2, lty=1, col=COLS[4])
        lines(C0.4.h.25$PrP / C0.4.h.25$slPrP ~ C0.4.h.25$r, lwd=2, lty=1, col=COLS[5])
        lines(C0.5.h.25$PrP / C0.5.h.25$slPrP ~ C0.5.h.25$r, lwd=2, lty=1, col=COLS[6])
        lines(C0.6.h.25$PrP / C0.6.h.25$slPrP ~ C0.6.h.25$r, lwd=2, lty=1, col=COLS[7])
        lines(C0.7.h.25$PrP / C0.7.h.25$slPrP ~ C0.7.h.25$r, lwd=2, lty=1, col=COLS[8])
        lines(C0.8.h.25$PrP / C0.8.h.25$slPrP ~ C0.8.h.25$r, lwd=2, lty=1, col=COLS[9])
        lines(C0.9.h.25$PrP / C0.9.h.25$slPrP ~ C0.9.h.25$r, lwd=2, lty=1, col=COLS[10])
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(0.03, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.25, expression(paste(italic(r))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)

}











#' Fig.S1: Supplementary figure showing comparison between deterministic
#'         recursion simulations and invasion analysis based on eigenvalues
#' 
#'
#' @title Fig.S1: Supplementary figure showing comparison between deterministic
#'                recursion simulations and invasion analysis based on eigenvalues
#' @export

Fig.S1_add  <-  function() {

    ## Read data files for plotting
        C0.0.h.5   <-  read.table('./output/data/simResults/recFwdSimLoop.out_C0_hf0.5_hm0.5_sMax1.txt', head=TRUE)
        C0.25.h.5  <-  read.table('./output/data/simResults/recFwdSimLoop.out_C0.25_hf0.5_hm0.5_sMax1.txt', head=TRUE)
        C0.5.h.5   <-  read.table('./output/data/simResults/recFwdSimLoop.out_C0.5_hf0.5_hm0.5_sMax1.txt', head=TRUE)
        C0.75.h.5  <-  read.table('./output/data/simResults/recFwdSimLoop.out_C0.75_hf0.5_hm0.5_sMax1.txt', head=TRUE)

    # Color scheme
    COLS  <-  c(transparentColor('seagreen3', opacity=0.2), 
                transparentColor('dodgerblue2', opacity=0.2), 
                transparentColor('tomato2', opacity=0.2),
                'black')

    #  Create vector of male selection coefficients for invasion functions
    sm  <-  seq(0,1,by=0.0001)

    # Set plot layout
    layout.mat <- matrix(c(1:16), nrow=4, ncol=4, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

##  Row 1: r = 0
    ##  Panel One: C = 0
        
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi.obOut   <-  inv.lab1.add.obOut(hf=0.5, hm=0.5, sm=sm)
        twoLoc.Hi.obOut[twoLoc.Hi.obOut > 1]  <-  1.00000001
        twoLoc.Hi.obOut[10001]  <-  1.00000001
        twoLoc.Lo.obOut  <-  inv.lAB1.add.obOut(hf=0.5, hm=0.5, sm=sm)
        
        r0.Hi  <-  inv.lab2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0)
        r0.Hi[r0.Hi > 1]  <-  1.00000001
        r0.Hi[r0.Hi == 'NaN']  <-  1.00000001
        r0.Lo  <-  inv.lAB2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0)

        pAgree  <-  rounded(sum(C0.0.h.5$agree[C0.0.h.5$r == 0.0])/length(C0.0.h.5$agree[C0.0.h.5$r == 0.0]), precision=3)
        pSim    <-  rounded(length(C0.0.h.5$sf[C0.0.h.5$Poly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.0]) / 
                            length(C0.0.h.5$sf[C0.0.h.5$r  ==  0.0]), precision=3)
        pEig    <-  rounded(length(C0.0.h.5$sf[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.0]) / 
                            length(C0.0.h.5$sf[C0.0.h.5$r  ==  0.0]), precision=3)

        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.0.h.5$sf[C0.0.h.5$Poly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.0] ~
               C0.0.h.5$sm[C0.0.h.5$Poly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.0.h.5$sf[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.0] ~
               C0.0.h.5$sm[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[5])
        points(C0.0.h.5$sf[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 1 & C0.0.h.5$r  ==  0.0] ~
               C0.0.h.5$sm[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 1 & C0.0.h.5$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.Hi[r0.Hi > twoLoc.Hi.obOut & r0.Hi <= 1] ~ sm[r0.Hi > twoLoc.Hi.obOut & r0.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.Lo[r0.Lo < twoLoc.Lo.obOut] ~ sm[r0.Lo < twoLoc.Lo.obOut], lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.5, 1.25, expression(paste(italic(C), ' = ', 0)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = 0")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.05, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.5, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi.obOut)
        rm(twoLoc.Lo.obOut)
        rm(r0.Hi)
        rm(r0.Lo)


    ##  Panel Two: C = 0.25
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.add(hf=0.5, hm=0.5, sm=sm, C=0.25)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.add(hf=0.5, hm=0.5, sm=sm, C=0.25)
        
        r0.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0, C=0.25)
        r0.Hi[r0.Hi > 1]  <-  1.00000001
        r0.Hi[r0.Hi == 'NaN']  <-  1.00000001
        r0.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0, C=0.25)

        pAgree  <-  rounded(sum(C0.25.h.5$agree[C0.25.h.5$r == 0.0])/length(C0.25.h.5$agree[C0.25.h.5$r == 0.0]), precision=3)
        pSim    <-  rounded(length(C0.25.h.5$sf[C0.25.h.5$Poly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.0]) / 
                            length(C0.25.h.5$sf[C0.25.h.5$r  ==  0.0]), precision=3)
        pEig    <-  rounded(length(C0.25.h.5$sf[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.0]) / 
                            length(C0.25.h.5$sf[C0.25.h.5$r  ==  0.0]), precision=3)
        
        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.25.h.5$sf[C0.25.h.5$Poly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.0] ~
               C0.25.h.5$sm[C0.25.h.5$Poly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.25.h.5$sf[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.0] ~
               C0.25.h.5$sm[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.25.h.5$sf[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 1 & C0.25.h.5$r  ==  0.0] ~
               C0.25.h.5$sm[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 1 & C0.25.h.5$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.Hi[r0.Hi > twoLoc.Hi & r0.Hi <= 1] ~ sm[r0.Hi > twoLoc.Hi & r0.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.Lo[r0.Lo < twoLoc.Lo] ~ sm[r0.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.25, expression(paste(italic(C), ' = ',0.25)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.5, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.Hi)
        rm(r0.Lo)



    ##  Panel Three: C = 0.5
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.add(hf=0.5, hm=0.5, sm=sm, C=0.5)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.add(hf=0.5, hm=0.5, sm=sm, C=0.5)
        
        r0.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0, C=0.5)
        r0.Hi[r0.Hi > 1]  <-  1.00000001
        r0.Hi[r0.Hi == 'NaN']  <-  1.00000001
        r0.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0, C=0.5)

        pAgree  <-  rounded(sum(C0.5.h.5$agree[C0.5.h.5$r == 0.0])/length(C0.5.h.5$agree[C0.5.h.5$r == 0.0]), precision=3)
        pSim    <-  rounded(length(C0.5.h.5$sf[C0.5.h.5$Poly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.0]) / 
                            length(C0.5.h.5$sf[C0.5.h.5$r  ==  0.0]), precision=3)
        pEig    <-  rounded(length(C0.5.h.5$sf[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.0]) / 
                            length(C0.5.h.5$sf[C0.5.h.5$r  ==  0.0]), precision=3)
        
        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.5.h.5$sf[C0.5.h.5$Poly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.0] ~
               C0.5.h.5$sm[C0.5.h.5$Poly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.5.h.5$sf[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.0] ~
               C0.5.h.5$sm[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.5.h.5$sf[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 1 & C0.5.h.5$r  ==  0.0] ~
               C0.5.h.5$sm[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 1 & C0.5.h.5$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.Hi[r0.Hi > twoLoc.Hi & r0.Hi <= 1] ~ sm[r0.Hi > twoLoc.Hi & r0.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.Lo[r0.Lo < twoLoc.Lo] ~ sm[r0.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.25, expression(paste(italic(C), ' = ',0.5)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.5, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.Hi)
        rm(r0.Lo)


    ##  Panel Four: C = 0.75
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.add(hf=0.5, hm=0.5, sm=sm, C=0.75)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.add(hf=0.5, hm=0.5, sm=sm, C=0.75)
        
        r0.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0, C=0.75)
        r0.Hi[r0.Hi > 1]  <-  1.00000001
        r0.Hi[r0.Hi == 'NaN']  <-  1.00000001
        r0.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0, C=0.75)

        pAgree  <-  rounded(sum(C0.75.h.5$agree[C0.75.h.5$r == 0.0])/length(C0.75.h.5$agree[C0.75.h.5$r == 0.0]), precision=3)
        pSim    <-  rounded(length(C0.75.h.5$sf[C0.75.h.5$Poly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.0]) / 
                            length(C0.75.h.5$sf[C0.75.h.5$r  ==  0.0]), precision=3)
        pEig    <-  rounded(length(C0.75.h.5$sf[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.0]) / 
                            length(C0.75.h.5$sf[C0.75.h.5$r  ==  0.0]), precision=3)
        
        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.75.h.5$sf[C0.75.h.5$Poly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.0] ~
               C0.75.h.5$sm[C0.75.h.5$Poly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.75.h.5$sf[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.0] ~
               C0.75.h.5$sm[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.75.h.5$sf[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 1 & C0.75.h.5$r  ==  0.0] ~
               C0.75.h.5$sm[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 1 & C0.75.h.5$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[1])
        points(0.02,0.99, pch=21, col=NA, cex=0.8, bg='seagreen3')
        points(0.02,0.91, pch=21, col=NA, cex=0.8, bg='tomato2')
        points(0.02,0.83, pch=21, col=NA, cex=0.8, bg='dodgerblue2')
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.Hi[r0.Hi > twoLoc.Hi & r0.Hi <= 1] ~ sm[r0.Hi > twoLoc.Hi & r0.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.Lo[r0.Lo < twoLoc.Lo] ~ sm[r0.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.25, expression(paste(italic(C), ' = ',0.75)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.6, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.88, substitute(p~" Sim.", list(p = pSim)),    cex=0.6, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.80, substitute(p~" Eig.", list(p = pEig)),    cex=0.6, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.Hi)
        rm(r0.Lo)



##  Row 2: r = 0.1
    ##  Panel 5: C = 0
        
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi.obOut   <-  inv.lab1.add.obOut(hf=0.5, hm=0.5, sm=sm)
        twoLoc.Hi.obOut[twoLoc.Hi.obOut > 1]  <-  1.00000001
        twoLoc.Hi.obOut[10001]  <-  1.00000001
        twoLoc.Lo.obOut  <-  inv.lAB1.add.obOut(hf=0.5, hm=0.5, sm=sm)
        
        r0.1.Hi  <-  inv.lab2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0.1)
        r0.1.Hi[r0.1.Hi > 1]  <-  1.00000001
        r0.1.Hi[r0.1.Hi == 'NaN']  <-  1.00000001
        r0.1.Lo  <-  inv.lAB2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0.1)

        pAgree  <-  rounded(sum(C0.0.h.5$agree[C0.0.h.5$r == 0.1])/length(C0.0.h.5$agree[C0.0.h.5$r == 0.1]), precision=3)
        pSim    <-  rounded(length(C0.0.h.5$sf[C0.0.h.5$Poly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.1]) / 
                            length(C0.0.h.5$sf[C0.0.h.5$r  ==  0.1]), precision=3)
        pEig    <-  rounded(length(C0.0.h.5$sf[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.1]) / 
                            length(C0.0.h.5$sf[C0.0.h.5$r  ==  0.1]), precision=3)

        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.0.h.5$sf[C0.0.h.5$Poly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.1] ~
               C0.0.h.5$sm[C0.0.h.5$Poly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.0.h.5$sf[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.1] ~
               C0.0.h.5$sm[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[5])
        points(C0.0.h.5$sf[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 1 & C0.0.h.5$r  ==  0.1] ~
               C0.0.h.5$sm[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 1 & C0.0.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.1.Hi[r0.1.Hi > twoLoc.Hi.obOut & r0.1.Hi <= 1] ~ sm[r0.1.Hi > twoLoc.Hi.obOut & r0.1.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.1.Lo[r0.1.Lo < twoLoc.Lo.obOut] ~ sm[r0.1.Lo < twoLoc.Lo.obOut], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi.obOut[twoLoc.Hi.obOut < 1] ~ sm[twoLoc.Hi.obOut < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo.obOut ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = 0.1")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.05, 1.075, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.5, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi.obOut)
        rm(twoLoc.Lo.obOut)
        rm(r0.1.Hi)
        rm(r0.1.Lo)


    ##  Panel 6: C = 0.25
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.add(hf=0.5, hm=0.5, sm=sm, C=0.25)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.add(hf=0.5, hm=0.5, sm=sm, C=0.25)
        
        r0.1.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.1, C=0.25)
        r0.1.Hi[r0.1.Hi > 1]  <-  1.00000001
        r0.1.Hi[r0.1.Hi == 'NaN']  <-  1.00000001
        r0.1.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.1, C=0.25)
        
        pAgree  <-  rounded(sum(C0.25.h.5$agree[C0.25.h.5$r == 0.1])/length(C0.25.h.5$agree[C0.25.h.5$r == 0.1]), precision=3)
        pSim    <-  rounded(length(C0.25.h.5$sf[C0.25.h.5$Poly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.1]) / 
                            length(C0.25.h.5$sf[C0.25.h.5$r  ==  0.1]), precision=3)
        pEig    <-  rounded(length(C0.25.h.5$sf[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.1]) / 
                            length(C0.25.h.5$sf[C0.25.h.5$r  ==  0.1]), precision=3)
        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.25.h.5$sf[C0.25.h.5$Poly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.1] ~
               C0.25.h.5$sm[C0.25.h.5$Poly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.25.h.5$sf[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.1] ~
               C0.25.h.5$sm[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.25.h.5$sf[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 1 & C0.25.h.5$r  ==  0.1] ~
               C0.25.h.5$sm[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 1 & C0.25.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.1.Hi[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1] ~ sm[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.1.Lo[r0.1.Lo < twoLoc.Lo] ~ sm[r0.1.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi[twoLoc.Hi < 1] ~ sm[twoLoc.Hi < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.1.Hi)
        rm(r0.1.Lo)



    ##  Panel 7: C = 0.5
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.add(hf=0.5, hm=0.5, sm=sm, C=0.5)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.add(hf=0.5, hm=0.5, sm=sm, C=0.5)
        
        r0.1.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.1, C=0.5)
        r0.1.Hi[r0.1.Hi > 1]  <-  1.00000001
        r0.1.Hi[r0.1.Hi == 'NaN']  <-  1.00000001
        r0.1.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.1, C=0.5)
        
        pAgree  <-  rounded(sum(C0.5.h.5$agree[C0.5.h.5$r == 0.1])/length(C0.5.h.5$agree[C0.5.h.5$r == 0.1]), precision=3)
        pSim    <-  rounded(length(C0.5.h.5$sf[C0.5.h.5$Poly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.1]) / 
                            length(C0.5.h.5$sf[C0.5.h.5$r  ==  0.1]), precision=3)
        pEig    <-  rounded(length(C0.5.h.5$sf[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.1]) / 
                            length(C0.5.h.5$sf[C0.5.h.5$r  ==  0.1]), precision=3)
        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.5.h.5$sf[C0.5.h.5$Poly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.1] ~
               C0.5.h.5$sm[C0.5.h.5$Poly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.5.h.5$sf[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.1] ~
               C0.5.h.5$sm[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.5.h.5$sf[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 1 & C0.5.h.5$r  ==  0.1] ~
               C0.5.h.5$sm[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 1 & C0.5.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.1.Hi[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1] ~ sm[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.1.Lo[r0.1.Lo < twoLoc.Lo] ~ sm[r0.1.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi[twoLoc.Hi < 1] ~ sm[twoLoc.Hi < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.1.Hi)
        rm(r0.1.Lo)


    ##  Panel 8: C = 0.75
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.add(hf=0.5, hm=0.5, sm=sm, C=0.75)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.add(hf=0.5, hm=0.5, sm=sm, C=0.75)
        
        r0.1.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.1, C=0.75)
        r0.1.Hi[r0.1.Hi > 1]  <-  1.00000001
        r0.1.Hi[r0.1.Hi == 'NaN']  <-  1.00000001
        r0.1.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.1, C=0.75)
        
        pAgree  <-  rounded(sum(C0.75.h.5$agree[C0.75.h.5$r == 0.1])/length(C0.75.h.5$agree[C0.75.h.5$r == 0.1]), precision=3)
        pSim    <-  rounded(length(C0.75.h.5$sf[C0.75.h.5$Poly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.1]) / 
                            length(C0.75.h.5$sf[C0.75.h.5$r  ==  0.1]), precision=3)
        pEig    <-  rounded(length(C0.75.h.5$sf[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.1]) / 
                            length(C0.75.h.5$sf[C0.75.h.5$r  ==  0.1]), precision=3)

        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.75.h.5$sf[C0.75.h.5$Poly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.1] ~
               C0.75.h.5$sm[C0.75.h.5$Poly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.75.h.5$sf[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.1] ~
               C0.75.h.5$sm[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.75.h.5$sf[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 1 & C0.75.h.5$r  ==  0.1] ~
               C0.75.h.5$sm[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 1 & C0.75.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.1.Hi[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1] ~ sm[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.1.Lo[r0.1.Lo < twoLoc.Lo] ~ sm[r0.1.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi[twoLoc.Hi < 1] ~ sm[twoLoc.Hi < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'H', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.1.Hi)
        rm(r0.1.Lo)






##  Row 3: r = 0.2
    ##  Panel 9: C = 0
        
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi.obOut   <-  inv.lab1.add.obOut(hf=0.5, hm=0.5, sm=sm)
        twoLoc.Hi.obOut[twoLoc.Hi.obOut > 1]  <-  1.00000001
        twoLoc.Hi.obOut[10001]  <-  1.00000001
        twoLoc.Lo.obOut  <-  inv.lAB1.add.obOut(hf=0.5, hm=0.5, sm=sm)
        
        r0.2.Hi  <-  inv.lab2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0.2)
        r0.2.Hi[r0.2.Hi > 1]  <-  1.00000001
        r0.2.Hi[r0.2.Hi == 'NaN']  <-  1.00000001
        r0.2.Lo  <-  inv.lAB2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0.2)

        pAgree  <-  rounded(sum(C0.0.h.5$agree[C0.0.h.5$r == 0.2])/length(C0.0.h.5$agree[C0.0.h.5$r == 0.2]), precision=3)
        pSim    <-  rounded(length(C0.0.h.5$sf[C0.0.h.5$Poly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.2]) / 
                            length(C0.0.h.5$sf[C0.0.h.5$r  ==  0.2]), precision=3)
        pEig    <-  rounded(length(C0.0.h.5$sf[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.2]) / 
                            length(C0.0.h.5$sf[C0.0.h.5$r  ==  0.2]), precision=3)

        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.0.h.5$sf[C0.0.h.5$Poly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.1] ~
               C0.0.h.5$sm[C0.0.h.5$Poly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.0.h.5$sf[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.1] ~
               C0.0.h.5$sm[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[5])
        points(C0.0.h.5$sf[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 1 & C0.0.h.5$r  ==  0.1] ~
               C0.0.h.5$sm[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 1 & C0.0.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.2.Hi[r0.2.Hi > twoLoc.Hi.obOut & r0.2.Hi <= 1] ~ sm[r0.2.Hi > twoLoc.Hi.obOut & r0.2.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.2.Lo[r0.2.Lo < twoLoc.Lo.obOut] ~ sm[r0.2.Lo < twoLoc.Lo.obOut], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi.obOut[twoLoc.Hi.obOut < 1] ~ sm[twoLoc.Hi.obOut < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo.obOut ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = 0.2")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.05, 1.075, 'I', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi.obOut)
        rm(twoLoc.Lo.obOut)
        rm(r0.2.Hi)
        rm(r0.2.Lo)


    ##  Panel 10: C = 0.25
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.add(hf=0.5, hm=0.5, sm=sm, C=0.25)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.add(hf=0.5, hm=0.5, sm=sm, C=0.25)
        
        r0.2.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.2, C=0.25)
        r0.2.Hi[r0.2.Hi > 1]  <-  1.00000001
        r0.2.Hi[r0.2.Hi == 'NaN']  <-  1.00000001
        r0.2.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.2, C=0.25)
        
        pAgree  <-  rounded(sum(C0.25.h.5$agree[C0.25.h.5$r == 0.2])/length(C0.25.h.5$agree[C0.25.h.5$r == 0.2]), precision=3)
        pSim    <-  rounded(length(C0.25.h.5$sf[C0.25.h.5$Poly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.2]) / 
                            length(C0.25.h.5$sf[C0.25.h.5$r  ==  0.2]), precision=3)
        pEig    <-  rounded(length(C0.25.h.5$sf[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.2]) / 
                            length(C0.25.h.5$sf[C0.25.h.5$r  ==  0.2]), precision=3)

        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.25.h.5$sf[C0.25.h.5$Poly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.1] ~
               C0.25.h.5$sm[C0.25.h.5$Poly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.25.h.5$sf[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.1] ~
               C0.25.h.5$sm[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.25.h.5$sf[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 1 & C0.25.h.5$r  ==  0.1] ~
               C0.25.h.5$sm[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 1 & C0.25.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.2.Hi[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1] ~ sm[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.2.Lo[r0.2.Lo < twoLoc.Lo] ~ sm[r0.2.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi[twoLoc.Hi < 1] ~ sm[twoLoc.Hi < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'J', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.2.Hi)
        rm(r0.2.Lo)



    ##  Panel 11: C = 0.5
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.add(hf=0.5, hm=0.5, sm=sm, C=0.5)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.add(hf=0.5, hm=0.5, sm=sm, C=0.5)
        
        r0.2.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.2, C=0.5)
        r0.2.Hi[r0.2.Hi > 1]  <-  1.00000001
        r0.2.Hi[r0.2.Hi == 'NaN']  <-  1.00000001
        r0.2.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.2, C=0.5)
        
        pAgree  <-  rounded(sum(C0.5.h.5$agree[C0.5.h.5$r == 0.2])/length(C0.5.h.5$agree[C0.5.h.5$r == 0.2]), precision=3)
        pSim    <-  rounded(length(C0.5.h.5$sf[C0.5.h.5$Poly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.2]) / 
                            length(C0.5.h.5$sf[C0.5.h.5$r  ==  0.2]), precision=3)
        pEig    <-  rounded(length(C0.5.h.5$sf[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.2]) / 
                            length(C0.5.h.5$sf[C0.5.h.5$r  ==  0.2]), precision=3)

        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.5.h.5$sf[C0.5.h.5$Poly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.1] ~
               C0.5.h.5$sm[C0.5.h.5$Poly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.5.h.5$sf[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.1] ~
               C0.5.h.5$sm[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.5.h.5$sf[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 1 & C0.5.h.5$r  ==  0.1] ~
               C0.5.h.5$sm[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 1 & C0.5.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.2.Hi[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1] ~ sm[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.2.Lo[r0.2.Lo < twoLoc.Lo] ~ sm[r0.2.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi[twoLoc.Hi < 1] ~ sm[twoLoc.Hi < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'K', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.2.Hi)
        rm(r0.2.Lo)


    ##  Panel 12: C = 0.75
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.add(hf=0.5, hm=0.5, sm=sm, C=0.75)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.add(hf=0.5, hm=0.5, sm=sm, C=0.75)
        
        r0.2.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.2, C=0.75)
        r0.2.Hi[r0.2.Hi > 1]  <-  1.00000001
        r0.2.Hi[r0.2.Hi == 'NaN']  <-  1.00000001
        r0.2.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.2, C=0.75)
        
        pAgree  <-  rounded(sum(C0.75.h.5$agree[C0.75.h.5$r == 0.2])/length(C0.75.h.5$agree[C0.75.h.5$r == 0.2]), precision=3)
        pSim    <-  rounded(length(C0.75.h.5$sf[C0.75.h.5$Poly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.2]) / 
                            length(C0.75.h.5$sf[C0.75.h.5$r  ==  0.2]), precision=3)
        pEig    <-  rounded(length(C0.75.h.5$sf[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.2]) / 
                            length(C0.75.h.5$sf[C0.75.h.5$r  ==  0.2]), precision=3)

        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.75.h.5$sf[C0.75.h.5$Poly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.1] ~
               C0.75.h.5$sm[C0.75.h.5$Poly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.75.h.5$sf[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.1] ~
               C0.75.h.5$sm[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.75.h.5$sf[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 1 & C0.75.h.5$r  ==  0.1] ~
               C0.75.h.5$sm[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 1 & C0.75.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.2.Hi[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1] ~ sm[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.2.Lo[r0.2.Lo < twoLoc.Lo] ~ sm[r0.2.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi[twoLoc.Hi < 1] ~ sm[twoLoc.Hi < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'L', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.2.Hi)
        rm(r0.2.Lo)






##  Row 4: r = 0.2
    ##  Panel 13: C = 0
        
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi.obOut   <-  inv.lab1.add.obOut(hf=0.5, hm=0.5, sm=sm)
        twoLoc.Hi.obOut[twoLoc.Hi.obOut > 1]  <-  1.00000001
        twoLoc.Hi.obOut[10001]  <-  1.00000001
        twoLoc.Lo.obOut  <-  inv.lAB1.add.obOut(hf=0.5, hm=0.5, sm=sm)
        
        r0.5.Hi  <-  inv.lab2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0.5)
        r0.5.Hi[r0.5.Hi > 1]  <-  1.00000001
        r0.5.Hi[r0.5.Hi == 'NaN']  <-  1.00000001
        r0.5.Lo  <-  inv.lAB2.add.obOut(hf = 0.5, hm = 0.5, sm=sm, r=0.5)

        pAgree  <-  rounded(sum(C0.0.h.5$agree[C0.0.h.5$r == 0.2])/length(C0.0.h.5$agree[C0.0.h.5$r == 0.5]), precision=3)
        pSim    <-  rounded(length(C0.0.h.5$sf[C0.0.h.5$Poly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.5]) / 
                            length(C0.0.h.5$sf[C0.0.h.5$r  ==  0.5]), precision=3)
        pEig    <-  rounded(length(C0.0.h.5$sf[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.5]) / 
                            length(C0.0.h.5$sf[C0.0.h.5$r  ==  0.5]), precision=3)

        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.0.h.5$sf[C0.0.h.5$Poly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.1] ~
               C0.0.h.5$sm[C0.0.h.5$Poly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.0.h.5$sf[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.1] ~
               C0.0.h.5$sm[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 0 & C0.0.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[5])
        points(C0.0.h.5$sf[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 1 & C0.0.h.5$r  ==  0.1] ~
               C0.0.h.5$sm[C0.0.h.5$eigPoly == 1 & C0.0.h.5$agree == 1 & C0.0.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.5.Hi[r0.5.Hi > twoLoc.Hi.obOut & r0.5.Hi <= 1] ~ sm[r0.5.Hi > twoLoc.Hi.obOut & r0.5.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.5.Lo[r0.5.Lo < twoLoc.Lo.obOut] ~ sm[r0.5.Lo < twoLoc.Lo.obOut], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi.obOut[twoLoc.Hi.obOut < 1] ~ sm[twoLoc.Hi.obOut < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo.obOut ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(0.05, 1.075, 'M', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = 0.5")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.4, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi.obOut)
        rm(twoLoc.Lo.obOut)
        rm(r0.5.Hi)
        rm(r0.5.Lo)


    ##  Panel 14: C = 0.25
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.add(hf=0.5, hm=0.5, sm=sm, C=0.25)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.add(hf=0.5, hm=0.5, sm=sm, C=0.25)
        
        r0.5.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.5, C=0.25)
        r0.5.Hi[r0.5.Hi > 1]  <-  1.00000001
        r0.5.Hi[r0.5.Hi == 'NaN']  <-  1.00000001
        r0.5.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.5, C=0.25)
        
        pAgree  <-  rounded(sum(C0.25.h.5$agree[C0.25.h.5$r == 0.2])/length(C0.25.h.5$agree[C0.25.h.5$r == 0.5]), precision=3)
        pSim    <-  rounded(length(C0.25.h.5$sf[C0.25.h.5$Poly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.5]) / 
                            length(C0.25.h.5$sf[C0.25.h.5$r  ==  0.5]), precision=3)
        pEig    <-  rounded(length(C0.25.h.5$sf[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.5]) / 
                            length(C0.25.h.5$sf[C0.25.h.5$r  ==  0.5]), precision=3)

        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.25.h.5$sf[C0.25.h.5$Poly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.1] ~
               C0.25.h.5$sm[C0.25.h.5$Poly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.25.h.5$sf[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.1] ~
               C0.25.h.5$sm[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 0 & C0.25.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.25.h.5$sf[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 1 & C0.25.h.5$r  ==  0.1] ~
               C0.25.h.5$sm[C0.25.h.5$eigPoly == 1 & C0.25.h.5$agree == 1 & C0.25.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.5.Hi[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1] ~ sm[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.5.Lo[r0.5.Lo < twoLoc.Lo] ~ sm[r0.5.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi[twoLoc.Hi < 1] ~ sm[twoLoc.Hi < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'N', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.4, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.5.Hi)
        rm(r0.5.Lo)



    ##  Panel 15: C = 0.5
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.add(hf=0.5, hm=0.5, sm=sm, C=0.5)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.add(hf=0.5, hm=0.5, sm=sm, C=0.5)
        
        r0.5.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.5, C=0.5)
        r0.5.Hi[r0.5.Hi > 1]  <-  1.00000001
        r0.5.Hi[r0.5.Hi == 'NaN']  <-  1.00000001
        r0.5.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.5, C=0.5)
        
        pAgree  <-  rounded(sum(C0.5.h.5$agree[C0.5.h.5$r == 0.2])/length(C0.5.h.5$agree[C0.5.h.5$r == 0.5]), precision=3)
        pSim    <-  rounded(length(C0.5.h.5$sf[C0.5.h.5$Poly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.5]) / 
                            length(C0.5.h.5$sf[C0.5.h.5$r  ==  0.5]), precision=3)
        pEig    <-  rounded(length(C0.5.h.5$sf[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.5]) / 
                            length(C0.5.h.5$sf[C0.5.h.5$r  ==  0.5]), precision=3)

        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.5.h.5$sf[C0.5.h.5$Poly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.1] ~
               C0.5.h.5$sm[C0.5.h.5$Poly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.5.h.5$sf[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.1] ~
               C0.5.h.5$sm[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 0 & C0.5.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.5.h.5$sf[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 1 & C0.5.h.5$r  ==  0.1] ~
               C0.5.h.5$sm[C0.5.h.5$eigPoly == 1 & C0.5.h.5$agree == 1 & C0.5.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.5.Hi[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1] ~ sm[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.5.Lo[r0.5.Lo < twoLoc.Lo] ~ sm[r0.5.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi[twoLoc.Hi < 1] ~ sm[twoLoc.Hi < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'O', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.4, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.5.Hi)
        rm(r0.5.Lo)


    ##  Panel 16: C = 0.75
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.add(hf=0.5, hm=0.5, sm=sm, C=0.75)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.add(hf=0.5, hm=0.5, sm=sm, C=0.75)
        
        r0.5.Hi  <-  inv.lab2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.5, C=0.75)
        r0.5.Hi[r0.5.Hi > 1]  <-  1.00000001
        r0.5.Hi[r0.5.Hi == 'NaN']  <-  1.00000001
        r0.5.Lo  <-  inv.lAB2.add(hf = 0.5, hm = 0.5, sm=sm, r=0.5, C=0.75)
        
        pAgree  <-  rounded(sum(C0.75.h.5$agree[C0.75.h.5$r == 0.2])/length(C0.75.h.5$agree[C0.75.h.5$r == 0.5]), precision=3)
        pSim    <-  rounded(length(C0.75.h.5$sf[C0.75.h.5$Poly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.5]) / 
                            length(C0.75.h.5$sf[C0.75.h.5$r  ==  0.5]), precision=3)
        pEig    <-  rounded(length(C0.75.h.5$sf[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.5]) / 
                            length(C0.75.h.5$sf[C0.75.h.5$r  ==  0.5]), precision=3)

        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.75.h.5$sf[C0.75.h.5$Poly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.1] ~
               C0.75.h.5$sm[C0.75.h.5$Poly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.75.h.5$sf[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.1] ~
               C0.75.h.5$sm[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 0 & C0.75.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.75.h.5$sf[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 1 & C0.75.h.5$r  ==  0.1] ~
               C0.75.h.5$sm[C0.75.h.5$eigPoly == 1 & C0.75.h.5$agree == 1 & C0.75.h.5$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.5.Hi[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1] ~ sm[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.5.Lo[r0.5.Lo < twoLoc.Lo] ~ sm[r0.5.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi[twoLoc.Hi < 1] ~ sm[twoLoc.Hi < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'P', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.4, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.5.Hi)
        rm(r0.5.Lo)


}










































#' Fig.S1: Supplementary figure showing comparison between deterministic
#'         recursion simulations and invasion analysis based on eigenvalues
#' 
#'
#' @title Fig.S1: Supplementary figure showing comparison between deterministic
#'                recursion simulations and invasion analysis based on eigenvalues
#' @export

Fig.S2_domRev  <-  function() {

    ## Read data files for plotting
        C0.0.h.25   <-  read.table('./output/data/simResults/recFwdSimLoop.out_C0_hf0.25_hm0.25_sMax1.txt', head=TRUE)
        C0.25.h.25  <-  read.table('./output/data/simResults/recFwdSimLoop.out_C0.25_hf0.25_hm0.25_sMax1.txt', head=TRUE)
        C0.5.h.25   <-  read.table('./output/data/simResults/recFwdSimLoop.out_C0.5_hf0.25_hm0.25_sMax1.txt', head=TRUE)
        C0.75.h.25  <-  read.table('./output/data/simResults/recFwdSimLoop.out_C0.75_hf0.25_hm0.25_sMax1.txt', head=TRUE)

    # Color scheme
    COLS  <-  c(transparentColor('seagreen3', opacity=0.2), 
                transparentColor('dodgerblue2', opacity=0.2), 
                transparentColor('tomato2', opacity=0.2),
                'black')

    #  Create vector of male selection coefficients for invasion functions
    sm  <-  seq(0,1,by=0.0001)

    # Set plot layout
    layout.mat <- matrix(c(1:16), nrow=4, ncol=4, byrow=TRUE)
    layout <- layout(layout.mat,respect=TRUE)

##  Row 1: r = 0
    ##  Panel One: C = 0
        
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi.obOut   <-  inv.lab1.domRev.obOut(hf=0.25, hm=0.25, sm=sm)
        twoLoc.Hi.obOut[twoLoc.Hi.obOut > 1]  <-  1.00000001
        twoLoc.Hi.obOut[10001]  <-  1.00000001
        twoLoc.Lo.obOut  <-  inv.lAB1.domRev.obOut(hf=0.25, hm=0.25, sm=sm)
        
        r0.Hi  <-  inv.lab2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0)
        r0.Hi[r0.Hi > 1]  <-  1.00000001
        r0.Hi[r0.Hi == 'NaN']  <-  1.00000001
        r0.Lo  <-  inv.lAB2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0)

        pAgree  <-  rounded(sum(C0.0.h.25$agree[C0.0.h.25$r == 0.0])/length(C0.0.h.25$agree[C0.0.h.25$r == 0.0]), precision=3)
        pSim    <-  rounded(length(C0.0.h.25$sf[C0.0.h.25$Poly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.0]) / 
                            length(C0.0.h.25$sf[C0.0.h.25$r  ==  0.0]), precision=3)
        pEig    <-  rounded(length(C0.0.h.25$sf[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.0]) / 
                            length(C0.0.h.25$sf[C0.0.h.25$r  ==  0.0]), precision=3)

        # Make the plot
        par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.0.h.25$sf[C0.0.h.25$Poly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.0] ~
               C0.0.h.25$sm[C0.0.h.25$Poly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.0.h.25$sf[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.0] ~
               C0.0.h.25$sm[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[5])
        points(C0.0.h.25$sf[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 1 & C0.0.h.25$r  ==  0.0] ~
               C0.0.h.25$sm[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 1 & C0.0.h.25$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.Hi[r0.Hi > twoLoc.Hi.obOut & r0.Hi <= 1] ~ sm[r0.Hi > twoLoc.Hi.obOut & r0.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.Lo[r0.Lo < twoLoc.Lo.obOut] ~ sm[r0.Lo < twoLoc.Lo.obOut], lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(0.5, 1.25, expression(paste(italic(C), ' = ', 0)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = 0")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.05, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.68, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.68, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.68, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.5, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi.obOut)
        rm(twoLoc.Lo.obOut)
        rm(r0.Hi)
        rm(r0.Lo)


    ##  Panel Two: C = 0.25
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.25)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.25)
        
        r0.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0, C=0.25)
        r0.Hi[r0.Hi > 1]  <-  1.00000001
        r0.Hi[r0.Hi == 'NaN']  <-  1.00000001
        r0.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0, C=0.25)

        pAgree  <-  rounded(sum(C0.25.h.25$agree[C0.25.h.25$r == 0.0])/length(C0.25.h.25$agree[C0.25.h.25$r == 0.0]), precision=3)
        pSim    <-  rounded(length(C0.25.h.25$sf[C0.25.h.25$Poly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.0]) / 
                            length(C0.25.h.25$sf[C0.25.h.25$r  ==  0.0]), precision=3)
        pEig    <-  rounded(length(C0.25.h.25$sf[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.0]) / 
                            length(C0.25.h.25$sf[C0.25.h.25$r  ==  0.0]), precision=3)
        
        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.25.h.25$sf[C0.25.h.25$Poly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.0] ~
               C0.25.h.25$sm[C0.25.h.25$Poly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.25.h.25$sf[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.0] ~
               C0.25.h.25$sm[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.25.h.25$sf[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 1 & C0.25.h.25$r  ==  0.0] ~
               C0.25.h.25$sm[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 1 & C0.25.h.25$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.Hi[r0.Hi > twoLoc.Hi & r0.Hi <= 1] ~ sm[r0.Hi > twoLoc.Hi & r0.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.Lo[r0.Lo < twoLoc.Lo] ~ sm[r0.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.25, expression(paste(italic(C), ' = ',0.25)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.5, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.Hi)
        rm(r0.Lo)



    ##  Panel Three: C = 0.5
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.5)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.5)
        
        r0.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0, C=0.5)
        r0.Hi[r0.Hi > 1]  <-  1.00000001
        r0.Hi[r0.Hi == 'NaN']  <-  1.00000001
        r0.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0, C=0.5)

        pAgree  <-  rounded(sum(C0.5.h.25$agree[C0.5.h.25$r == 0.0])/length(C0.5.h.25$agree[C0.5.h.25$r == 0.0]), precision=3)
        pSim    <-  rounded(length(C0.5.h.25$sf[C0.5.h.25$Poly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.0]) / 
                            length(C0.5.h.25$sf[C0.5.h.25$r  ==  0.0]), precision=3)
        pEig    <-  rounded(length(C0.5.h.25$sf[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.0]) / 
                            length(C0.5.h.25$sf[C0.5.h.25$r  ==  0.0]), precision=3)
        
        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.5.h.25$sf[C0.5.h.25$Poly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.0] ~
               C0.5.h.25$sm[C0.5.h.25$Poly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.5.h.25$sf[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.0] ~
               C0.5.h.25$sm[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.5.h.25$sf[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 1 & C0.5.h.25$r  ==  0.0] ~
               C0.5.h.25$sm[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 1 & C0.5.h.25$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.Hi[r0.Hi > twoLoc.Hi & r0.Hi <= 1] ~ sm[r0.Hi > twoLoc.Hi & r0.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.Lo[r0.Lo < twoLoc.Lo] ~ sm[r0.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.25, expression(paste(italic(C), ' = ',0.5)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.5, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.Hi)
        rm(r0.Lo)


    ##  Panel Four: C = 0.75
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.75)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.75)
        
        r0.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0, C=0.75)
        r0.Hi[r0.Hi > 1]  <-  1.00000001
        r0.Hi[r0.Hi == 'NaN']  <-  1.00000001
        r0.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0, C=0.75)

        pAgree  <-  rounded(sum(C0.75.h.25$agree[C0.75.h.25$r == 0.0])/length(C0.75.h.25$agree[C0.75.h.25$r == 0.0]), precision=3)
        pSim    <-  rounded(length(C0.75.h.25$sf[C0.75.h.25$Poly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.0]) / 
                            length(C0.75.h.25$sf[C0.75.h.25$r  ==  0.0]), precision=3)
        pEig    <-  rounded(length(C0.75.h.25$sf[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.0]) / 
                            length(C0.75.h.25$sf[C0.75.h.25$r  ==  0.0]), precision=3)
        
        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.75.h.25$sf[C0.75.h.25$Poly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.0] ~
               C0.75.h.25$sm[C0.75.h.25$Poly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.75.h.25$sf[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.0] ~
               C0.75.h.25$sm[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.75.h.25$sf[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 1 & C0.75.h.25$r  ==  0.0] ~
               C0.75.h.25$sm[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 1 & C0.75.h.25$r  ==  0.0], pch=21, col=NA, cex=0.75, bg=COLS[1])
        points(0.02,0.99, pch=21, col=NA, cex=0.8, bg='seagreen3')
        points(0.02,0.91, pch=21, col=NA, cex=0.8, bg='tomato2')
        points(0.02,0.83, pch=21, col=NA, cex=0.8, bg='dodgerblue2')
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.Hi[r0.Hi > twoLoc.Hi & r0.Hi <= 1] ~ sm[r0.Hi > twoLoc.Hi & r0.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.Lo[r0.Lo < twoLoc.Lo] ~ sm[r0.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.5, 1.25, expression(paste(italic(C), ' = ',0.75)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 1.075, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.6, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.88, substitute(p~" Sim.", list(p = pSim)),    cex=0.6, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.1, 0.80, substitute(p~" Eig.", list(p = pEig)),    cex=0.6, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.Hi)
        rm(r0.Lo)



##  Row 2: r = 0.1
    ##  Panel 5: C = 0
        
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi.obOut   <-  inv.lab1.domRev.obOut(hf=0.25, hm=0.25, sm=sm)
        twoLoc.Hi.obOut[twoLoc.Hi.obOut > 1]  <-  1.00000001
        twoLoc.Hi.obOut[10001]  <-  1.00000001
        twoLoc.Lo.obOut  <-  inv.lAB1.domRev.obOut(hf=0.25, hm=0.25, sm=sm)
        
        r0.1.Hi  <-  inv.lab2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0.1)
        r0.1.Hi[r0.1.Hi > 1]  <-  1.00000001
        r0.1.Hi[r0.1.Hi == 'NaN']  <-  1.00000001
        r0.1.Lo  <-  inv.lAB2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0.1)

        pAgree  <-  rounded(sum(C0.0.h.25$agree[C0.0.h.25$r == 0.1])/length(C0.0.h.25$agree[C0.0.h.25$r == 0.1]), precision=3)
        pSim    <-  rounded(length(C0.0.h.25$sf[C0.0.h.25$Poly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.1]) / 
                            length(C0.0.h.25$sf[C0.0.h.25$r  ==  0.1]), precision=3)
        pEig    <-  rounded(length(C0.0.h.25$sf[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.1]) / 
                            length(C0.0.h.25$sf[C0.0.h.25$r  ==  0.1]), precision=3)

        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.0.h.25$sf[C0.0.h.25$Poly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.1] ~
               C0.0.h.25$sm[C0.0.h.25$Poly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.0.h.25$sf[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.1] ~
               C0.0.h.25$sm[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[5])
        points(C0.0.h.25$sf[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 1 & C0.0.h.25$r  ==  0.1] ~
               C0.0.h.25$sm[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 1 & C0.0.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.1.Hi[r0.1.Hi > twoLoc.Hi.obOut & r0.1.Hi <= 1] ~ sm[r0.1.Hi > twoLoc.Hi.obOut & r0.1.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.1.Lo[r0.1.Lo < twoLoc.Lo.obOut] ~ sm[r0.1.Lo < twoLoc.Lo.obOut], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi.obOut[twoLoc.Hi.obOut < 1] ~ sm[twoLoc.Hi.obOut < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo.obOut ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = 0.1")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.05, 1.075, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.68, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.68, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.5, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.68, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.5, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi.obOut)
        rm(twoLoc.Lo.obOut)
        rm(r0.1.Hi)
        rm(r0.1.Lo)


    ##  Panel 6: C = 0.25
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.25)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.25)
        
        r0.1.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.1, C=0.25)
        r0.1.Hi[r0.1.Hi > 1]  <-  1.00000001
        r0.1.Hi[r0.1.Hi == 'NaN']  <-  1.00000001
        r0.1.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.1, C=0.25)
        
        pAgree  <-  rounded(sum(C0.25.h.25$agree[C0.25.h.25$r == 0.1])/length(C0.25.h.25$agree[C0.25.h.25$r == 0.1]), precision=3)
        pSim    <-  rounded(length(C0.25.h.25$sf[C0.25.h.25$Poly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.1]) / 
                            length(C0.25.h.25$sf[C0.25.h.25$r  ==  0.1]), precision=3)
        pEig    <-  rounded(length(C0.25.h.25$sf[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.1]) / 
                            length(C0.25.h.25$sf[C0.25.h.25$r  ==  0.1]), precision=3)
        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.25.h.25$sf[C0.25.h.25$Poly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.1] ~
               C0.25.h.25$sm[C0.25.h.25$Poly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.25.h.25$sf[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.1] ~
               C0.25.h.25$sm[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.25.h.25$sf[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 1 & C0.25.h.25$r  ==  0.1] ~
               C0.25.h.25$sm[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 1 & C0.25.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.1.Hi[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1] ~ sm[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.1.Lo[r0.1.Lo < twoLoc.Lo] ~ sm[r0.1.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi[twoLoc.Hi < 1] ~ sm[twoLoc.Hi < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.1.Hi)
        rm(r0.1.Lo)



    ##  Panel 7: C = 0.5
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.5)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.5)
        
        r0.1.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.1, C=0.5)
        r0.1.Hi[r0.1.Hi > 1]  <-  1.00000001
        r0.1.Hi[r0.1.Hi == 'NaN']  <-  1.00000001
        r0.1.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.1, C=0.5)
        
        pAgree  <-  rounded(sum(C0.5.h.25$agree[C0.5.h.25$r == 0.1])/length(C0.5.h.25$agree[C0.5.h.25$r == 0.1]), precision=3)
        pSim    <-  rounded(length(C0.5.h.25$sf[C0.5.h.25$Poly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.1]) / 
                            length(C0.5.h.25$sf[C0.5.h.25$r  ==  0.1]), precision=3)
        pEig    <-  rounded(length(C0.5.h.25$sf[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.1]) / 
                            length(C0.5.h.25$sf[C0.5.h.25$r  ==  0.1]), precision=3)
        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.5.h.25$sf[C0.5.h.25$Poly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.1] ~
               C0.5.h.25$sm[C0.5.h.25$Poly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.5.h.25$sf[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.1] ~
               C0.5.h.25$sm[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.5.h.25$sf[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 1 & C0.5.h.25$r  ==  0.1] ~
               C0.5.h.25$sm[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 1 & C0.5.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.1.Hi[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1] ~ sm[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.1.Lo[r0.1.Lo < twoLoc.Lo] ~ sm[r0.1.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi[twoLoc.Hi < 1] ~ sm[twoLoc.Hi < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.1.Hi)
        rm(r0.1.Lo)


    ##  Panel 8: C = 0.75
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.75)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.75)
        
        r0.1.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.1, C=0.75)
        r0.1.Hi[r0.1.Hi > 1]  <-  1.00000001
        r0.1.Hi[r0.1.Hi == 'NaN']  <-  1.00000001
        r0.1.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.1, C=0.75)
        
        pAgree  <-  rounded(sum(C0.75.h.25$agree[C0.75.h.25$r == 0.1])/length(C0.75.h.25$agree[C0.75.h.25$r == 0.1]), precision=3)
        pSim    <-  rounded(length(C0.75.h.25$sf[C0.75.h.25$Poly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.1]) / 
                            length(C0.75.h.25$sf[C0.75.h.25$r  ==  0.1]), precision=3)
        pEig    <-  rounded(length(C0.75.h.25$sf[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.1]) / 
                            length(C0.75.h.25$sf[C0.75.h.25$r  ==  0.1]), precision=3)

        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.75.h.25$sf[C0.75.h.25$Poly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.1] ~
               C0.75.h.25$sm[C0.75.h.25$Poly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.75.h.25$sf[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.1] ~
               C0.75.h.25$sm[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.75.h.25$sf[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 1 & C0.75.h.25$r  ==  0.1] ~
               C0.75.h.25$sm[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 1 & C0.75.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.1.Hi[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1] ~ sm[r0.1.Hi > twoLoc.Hi & r0.1.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.1.Lo[r0.1.Lo < twoLoc.Lo] ~ sm[r0.1.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi[twoLoc.Hi < 1] ~ sm[twoLoc.Hi < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'H', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.1.Hi)
        rm(r0.1.Lo)






##  Row 3: r = 0.2
    ##  Panel 9: C = 0
        
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi.obOut   <-  inv.lab1.domRev.obOut(hf=0.25, hm=0.25, sm=sm)
        twoLoc.Hi.obOut[twoLoc.Hi.obOut > 1]  <-  1.00000001
        twoLoc.Hi.obOut[10001]  <-  1.00000001
        twoLoc.Lo.obOut  <-  inv.lAB1.domRev.obOut(hf=0.25, hm=0.25, sm=sm)
        
        r0.2.Hi  <-  inv.lab2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0.2)
        r0.2.Hi[r0.2.Hi > 1]  <-  1.00000001
        r0.2.Hi[r0.2.Hi == 'NaN']  <-  1.00000001
        r0.2.Lo  <-  inv.lAB2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0.2)

        pAgree  <-  rounded(sum(C0.0.h.25$agree[C0.0.h.25$r == 0.2])/length(C0.0.h.25$agree[C0.0.h.25$r == 0.2]), precision=3)
        pSim    <-  rounded(length(C0.0.h.25$sf[C0.0.h.25$Poly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.2]) / 
                            length(C0.0.h.25$sf[C0.0.h.25$r  ==  0.2]), precision=3)
        pEig    <-  rounded(length(C0.0.h.25$sf[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.2]) / 
                            length(C0.0.h.25$sf[C0.0.h.25$r  ==  0.2]), precision=3)

        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.0.h.25$sf[C0.0.h.25$Poly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.1] ~
               C0.0.h.25$sm[C0.0.h.25$Poly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.0.h.25$sf[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.1] ~
               C0.0.h.25$sm[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[5])
        points(C0.0.h.25$sf[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 1 & C0.0.h.25$r  ==  0.1] ~
               C0.0.h.25$sm[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 1 & C0.0.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.2.Hi[r0.2.Hi > twoLoc.Hi.obOut & r0.2.Hi <= 1] ~ sm[r0.2.Hi > twoLoc.Hi.obOut & r0.2.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.2.Lo[r0.2.Lo < twoLoc.Lo.obOut] ~ sm[r0.2.Lo < twoLoc.Lo.obOut], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi.obOut[twoLoc.Hi.obOut < 1] ~ sm[twoLoc.Hi.obOut < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo.obOut ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = 0.2")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.05, 1.075, 'I', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.68, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.68, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.68, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi.obOut)
        rm(twoLoc.Lo.obOut)
        rm(r0.2.Hi)
        rm(r0.2.Lo)


    ##  Panel 10: C = 0.25
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.25)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.25)
        
        r0.2.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.2, C=0.25)
        r0.2.Hi[r0.2.Hi > 1]  <-  1.00000001
        r0.2.Hi[r0.2.Hi == 'NaN']  <-  1.00000001
        r0.2.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.2, C=0.25)
        
        pAgree  <-  rounded(sum(C0.25.h.25$agree[C0.25.h.25$r == 0.2])/length(C0.25.h.25$agree[C0.25.h.25$r == 0.2]), precision=3)
        pSim    <-  rounded(length(C0.25.h.25$sf[C0.25.h.25$Poly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.2]) / 
                            length(C0.25.h.25$sf[C0.25.h.25$r  ==  0.2]), precision=3)
        pEig    <-  rounded(length(C0.25.h.25$sf[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.2]) / 
                            length(C0.25.h.25$sf[C0.25.h.25$r  ==  0.2]), precision=3)

        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.25.h.25$sf[C0.25.h.25$Poly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.1] ~
               C0.25.h.25$sm[C0.25.h.25$Poly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.25.h.25$sf[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.1] ~
               C0.25.h.25$sm[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.25.h.25$sf[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 1 & C0.25.h.25$r  ==  0.1] ~
               C0.25.h.25$sm[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 1 & C0.25.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.2.Hi[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1] ~ sm[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.2.Lo[r0.2.Lo < twoLoc.Lo] ~ sm[r0.2.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi[twoLoc.Hi < 1] ~ sm[twoLoc.Hi < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'J', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.2.Hi)
        rm(r0.2.Lo)



    ##  Panel 11: C = 0.5
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.5)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.5)
        
        r0.2.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.2, C=0.5)
        r0.2.Hi[r0.2.Hi > 1]  <-  1.00000001
        r0.2.Hi[r0.2.Hi == 'NaN']  <-  1.00000001
        r0.2.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.2, C=0.5)
        
        pAgree  <-  rounded(sum(C0.5.h.25$agree[C0.5.h.25$r == 0.2])/length(C0.5.h.25$agree[C0.5.h.25$r == 0.2]), precision=3)
        pSim    <-  rounded(length(C0.5.h.25$sf[C0.5.h.25$Poly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.2]) / 
                            length(C0.5.h.25$sf[C0.5.h.25$r  ==  0.2]), precision=3)
        pEig    <-  rounded(length(C0.5.h.25$sf[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.2]) / 
                            length(C0.5.h.25$sf[C0.5.h.25$r  ==  0.2]), precision=3)

        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.5.h.25$sf[C0.5.h.25$Poly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.1] ~
               C0.5.h.25$sm[C0.5.h.25$Poly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.5.h.25$sf[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.1] ~
               C0.5.h.25$sm[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.5.h.25$sf[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 1 & C0.5.h.25$r  ==  0.1] ~
               C0.5.h.25$sm[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 1 & C0.5.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.2.Hi[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1] ~ sm[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.2.Lo[r0.2.Lo < twoLoc.Lo] ~ sm[r0.2.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi[twoLoc.Hi < 1] ~ sm[twoLoc.Hi < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'K', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.2.Hi)
        rm(r0.2.Lo)


    ##  Panel 12: C = 0.75
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.75)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.75)
        
        r0.2.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.2, C=0.75)
        r0.2.Hi[r0.2.Hi > 1]  <-  1.00000001
        r0.2.Hi[r0.2.Hi == 'NaN']  <-  1.00000001
        r0.2.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.2, C=0.75)
        
        pAgree  <-  rounded(sum(C0.75.h.25$agree[C0.75.h.25$r == 0.2])/length(C0.75.h.25$agree[C0.75.h.25$r == 0.2]), precision=3)
        pSim    <-  rounded(length(C0.75.h.25$sf[C0.75.h.25$Poly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.2]) / 
                            length(C0.75.h.25$sf[C0.75.h.25$r  ==  0.2]), precision=3)
        pEig    <-  rounded(length(C0.75.h.25$sf[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.2]) / 
                            length(C0.75.h.25$sf[C0.75.h.25$r  ==  0.2]), precision=3)

        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.75.h.25$sf[C0.75.h.25$Poly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.1] ~
               C0.75.h.25$sm[C0.75.h.25$Poly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.75.h.25$sf[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.1] ~
               C0.75.h.25$sm[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.75.h.25$sf[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 1 & C0.75.h.25$r  ==  0.1] ~
               C0.75.h.25$sm[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 1 & C0.75.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.2.Hi[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1] ~ sm[r0.2.Hi > twoLoc.Hi & r0.2.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.2.Lo[r0.2.Lo < twoLoc.Lo] ~ sm[r0.2.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi[twoLoc.Hi < 1] ~ sm[twoLoc.Hi < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'L', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.2.Hi)
        rm(r0.2.Lo)






##  Row 4: r = 0.2
    ##  Panel 13: C = 0
        
        # Calculate plotting lines for solutions not involving recombination
        twoLoc.Hi.obOut   <-  inv.lab1.domRev.obOut(hf=0.25, hm=0.25, sm=sm)
        twoLoc.Hi.obOut[twoLoc.Hi.obOut > 1]  <-  1.00000001
        twoLoc.Hi.obOut[10001]  <-  1.00000001
        twoLoc.Lo.obOut  <-  inv.lAB1.domRev.obOut(hf=0.25, hm=0.25, sm=sm)
        
        r0.5.Hi  <-  inv.lab2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0.5)
        r0.5.Hi[r0.5.Hi > 1]  <-  1.00000001
        r0.5.Hi[r0.5.Hi == 'NaN']  <-  1.00000001
        r0.5.Lo  <-  inv.lAB2.domRev.obOut(hf = 0.25, hm = 0.25, sm=sm, r=0.5)

        pAgree  <-  rounded(sum(C0.0.h.25$agree[C0.0.h.25$r == 0.2])/length(C0.0.h.25$agree[C0.0.h.25$r == 0.5]), precision=3)
        pSim    <-  rounded(length(C0.0.h.25$sf[C0.0.h.25$Poly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.5]) / 
                            length(C0.0.h.25$sf[C0.0.h.25$r  ==  0.5]), precision=3)
        pEig    <-  rounded(length(C0.0.h.25$sf[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.5]) / 
                            length(C0.0.h.25$sf[C0.0.h.25$r  ==  0.5]), precision=3)

        # Make the plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.0.h.25$sf[C0.0.h.25$Poly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.1] ~
               C0.0.h.25$sm[C0.0.h.25$Poly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.0.h.25$sf[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.1] ~
               C0.0.h.25$sm[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 0 & C0.0.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[5])
        points(C0.0.h.25$sf[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 1 & C0.0.h.25$r  ==  0.1] ~
               C0.0.h.25$sm[C0.0.h.25$eigPoly == 1 & C0.0.h.25$agree == 1 & C0.0.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.5.Hi[r0.5.Hi > twoLoc.Hi.obOut & r0.5.Hi <= 1] ~ sm[r0.5.Hi > twoLoc.Hi.obOut & r0.5.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.5.Lo[r0.5.Lo < twoLoc.Lo.obOut] ~ sm[r0.5.Lo < twoLoc.Lo.obOut], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi.obOut[twoLoc.Hi.obOut < 1] ~ sm[twoLoc.Hi.obOut < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo.obOut ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(0.05, 1.075, 'M', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.6, 0.5, expression(paste(italic(r), " = 0.5")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(-0.4, 0.5, expression(paste(italic(s[f]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.5, -0.4, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.68, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.68, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.68, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi.obOut)
        rm(twoLoc.Lo.obOut)
        rm(r0.5.Hi)
        rm(r0.5.Lo)


    ##  Panel 14: C = 0.25
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.25)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.25)
        
        r0.5.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.5, C=0.25)
        r0.5.Hi[r0.5.Hi > 1]  <-  1.00000001
        r0.5.Hi[r0.5.Hi == 'NaN']  <-  1.00000001
        r0.5.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.5, C=0.25)
        
        pAgree  <-  rounded(sum(C0.25.h.25$agree[C0.25.h.25$r == 0.2])/length(C0.25.h.25$agree[C0.25.h.25$r == 0.5]), precision=3)
        pSim    <-  rounded(length(C0.25.h.25$sf[C0.25.h.25$Poly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.5]) / 
                            length(C0.25.h.25$sf[C0.25.h.25$r  ==  0.5]), precision=3)
        pEig    <-  rounded(length(C0.25.h.25$sf[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.5]) / 
                            length(C0.25.h.25$sf[C0.25.h.25$r  ==  0.5]), precision=3)

        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.25.h.25$sf[C0.25.h.25$Poly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.1] ~
               C0.25.h.25$sm[C0.25.h.25$Poly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.25.h.25$sf[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.1] ~
               C0.25.h.25$sm[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 0 & C0.25.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.25.h.25$sf[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 1 & C0.25.h.25$r  ==  0.1] ~
               C0.25.h.25$sm[C0.25.h.25$eigPoly == 1 & C0.25.h.25$agree == 1 & C0.25.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.5.Hi[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1] ~ sm[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.5.Lo[r0.5.Lo < twoLoc.Lo] ~ sm[r0.5.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi[twoLoc.Hi < 1] ~ sm[twoLoc.Hi < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'N', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.4, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.5.Hi)
        rm(r0.5.Lo)



    ##  Panel 15: C = 0.5
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.5)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.5)
        
        r0.5.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.5, C=0.5)
        r0.5.Hi[r0.5.Hi > 1]  <-  1.00000001
        r0.5.Hi[r0.5.Hi == 'NaN']  <-  1.00000001
        r0.5.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.5, C=0.5)
        
        pAgree  <-  rounded(sum(C0.5.h.25$agree[C0.5.h.25$r == 0.2])/length(C0.5.h.25$agree[C0.5.h.25$r == 0.5]), precision=3)
        pSim    <-  rounded(length(C0.5.h.25$sf[C0.5.h.25$Poly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.5]) / 
                            length(C0.5.h.25$sf[C0.5.h.25$r  ==  0.5]), precision=3)
        pEig    <-  rounded(length(C0.5.h.25$sf[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.5]) / 
                            length(C0.5.h.25$sf[C0.5.h.25$r  ==  0.5]), precision=3)

        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.5.h.25$sf[C0.5.h.25$Poly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.1] ~
               C0.5.h.25$sm[C0.5.h.25$Poly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.5.h.25$sf[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.1] ~
               C0.5.h.25$sm[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 0 & C0.5.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.5.h.25$sf[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 1 & C0.5.h.25$r  ==  0.1] ~
               C0.5.h.25$sm[C0.5.h.25$eigPoly == 1 & C0.5.h.25$agree == 1 & C0.5.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.5.Hi[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1] ~ sm[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.5.Lo[r0.5.Lo < twoLoc.Lo] ~ sm[r0.5.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi[twoLoc.Hi < 1] ~ sm[twoLoc.Hi < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'O', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.4, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.5.Hi)
        rm(r0.5.Lo)


    ##  Panel 16: C = 0.75
        
        # Calculate plotting lines for solutions not involving recombination

        twoLoc.Hi   <-  inv.lab1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.75)
        twoLoc.Hi[twoLoc.Hi > 1]  <-  1.00000001
        twoLoc.Hi[10001]  <-  1.00000001
        twoLoc.Lo  <-  inv.lAB1.domRev(hf=0.25, hm=0.25, sm=sm, C=0.75)
        
        r0.5.Hi  <-  inv.lab2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.5, C=0.75)
        r0.5.Hi[r0.5.Hi > 1]  <-  1.00000001
        r0.5.Hi[r0.5.Hi == 'NaN']  <-  1.00000001
        r0.5.Lo  <-  inv.lAB2.domRev(hf = 0.25, hm = 0.25, sm=sm, r=0.5, C=0.75)
        
        pAgree  <-  rounded(sum(C0.75.h.25$agree[C0.75.h.25$r == 0.2])/length(C0.75.h.25$agree[C0.75.h.25$r == 0.5]), precision=3)
        pSim    <-  rounded(length(C0.75.h.25$sf[C0.75.h.25$Poly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.5]) / 
                            length(C0.75.h.25$sf[C0.75.h.25$r  ==  0.5]), precision=3)
        pEig    <-  rounded(length(C0.75.h.25$sf[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.5]) / 
                            length(C0.75.h.25$sf[C0.75.h.25$r  ==  0.5]), precision=3)

        # Make the plot
        #par(omi=rep(0.5, 4), mar = c(3,3,0.5,0.5), bty='o', xaxt='s', yaxt='s')
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Simulation points
        points(C0.75.h.25$sf[C0.75.h.25$Poly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.1] ~
               C0.75.h.25$sm[C0.75.h.25$Poly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[2])
        points(C0.75.h.25$sf[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.1] ~
               C0.75.h.25$sm[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 0 & C0.75.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[3])
        points(C0.75.h.25$sf[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 1 & C0.75.h.25$r  ==  0.1] ~
               C0.75.h.25$sm[C0.75.h.25$eigPoly == 1 & C0.75.h.25$agree == 1 & C0.75.h.25$r  ==  0.1], pch=21, col=NA, cex=0.75, bg=COLS[1])
        # Overlay analytic solutions for single locus & r=0 cases
            #  w/ recombination
            lines(r0.5.Hi[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1] ~ sm[r0.5.Hi > twoLoc.Hi & r0.5.Hi <= 1], lwd=1, col=COLS[4], lty=1)
            lines(r0.5.Lo[r0.5.Lo < twoLoc.Lo] ~ sm[r0.5.Lo < twoLoc.Lo], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Hi[twoLoc.Hi < 1] ~ sm[twoLoc.Hi < 1], lwd=1, col=COLS[4], lty=1)
            lines(twoLoc.Lo ~ sm, lwd=1, col=COLS[4], lty=1)
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        proportionalLabel(0.05, 1.075, 'P', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, -0.4, expression(paste(italic(s[m]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.95, substitute(p~" Agree", list(p = pAgree)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.88, substitute(p~" Sim.", list(p = pSim)), cex=0.55, adj=c(0, 0.5), xpd=NA)
        proportionalLabel(0.05, 0.80, substitute(p~" Eig.", list(p = pEig)), cex=0.55, adj=c(0, 0.5), xpd=NA)

        # Garbage collection
        rm(twoLoc.Hi)
        rm(twoLoc.Lo)
        rm(r0.5.Hi)
        rm(r0.5.Lo)


}
