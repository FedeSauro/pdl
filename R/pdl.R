### DA FARE:
# - scrivere autonomamente codice per la stima HAC
# - funzione per pseudoprevisioni con possibilita' di impostare l'orizzonte


## auxfit, auxfitL1, auxfitL2: stima modello con
##   parametri esterni noti (funzioni ausiliarie)
##
## L, L1, L2: funzioni kernel
##
## pdl, pdl_L2, pdl_L2: stima modello compresi parametri esterni


##############################################
###  FUNZIONI PER POLINOMIO NON VINCOLATO  ###
##############################################

## FUNZIONE PER OTTENERE LA MATRICE DEI RITARDI
#
# x: dati della x
# from, to: ritardi considerati (from=to=0 => no ritardi)
# ndrop: numero di osservazioni iniziali da porre a missing
#  (serve se si confrontano modelli con p diverso)
#
L <- function (x, from=0, to=0, ndrop=0) {
  # error
  if (to<from) stop("'from' must be lower than 'to'")
  
  # body
  if(to>0) {
    n = length(x)
    M = x
    for (i in 1:to) {
      M <- cbind(M, c(rep(NA,i),x[1:(n-i)]))
      }
    colnames(M) <- 0:to
    if(from>0 & from<=to) M <- M[,(from+1):(to+1),drop=F]
    if(ndrop>0) M[1:ndrop,] <- NA
    return(M)
    } else {
    return(x)  
    }
  }

## FUNZIONE PER COSTRUIRE LA MATRICE H
#
# from, to: ritardi considerati
# d: grado del polinomio
#
Hmat <- function(from, to, d) {
  M <- c()
  for(i in 0:d) {
    M <- cbind(M,(from:to)^i)
    }
  return(M)
  }

## FUNZIONE PER OTTENERE IL KERNEL
#
# x: dati della x
# from, to: ritardi considerati (from=to=0 => no ritardi)
# d: grado del polinomio (d=0 => ritardi non vincolati)
# ndrop: numero di osservazioni iniziali da porre a missing
#  (serve se si confrontano modelli con p diverso)
#
L <- function(x, from=0, to=0, d=0, ndrop=0) {
  # error 
  if (d > to) stop("'d' must be lower than 'to'")
  
  # body
  M <- L(x=x, from=from, to=to, ndrop=ndrop)
  if(d>0 & to>0) {
    return(M%*%Hmat(from=from, to=to, d=d))
    } else {
    return(M)
    }
  }

## FUNZIONE PER STIMARE UN MODELLO PDL
#
# y.name: nome della y
# x.names: nomi delle x con lags
# z.names: nomi delle x senza lags
# data: dataset da cui prendere le variabili
# p: vettore contenente il numero di ritardi per ciascuna variabile in x.names
# d: vettore contenente il grado del polinomio per ciascuna variabile in x.names
#
auxfit <- function(y.name, x.names, z.names=NULL, data, from=NULL, to=NULL, d=NULL, ndrop=0) {
  # se 'from' o 'to' hanno lunghezza minore di x.names, riempio con 0
  if(length(from)<length(x.names)) from[(length(from)+1):length(x.names)] <- 0
  if(length(to)<length(x.names)) to[(length(to)+1):length(x.names)] <- 0
  # se d ha lunghezza minore di x.names, riempio con 0
  if(length(d)<length(x.names)) d[(length(d)+1):length(x.names)] <- 0
  # crea formula
  if(ndrop==0) {
    form <- paste(y.name,"~",paste("L(",x.names,",",from,",",to,",d=",d,")",sep="",collapse="+"),sep="")    
    } else {
    form <- paste(y.name,"~",paste("L(",x.names,",",from,",",to,",d=",d,",ndrop=",ndrop,")",sep="",collapse="+"),sep="")
    }
  if(!is.null(z.names)) {
    form <- paste(form,"+",paste(z.names,collapse="+"))
    }
  mod <- lm(formula(form),data=data)
  # aggiungo info
  mod$call$formula <- formula(form)
  mod$from <- from
  mod$to <- to
  mod$d <- d
  names(mod$from) <- names(mod$to) <- names(mod$d) <- x.names
  mod$ndrop <- ndrop
  mod$variables <- list(y.name=y.name,x.names=x.names,z.names=z.names)
  # aggiungo classe 'pdl'
  class(mod) <- c("pdl","lm")
  mod
  }

## METODO SUMMARY PER CLASSE 'pdl'
summary.pdl <- function(object, cumulative=FALSE, alpha_heterosch = 0.05, ...) {
  #HAC = FALSE
  #if (object$ylag > 0) {
  #  HAC = TRUE
  #  lag_HAC = object$ylag
  #} else if (heterosch(object)<alpha_heterosch) {
  #  HAC = TRUE
  #  lag_HAC = 0
  #}
  summ <- summary.lm(object)
  #if (HAC) {
  #  HAC_COV <- NeweyWest(object, lag = lag_HAC)
  #  tab <- coeftest(object, vcov. = HAC_COV)
  #} else {
    tab <- summ$coefficients
  #}  
  summ$coefficients <- tab
  S <- vcov(object)
  #
  xnomi <- object$variables$x.names
  znomi <- object$variables$z.names
  nomiTab <- c("(Intercept)",znomi)
  tabOK <- tab[nomiTab,,drop=F]
  AUX <- 1
  for(i in 1:length(xnomi)) {
    from_i <- object$from[xnomi[i]]
    to_i <- object$to[xnomi[i]]
    if(to_i>0) {
      d_i <- object$d[xnomi[i]]
      } else {
      d_i <- 0
      }
    if(d_i>0) {
      ind_i <- (AUX+1):(AUX+1+d_i)
      theta_i <- tab[ind_i,1]
      S_i <- S[ind_i,ind_i]
      H_i <- Hmat(from=from_i, to=to_i, d=d_i)
      beta_i <- H_i%*%theta_i
      V_i <- H_i%*%S_i%*%t(H_i)
      se_i <- sqrt(diag(V_i))
      if(cumulative) {  ## <-----
        ###
        beta_i <- cumsum(beta_i)
        se_i <- c()
        for(j in 1:length(beta_i)) se_i[j] <- sqrt(sum(V_i[1:j,1:j]))
        ###
        }
      } else {
      ind_i <- (AUX+1):(AUX+1+to_i-from_i)
      beta_i <- tab[ind_i,1,drop=F]
      V_i <- S[ind_i,ind_i,drop=F]
      se_i <- sqrt(diag(V_i))
      if(cumulative) {  ## <-----
        ###
        beta_i <- cumsum(beta_i)
        se_i <- c()
        for(j in 1:length(beta_i)) se_i[j] <- sqrt(sum(V_i[1:j,1:j]))
        ###
        }
      }
    AUX <- max(ind_i)
    if(object$ndrop==0) {
      nomiTab <- c(nomiTab, paste("L(",xnomi[i],", ",from_i,", ",to_i,", ",d_i,")",from_i:to_i,sep=""))
      } else {
      nomiTab <- c(nomiTab, paste("L(",xnomi[i],", ",from_i,", ",to_i,", ",d_i,", ndrop = ",object$ndrop,")",from_i:to_i,sep=""))
      }
    zval_i <- beta_i/se_i
    pval_i <- 2*pnorm(-abs(zval_i))
    tabOK <- rbind(tabOK,cbind(beta_i,se_i,zval_i,pval_i))
    rownames(tabOK) <- nomiTab
    }
  summ$coefficients <- tabOK
  summ
  }


##  METODO CONFINT PER CLASSE 'pdl'
confint.pdl <- function(object, parm, level=0.95, cumulative=FALSE, ...) {
  ttab <- summary(object, cumulative=cumulative, ...)$coefficients
  tquan <- qt((1+level)/2, object$df.residual)
  ci <- data.frame(Estimate=ttab[,1], ttab[,1]-ttab[,2]*tquan, ttab[,1]+ttab[,2]*tquan)
  colnames(ci)[2:3] <- paste(100*c(1-level,1+level)/2,"%",sep="")
  ci
  }

##  METODO PLOT PER CLASSE 'pdl'
plot.pdl <- function(x, x.name=NULL, conf=0.95, interpolation=TRUE,
                     xlim=NULL, ylim=NULL, xlab="", ylab="", main = "", 
                     legend = "", add.grid=TRUE, ggp = FALSE, ...) {
  if(is.null(x.name)) x.name <- x$variables$x.names[1]  ## prendo la prima X
  lagsx <- x$from[x.name]
  lagdx <- x$to[x.name]
  ttab <- summary(x)$coefficients  ## coefficienti stimati
  #
  blab <- paste("L\\(",x.name,",",sep="")
  ttabOK <- ttab[grep(blab,rownames(ttab)),,drop=F]  ## matching parziale dei nomi di riga
  #
  bhat <- sehat <- rep(0,lagdx+1)
  bhat[(lagsx+1):(lagdx+1)] <- ttabOK[,1]
  sehat[(lagsx+1):(lagdx+1)] <- ttabOK[,2]
  tquan <- -qt((1-conf)/2,x$df.residual)
  
  if (main == "") main <- paste(x.name,"- degree ", x$d[x.name])
  if (xlab == "") xlab <- "Lags"
  if (ylab == "") ylab <- "Coefficients"
  if (legend == "") legend <- ""
  
  if (!ggp) {
    if(interpolation) {
      lagseq <- seq(0,lagdx,length=100)  ## griglia di ritardi
      bseq <- spline(0:(lagdx), bhat, xout=lagseq)$y  ## interpolazione
      bseq_sx <- spline(0:(lagdx), bhat-tquan*sehat, xout=lagseq)$y
      bseq_dx <- spline(0:(lagdx), bhat+tquan*sehat, xout=lagseq)$y
    } else {
      lagseq <- 0:(lagdx)
      bseq <- bhat
      bseq_sx <- bhat-tquan*sehat
      bseq_dx <- bhat+tquan*sehat
    }
    
    if(is.null(xlim)) xlim <- c(0,lagdx)
    if(is.null(ylim)) ylim <- range(c(bseq_sx,bseq_dx))
    plot(lagseq, bseq, type="n", ylim=ylim, xlim=xlim, main = main, xlab=xlab, ylab=ylab, ...)
    if(add.grid) grid()
    abline(h=0, lty=2)
    
    if(interpolation) {
      lines(lagseq, bseq, lwd=2)  ## curva interpolante
      #points(lagsx:lagdx, bhat[(lagsx+1):(lagdx+1)], cex=1.3, lwd=2)  ## punti sui ritardi interi
      polygon(c(lagseq,rev(lagseq)), c(bseq_sx,rev(bseq_dx)),
              border=NA, col=adjustcolor("grey75",alpha.f=0.5))  ## regione di confidenza
    } else {
      ind <- which(bseq!=0)
      arrows(lagseq[ind], bseq_sx[ind], lagseq[ind], bseq_dx[ind],
             code=3, angle=90, length=0.05)
      points(lagseq[ind], bseq[ind])
    }
    legend("topright", legend = legend, cex=0.9, bty="n")
    box()

  } else {
    if(is.null(xlim)) xlim <- c(0,lagdx)
    
    if(interpolation) {
      dataplot <- data.frame(
        lagseq <- seq(0,lagdx,length=100),  ## griglia di ritardi
        bseq <- spline(0:lagdx, bhat, xout=lagseq)$y,  ## interpolazione
        bseq_sx <- spline(0:lagdx, bhat-tquan*sehat, xout=lagseq)$y,
        bseq_dx <- spline(0:lagdx, bhat+tquan*sehat, xout=lagseq)$y
      )
      ggplot(dataplot, aes(x=lagseq,y=bseq, group = 1)) +
        geom_line(col = 'darkblue')+
        geom_line(y=0, col = 'black', lty = 2)+
        geom_ribbon(aes(ymin=bseq_sx, ymax=bseq_dx), alpha=0.5)+
        xlim(xlim)+
        ylim(range(c(bseq_sx,bseq_dx)))+
        labs(title = main, x = xlab, y = ylab)+
        theme(plot.title = element_text(hjust = 0.5))+
        annotate(geom = "text", x = (xlim[2]+xlim[1])*0.7, y = range(c(bseq_sx,bseq_dx))[2], 
                 label = legend, hjust = 0, vjust = 1, size = 4)
    } else {
      
      dataplot <- data.frame(
        lagseq <- 0:lagdx,
        bseq <- bhat,
        bseq_sx <- bhat-tquan*sehat,
        bseq_dx <- bhat+tquan*sehat
      )  
      ggplot(dataplot, aes(x=factor(lagseq), y=bseq, group = 1)) +
        geom_point(size = 3, color = "darkblue")+
        geom_line(y=0, col = 'black', lty = 2)+
        geom_errorbar( aes(ymin= bseq_sx, ymax=bseq_dx), width=0.4, colour="darkblue", size=1)+
        labs(title = main, x = xlab, y = ylab)+
        theme(plot.title = element_text(hjust = 0.5))+
        annotate(geom = "text", x = (xlim[2]+xlim[1])*0.7, y = range(c(bseq_sx,bseq_dx))[2], 
                 label = legend, hjust = 0, vjust = 1, size = 4)
      }
    }
  }

## funzione per calcolo BIC esatto
# N.B. la funzione BIC() non considera la costante di normalizzazione
bicCalc <- function(x) {
  extractAIC(x, k=log(nobs(x)))[2]
  }

### SELEZIONE LAG DI X E DI Y  (VISIBILE) <---------------------
# in questa funzione inserire x.name solo 1 variabile alla volta
#
pdl <- function(data, y.name, x.name, z.names=NULL, from.max=0, to.max, d=0, ylag.max = 0, quiet=FALSE) {
  ### modello nullo (senza X)
  #
  #if(length(z.names)>0) {
  #  formNull <- paste(y.name,"~",paste(z.names,collapse="+"),sep="")
  #  } else {
  #  formNull <- paste(y.name,"~1",sep="")
  #  }
  #actualMod <- lm(formula(formNull), data=data)
  #
  actualBic <- Inf
  for(i in 0:from.max) {
    for(j in i:to.max) {
      dnew <- ifelse(d > j, j, d)
      ijmod <- auxfit(y.name=y.name, x.names=x.name, z.names=z.names,
                      from=i, to=j, d=dnew, ndrop=to.max, data=data)
      ijbic <- bicCalc(ijmod)
      if(ijbic<actualBic) {
        xlag <- c(i,j, dnew)
        actualBic <- ijbic
        actualMod <- ijmod
        }
      }
    }
  #
  ### selezione ordine di autocorrelazione residua
  if (ylag.max > 0) {
    ybic <- sapply(1:to.max, function(x){
      bicCalc(lm(actualMod$residual~L(actualMod$residual, from=1, to=x)))
      })
    ybic <- c(bicCalc(lm(actualMod$residual~1)),ybic)
    ylag.max <- which.min(ybic)-1
  }
    #
  ### modello finale con ndrop=0
  if(ylag.max>0) {
    mod <- auxfit(y.name=y.name, x.names=c(y.name,x.name), z.names=z.names,
                  from=c(0,xlag[1]), to=c(ylag.max,xlag[2]), d=c(0,xlag[3]), ndrop=0, data=data) 
    } else {
    mod <- auxfit(y.name=y.name, x.names=x.name, z.names=z.names,
                  from=xlag[1], to=xlag[2], d=xlag[3], ndrop=0, data=data) 
    }
  mod$ylag <- ylag.max
  return(mod)
  }


##############################################
###  FUNZIONI PER LAG LINEARE DECRESCENTE  ###
##############################################

## FUNZIONE per costruire i pesi
# c e' l'ultimo ritardo che vogliamo considerare
L1_weights <- function(lag, c) {
  kost <- 1+c/2
  wei <- sapply(lag, function(k){
    (1-k/(c+1))/kost*(k>=0&k<(c+1))
  })
  names(wei) <- lag
  wei
}

# L1 kernel projection
L1 <- function(x, c, ndrop=0) {
  n <- length(x)
  wei <- L1_weights(0:n, c=c)
  res <- c()
  for(i in (c+1):n) {
    res[i] <- sum(wei[1:i]*x[i-(0:(i-1))])
  }
  if(ndrop>0) res[1:ndrop] <- NA
  res
}

##  FUNZIONE PER STIMARE UN MODELLO LAG LINEARE DECRESCENTE
# y.name: nome della y
# x.names: nomi delle x con lags
# z.names: nomi delle x senza lags
# data: dataset da cui prendere le variabili
# c: estremo superiore vincolato
#
auxfitL1 <- function(y.name, x.names, z.names=NULL, c, data, ndrop=0) {
  # crea formula
  if(ndrop==0) {
    form <- paste(y.name,"~",paste("L1(",x.names,",",c,")",sep="",collapse="+"),sep="")    
  } else {
    form <- paste(y.name,"~",paste("L1(",x.names,",",c,",ndrop=",ndrop,")",sep="",collapse="+"),sep="")
  }
  if(!is.null(z.names)) {
    form <- paste(form,"+",paste(z.names,collapse="+"))
  }
  mod <- lm(formula(form),data=data)
  # aggiungo info
  mod$call$formula <- formula(form)
  mod$c <- c
  mod$ndrop <- ndrop
  mod$variables <- list(y.name=y.name,x.names=x.names,z.names=z.names)
  # aggiungo classe 'L1'
  class(mod) <- c("pdl_L1","lm")
  mod
}

##  METODO SUMMARY PER CLASSE 'L1'
summary.pdl_L1 <- function(object, cumulative=FALSE, alpha_heterosch = 0.05, ...) {
  #HAC = FALSE
  #if (object$ylag > 0) {
  #  HAC = TRUE
  #  lag_HAC = object$ylag
  #} else if (heterosch(object)<alpha_heterosch) {
  #  HAC = TRUE
  #  lag_HAC = 0
  #}
  summ <- summary.lm(object)
  #if (HAC) {
  #  HAC_COV <- NeweyWest(object, lag = lag_HAC)
  #  tab <- coeftest(object, vcov. = HAC_COV)
  #} else {
  tab <- summ$coefficients
  #}
  summ$coefficients <- tab
  #
  xnomi <- object$variables$x.names
  znomi <- object$variables$z.names
  nomiTab <- c("(Intercept)",znomi)
  tabOK <- tab[nomiTab,,drop=F]
  AUX <- 1
  for(i in 1:length(xnomi)) {
    c <- object$c
    ind_i <- (AUX+1):(AUX+1)
    theta_i <- tab[ind_i,1]
    h_i <- L1_weights(0:object$c, object$c)
    beta_i <- h_i*theta_i
    se_i <- h_i*tab[ind_i, 2]
    
    if(cumulative) {  ## <-----
      ###
      beta_i <- cumsum(beta_i)
      se_i <- cumsum(se_i)
    }
    AUX <- max(ind_i)
    if(object$ndrop==0) {
      nomiTab <- c(nomiTab, paste("L1(",xnomi[i],", ",object$c,")",0:object$c,sep=""))
    } else {
      nomiTab <- c(nomiTab, paste("L1(",xnomi[i],", ",object$c,", ndrop = ",object$ndrop,")",0:object$c,sep=""))
    }
    zval_i <- beta_i/se_i
    pval_i <- 2*pnorm(-abs(zval_i))
    tabOK <- rbind(tabOK,cbind(beta_i,se_i,zval_i,pval_i))
    rownames(tabOK) <- nomiTab
  }
  summ$coefficients <- tabOK
  summ
}

## METODO CONFINT PER CLASSE 'L1'
confint.pdl_L1 <- function(object, parm, level=0.95, cumulative=FALSE, ...) {
  ttab <- summary(object, cumulative=cumulative, ...)$coefficients
  tquan <- qt((1+level)/2, object$df.residual)
  ci <- data.frame(Estimate=ttab[,1], ttab[,1]-ttab[,2]*tquan, ttab[,1]+ttab[,2]*tquan)
  colnames(ci)[2:3] <- paste(100*c(1-level,1+level)/2,"%",sep="")
  ci
}

## METODO PLOT PER CLASSE 'L1'
plot.pdl_L1 <- function(x, x.name=NULL, conf=0.95, interpolation=TRUE,
                    xlim=NULL, ylim=NULL, xlab="", ylab="", main = "", 
                    legend = "", add.grid=TRUE, ggp = FALSE, ...) {
  if(is.null(x.name)) x.name <- x$variables$x.names[1]  ## prendo la prima X
  lagsx <- 0
  lagdx <- x$c
  ttab <- summary(x)$coefficients  ## coefficienti stimati
  #
  blab <- paste("L1\\(",x.name,",",sep="")
  ttabOK <- ttab[grep(blab,rownames(ttab)),,drop=F]  ## matching parziale dei nomi di riga
  #
  bhat <- sehat <- rep(0,(lagdx+1))
  bhat[0:(lagdx+1)] <- ttabOK[,1]
  sehat[0:(lagdx+1)] <- ttabOK[,2]
  tquan <- -qt((1-conf)/2,x$df.residual)
  
  #if (main == "") main <- paste(x.name,"- subjected to linear type")
  if (xlab == "") xlab <- "Lags"
  if (ylab == "") ylab <- "Coefficients"
  if (legend == "") legend <- ""
  
  if (!ggp) {
    if(interpolation) {
      lagseq <- seq(0,lagdx+1,length=100)  ## griglia di ritardi
      bseq <- spline(0:lagdx+1, bhat, xout=lagseq)$y  ## interpolazione
      bseq_sx <- spline(0:lagdx+1, bhat-tquan*sehat, xout=lagseq)$y
      bseq_dx <- spline(0:lagdx+1, bhat+tquan*sehat, xout=lagseq)$y
    } else {
      lagseq <- 0:(lagdx+1)
      bseq <- c(bhat, 0)
      bseq_sx <- c(bhat-tquan*sehat, 0)
      bseq_dx <- c(bhat+tquan*sehat, 0)
    }
    
    if(is.null(xlim)) xlim <- c(0,lagdx+1)
    if(is.null(ylim)) ylim <- range(c(bseq_sx,bseq_dx))
    plot(lagseq, bseq, type="n", ylim=ylim, xlim=xlim, main = main, xlab=xlab, ylab=ylab)
    if(add.grid) grid()
    abline(h=0, lty=2)
    
    if(interpolation) {
      lines(lagseq, bseq, lwd=2)  ## curva interpolante
      #points(lagsx:lagdx, bhat[(lagsx+1):(lagdx+1)], cex=1.3, lwd=2)  ## punti sui ritardi interi
      polygon(c(lagseq,rev(lagseq)), c(bseq_sx,rev(bseq_dx)),
              border=NA, col=adjustcolor("grey75",alpha.f=0.5))  ## regione di confidenza
    } else {
      ind <- which(bseq!=0)
      arrows(lagseq[ind], bseq_sx[ind], lagseq[ind], bseq_dx[ind],
             code=3, angle=90, length=0.05)
      points(lagseq[ind], bseq[ind])
    }
    legend("topright", legend = legend, cex=0.9, bty="n")
    box()
    
  } else {
    if(is.null(xlim)) xlim <- c((lagsx-1),(lagdx+1))
    
    if(interpolation) {
      dataplot <- data.frame(
        lagseq <- seq(0,lagdx+1,length=100),  ## griglia di ritardi
        bseq <- spline(0:lagdx+1, bhat, xout=lagseq)$y,  ## interpolazione
        bseq_sx <- spline(0:lagdx+1, bhat-tquan*sehat, xout=lagseq)$y,
        bseq_dx <- spline(0:lagdx+1, bhat+tquan*sehat, xout=lagseq)$y
      )
      ggplot(dataplot, aes(x=lagseq,y=bseq, group = 1)) +
        geom_line(col = 'darkblue')+
        geom_line(y=0, col = 'black', lty = 2)+
        geom_ribbon(aes(ymin=bseq_sx, ymax=bseq_dx), alpha=0.5)+
        xlim(xlim)+
        ylim(range(c(bseq_sx,bseq_dx)))+
        labs(title = main, x = xlab, y = ylab)+
        theme(plot.title = element_text(hjust = 0.5))+
        annotate(geom = "text", x = (xlim[2]+xlim[1])*0.7, y = range(c(bseq_sx,bseq_dx))[2], 
                 label = legend, hjust = 0, vjust = 1, size = 4)
    } else {
      
      dataplot <- data.frame(
        lagseq <- 0:(lagdx+1),
        bseq <- c(bhat, 0),
        bseq_sx <- c(bhat-tquan*sehat, 0),
        bseq_dx <- c(bhat+tquan*sehat, 0)
      )  
      ggplot(dataplot, aes(x=factor(lagseq), y=bseq, group = 1)) +
        geom_point(size = 3, color = "darkblue")+
        geom_line(y=0, col = 'black', lty = 2)+
        geom_errorbar( aes(ymin= bseq_sx, ymax=bseq_dx), width=0.4, colour="darkblue", size=1)+
        labs(title = main, x = xlab, y = ylab)+
        theme(plot.title = element_text(hjust = 0.5))+
        annotate(geom = "text", x = (xlim[2]+xlim[1])*0.7, y = range(c(bseq_sx,bseq_dx))[2], 
                 label = legend, hjust = 0, vjust = 1, size = 4)
    }
  }
}


### SELEZIONE LAG DI X E DI Y PER LAG LINEARE DECRESCENTE (VISIBILE) <--------
# in questa funzione inserire x.name solo 1 variabile alla volta
#
pdl_L1 <- function(data, y.name, x.name, z.names=NULL, c, ylag.max = 0, quiet=FALSE) {
  ### modello nullo (senza X)
  #
  #if(length(z.names)>0) {
  #  formNull <- paste(y.name,"~",paste(z.names,collapse="+"),sep="")
  #  } else {
  #  formNull <- paste(y.name,"~1",sep="")
  #  }
  #actualMod <- lm(formula(formNull), data=data)
  #
  actualBic <- Inf
  for(i in 0:c) {
      
      imod <- auxfitL1(y.name=y.name, x.names=x.name, z.names=z.names, 
                         c=i,ndrop=i, data=data)
      ibic <- bicCalc(imod)
      
      if(ibic<actualBic) {
        
        xlag <- i
        actualBic <- ibic
        actualMod <- imod
      }
  }
  #
  ### selezione ordine di autocorrelazione residua
  if (ylag.max>0) {
    ybic <- sapply(1:ylag.max, function(x){
      bicCalc(lm(actualMod$residual~L(actualMod$residual, from=1, to=x)))
    })
    ybic <- c(bicCalc(lm(actualMod$residual~1)),ybic)
    ylag.max <- which.min(ybic)-1
    
  }
  #
  ### modello finale con ndrop=0
  if(ylag.max>0) {
    mod <- auxfitL1(y.name=y.name, x.names=c(y.name,x.name), z.names=z.names,
                     c=c(ylag.max,xlag), ndrop=0, data=data) 
  } else {
    mod <- auxfitL1(y.name=y.name, x.names=x.name, z.names=z.names,
                     c=xlag, ndrop=0, data=data) 
  }
  mod$ylag <- ylag.max
  return(mod)
}


########################################################
###  FUNZIONE PER LAG PARABOLICO CON ESTREMI FISSATI ###
########################################################

# funzione dei pesi
L2_weights <- function(lag, a, b) {
  kost <- (-1/6)*(b*(b+1)*(2*b+1)-(a-1)*a*(2*a-1))+(a+b)/2*(b*(b+1)-(a-1)*a)-(b-a+1)*(a-1)*(b+1)
  ## kost e' una costante di normalizzazione dei pesi h_k: in questo modo
  ## essi sommeranno a 1, quindi theta si interpreta subito come coeff cumulato
  ## l'ho trovata con la formula di gauss e di sicuro e' semplificabile, ma non ci ho mai provato
  wei <- sapply(lag, function(k){
    (-k^2+(a+b)*k-(a-1)*(b+1))/kost*(k>a-1&k<b+1)
  })
  names(wei) <- lag
  wei
}

# L2 kernel projection
L2 <- function(x, a, b, ndrop=0) {
  n <- length(x)
  wei <- L2_weights(0:n, a=a, b=b)
  res <- c()
  for(i in (b+1):n) {
    res[i] <- sum(wei[1:i]*x[i-(0:(i-1))])
  }
  if(ndrop>0) res[1:ndrop] <- NA
  res
}


###  FUNZIONE PER STIMARE UN MODELLO CON LAG PARABOLICO
# y.name: nome della y
# x.names: nomi delle x con lags
# z.names: nomi delle x senza lags
# data: dataset da cui prendere le variabili
# a: estremo inferiore vincolato
# b: estremo superiore vincolato
#
auxfitL2 <- function(y.name, x.names, z.names=NULL, a, b, data, ndrop=0) {
  # crea formula
  if(ndrop==0) {
    form <- paste(y.name,"~",paste("L2(",x.names,",",a,",",b,")",sep="",collapse="+"),sep="")    
    } else {
    form <- paste(y.name,"~",paste("L2(",x.names,",",a,",",b,",ndrop=",ndrop,")",sep="",collapse="+"),sep="")
    }
  if(!is.null(z.names)) {
    form <- paste(form,"+",paste(z.names,collapse="+"))
    }
  mod <- lm(formula(form),data=data)
  # aggiungo info
  mod$call$formula <- formula(form)
  mod$a <- a
  mod$b <- b
  mod$ndrop <- ndrop
  mod$variables <- list(y.name=y.name,x.names=x.names,z.names=z.names)
  # aggiungo classe 'L2'
  class(mod) <- c("pdl_L2","lm")
  mod
  }

###  METODO SUMMARY PER CLASSE 'L2'
summary.pdl_L2 <- function(object, cumulative=FALSE, alpha_heterosch = 0.05, ...) {
  #HAC = FALSE
  #if (object$ylag > 0) {
  #  HAC = TRUE
  #  lag_HAC = object$ylag
  #} else if (heterosch(object)<alpha_heterosch) {
  #  HAC = TRUE
  #  lag_HAC = 0
  #}
  summ <- summary.lm(object)
  #if (HAC) {
  #  HAC_COV <- NeweyWest(object, lag = lag_HAC)
  #  tab <- coeftest(object, vcov. = HAC_COV)
  #} else {
  tab <- summ$coefficients
  #}
  summ$coefficients <- tab
  #
  xnomi <- object$variables$x.names
  znomi <- object$variables$z.names
  nomiTab <- c("(Intercept)",znomi)
  tabOK <- tab[nomiTab,,drop=F]
  AUX <- 1
  for(i in 1:length(xnomi)) {
    a <- object$a
    b <- object$b
    ind_i <- (AUX+1):(AUX+1)
    theta_i <- tab[ind_i, 1]
    h_i <- L2_weights(object$a:object$b, object$a, object$b)
    beta_i <- h_i*theta_i
    se_i <- h_i*tab[ind_i, 2]
    
    if(cumulative) {  ## <-----
        ###
        beta_i <- cumsum(beta_i)
        se_i <- cumsum(se_i)
    }
  AUX <- max(ind_i)
  if(object$ndrop==0) {
    nomiTab <- c(nomiTab, paste("L2(",xnomi[i],", ",object$a,", ",object$b,")",object$a:object$b,sep=""))
    } else {
    nomiTab <- c(nomiTab, paste("L2(",xnomi[i],", ",object$a,", ",object$b,", ndrop = ",object$ndrop,")",object$a:object$b,sep=""))
  }
  zval_i <- beta_i/se_i
  pval_i <- 2*pnorm(-abs(zval_i))
  tabOK <- rbind(tabOK,cbind(beta_i,se_i,zval_i,pval_i))
  rownames(tabOK) <- nomiTab
  }
  summ$coefficients <- tabOK
  summ
}


###  METODO CONFINT PER CLASSE 'L2'
confint.pdl_L2 <- function(object, parm, level=0.95, cumulative=FALSE, ...) {
  ttab <- summary(object, cumulative=cumulative, ...)$coefficients
  tquan <- qt((1+level)/2, object$df.residual)
  ci <- data.frame(Estimate=ttab[,1], ttab[,1]-ttab[,2]*tquan, ttab[,1]+ttab[,2]*tquan)
  colnames(ci)[2:3] <- paste(100*c(1-level,1+level)/2,"%",sep="")
  ci
}


##  METODO PLOT PER CLASSE 'L2'
plot.pdl_L2 <- function(x, x.name=NULL, conf=0.95, interpolation=TRUE,
                       xlim=NULL, ylim=NULL, xlab="", ylab="", main = "", 
                       legend = "", add.grid=TRUE, ggp = FALSE, ...) {
  if(is.null(x.name)) x.name <- x$variables$x.names[1]  ## prendo la prima X
  lagsx <- x$a
  lagdx <- x$b
  ttab <- summary(x)$coefficients  ## coefficienti stimati
  #
  blab <- paste("L2\\(",x.name,",",sep="")
  ttabOK <- ttab[grep(blab,rownames(ttab)),,drop=F]  ## matching parziale dei nomi di riga
  #
  bhat <- sehat <- rep(0,(lagdx-lagsx+3))
  bhat[2:(lagdx-lagsx+2)] <- ttabOK[,1]
  sehat[2:(lagdx-lagsx+2)] <- ttabOK[,2]
  tquan <- -qt((1-conf)/2,x$df.residual)
  
  #if (main == "") main <- paste(x.name,"- subjected to parabolic type")
  if (xlab == "") xlab <- "Lags"
  if (ylab == "") ylab <- "Coefficients"
  if (legend == "") legend <- ""
  
  if (!ggp) {
    if(interpolation) {
      lagseq <- seq(lagsx-1,lagdx+1,length=100)  ## griglia di ritardi
      bseq <- spline((lagsx-1):(lagdx+1), bhat, xout=lagseq)$y  ## interpolazione
      bseq_sx <- spline((lagsx-1):(lagdx+1), bhat-tquan*sehat, xout=lagseq)$y
      bseq_dx <- spline((lagsx-1):(lagdx+1), bhat+tquan*sehat, xout=lagseq)$y
    } else {
      lagseq <- (lagsx-1):(lagdx+1)
      bseq <- bhat
      bseq_sx <- bhat-tquan*sehat
      bseq_dx <- bhat+tquan*sehat
    }
    
    if(is.null(xlim)) xlim <- c(lagsx-1,lagdx+1)
    if(is.null(ylim)) ylim <- range(c(bseq_sx,bseq_dx))
    plot(lagseq, bseq, type="n", ylim=ylim, xlim=xlim, main = main, xlab=xlab, ylab=ylab)
    if(add.grid) grid()
    abline(h=0, lty=2)
    
    if(interpolation) {
      lines(lagseq, bseq, lwd=2)  ## curva interpolante
      #points(lagsx:lagdx, bhat[(lagsx+1):(lagdx+1)], cex=1.3, lwd=2)  ## punti sui ritardi interi
      polygon(c(lagseq,rev(lagseq)), c(bseq_sx,rev(bseq_dx)),
              border=NA, col=adjustcolor("grey75",alpha.f=0.5))  ## regione di confidenza
    } else {
      ind <- which(bseq!=0)
      arrows(lagseq[ind], bseq_sx[ind], lagseq[ind], bseq_dx[ind],
             code=3, angle=90, length=0.05)
      points(lagseq[ind], bseq[ind])
    }
    legend("topright", legend = legend, cex=0.9, bty="n")
    box()
    
  } else {
    if(is.null(xlim)) xlim <- c((lagsx-1),(lagdx+1))
    
    if(interpolation) {
      dataplot <- data.frame(
        lagseq <- seq((lagsx-1),(lagdx+1),length=100),  ## griglia di ritardi
        bseq <- spline((lagsx-1):(lagdx+1), bhat, xout=lagseq)$y,  ## interpolazione
        bseq_sx <- spline((lagsx-1):(lagdx+1), bhat-tquan*sehat, xout=lagseq)$y,
        bseq_dx <- spline((lagsx-1):(lagdx+1), bhat+tquan*sehat, xout=lagseq)$y
      )
      ggplot(dataplot, aes(x=lagseq,y=bseq, group = 1)) +
        geom_line(col = 'darkblue')+
        geom_line(y=0, col = 'black', lty = 2)+
        geom_ribbon(aes(ymin=bseq_sx, ymax=bseq_dx), alpha=0.5)+
        xlim(xlim)+
        ylim(range(c(bseq_sx,bseq_dx)))+
        labs(title = main, x = xlab, y = ylab)+
        theme(plot.title = element_text(hjust = 0.5))+
        annotate(geom = "text", x = (xlim[2]+xlim[1])*0.7, y = range(c(bseq_sx,bseq_dx))[2], 
                 label = legend, hjust = 0, vjust = 1, size = 4)
    } else {
      
      dataplot <- data.frame(
        lagseq <- (lagsx-1):(lagdx+1),
        bseq <- bhat,
        bseq_sx <- bhat-tquan*sehat,
        bseq_dx <- bhat+tquan*sehat
      )  
      ggplot(dataplot, aes(x=factor(lagseq), y=bseq, group = 1)) +
        geom_point(size = 3, color = "darkblue")+
        geom_line(y=0, col = 'black', lty = 2)+
        geom_errorbar( aes(ymin= bseq_sx, ymax=bseq_dx), width=0.4, colour="darkblue", size=1)+
        labs(title = main, x = xlab, y = ylab)+
        theme(plot.title = element_text(hjust = 0.5))+
        annotate(geom = "text", x = (xlim[2]+xlim[1])*0.7, y = range(c(bseq_sx,bseq_dx))[2], 
                 label = legend, hjust = 0, vjust = 1, size = 4)
    }
  }
}


### SELEZIONE LAG DI X E DI Y PER IL LAG PARABOLICO (VISIBILE) <-------------
# in questa funzione inserire x.name solo 1 variabile alla volta
#
pdl_L2 <- function(data, y.name, x.name, z.names=NULL, a = 0, b, min.width=0, ylag.max = 0, quiet=FALSE) {
  ### modello nullo (senza X)
  #
  #if(length(z.names)>0) {
  #  formNull <- paste(y.name,"~",paste(z.names,collapse="+"),sep="")
  #  } else {
  #  formNull <- paste(y.name,"~1",sep="")
  #  }
  #actualMod <- lm(formula(formNull), data=data)
  #
  if (min.width>b) min.width <- b
  actualBic <- Inf
  for(i in 0:min(a,b-min.width)) {
    for(j in (i+min.width):b) {
      
      ijmod <- auxfitL2(y.name=y.name, x.names=x.name, z.names=z.names, 
                        a=i, b=j,ndrop=b, data=data)
      ijbic <- bicCalc(ijmod)
      
      if(ijbic<actualBic) {
      
        xlag <- c(i,j)
        actualBic <- ijbic
        actualMod <- ijmod
        }
      }
    }
  #
  ### selezione ordine di autocorrelazione residua
  if (ylag.max>0) {
    ybic <- sapply(1:ylag.max, function(x){
      bicCalc(lm(actualMod$residual~L(actualMod$residual, from=1, to=x)))
    })
    ybic <- c(bicCalc(lm(actualMod$residual~1)),ybic)
    ylag.max <- which.min(ybic)-1
    
  }
  #
  ### modello finale con ndrop=0
  if(ylag.max>0) {
    mod <- auxfitL2(y.name=y.name, x.names=c(y.name,x.name), z.names=z.names,
                  a=c(0,xlag[1]), b=c(ylag.max,xlag[2]), ndrop=0, data=data) 
    } else {
    mod <- auxfitL2(y.name=y.name, x.names=x.name, z.names=z.names,
                  a=xlag[1], b=xlag[2], ndrop=0, data=data) 
    }
  mod$ylag <- ylag.max
  return(mod)
}


########################################################
# TEST DI WHITE
# verifica eteroschedasticita' dei residui
# restituisce il p-value della statistica chiquadro
# rifiutiamo H0 che ci sia omoschedasticita'
heterosch <- function(object, ...) {
  e2 <- object$residuals^2
  Xstar <- model.matrix(object)[,2]
  reg <- lm(e2 ~ Xstar + I(Xstar)^2)
  chisq <- length(Xstar)*summary(reg)$r.squared
  1-pchisq(chisq, 2)  ## pvalue
}
