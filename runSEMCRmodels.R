# Blackbird example
# Cubaynes, S, C. Doutrelant, A. Gregoire, P. Perret, B. Faivre, and O. Gimenez (2012). Testing hypotheses in evolutionary ecology with imperfect detection: Structural equation modeling of mark-recapture data. Ecology 93: 248-255.

#------ fit sem cr model

# covariates
# 'Tarsus','Phalanx','Rectrice','Wing'
covF=matrix(c(32.96,9.44,99.38,123.13,
31.51,8.72,100.75,123.25,
33.22,9.54,96.63,121.38,
33.69,9.83,98.63,121.63,
32.20,8.78,93.38,120.25,
33.27,9.20,95.50,117.88,
32.11,9.32,97.00,120.00,
33.96,9.91,100.00,123.50,
33.20,9.55,99.00,122.00,
35.20,9.46,92.63,118.50,
32.49,8.89,96.00,119.25,
33.74,9.86,99.50,122.13,
33.19,9.10,99.50,122.38,
33.04,9.41,100.13,123.50,
34.42,9.58,94.75,120.50,
33.20,9.35,98.00,121.88,
33.50,9.51,100.00,120.50,
33.13,9.21,99.50,119.75,
32.36,9.43,95.88,118.50,
35.14,9.90,92.75,122.00,
34.00,9.70,103.75,124.00,
33.71,9.75,99.50,122.63,
33.39,9.58,97.88,120.75,
32.82,9.18,97.75,124.13,
33.05,9.21,93.88,119.88,
30.97,8.81,87.75,114.13,
33.92,9.63,97.25,122.75,
35.74,10.17,95.50,122.25,
33.11,9.15,99.63,118.75,
33.79,9.58,95.75,122.00,
33.78,9.31,98.88,122.50,
31.86,9.31,94.50,119.75,
33.73,10.01,100.25,124.00,
32.32,9.52,98.13,121.50,
33.80,9.85,97.00,124.25,
33.89,9.22,91.75,116.13,
33.74,9.39,100.13,122.13,
34.05,9.66,97.38,121.75,
33.97,9.76,100.13,121.75,
32.68,9.88,97.00,122.00,
33.76,9.53,95.50,120.38,
32.42,9.18,101.25,120.50,
32.78,9.77,99.25,127.38,
30.22,8.56,97.13,119.75,
31.30,8.78,96.25,116.00,
31.93,8.95,96.00,120.13,
33.11,9.04,102.30,124.00,
31.36,9.30,104.50,122.50,
32.05,9.78,100.50,124.10,
34.22,9.61,99.30,124.50,
33.80,9.57,88.30,122.30,
34.44,9.98,95.25,121.50,
33.76,9.86,102.00,123.30,
33.96,9.82,105.40,126.60,
32.68,9.33,98.50,121.50,
29.39,8.20,88.00,117.50,
32.48,9.67,99.10,120.40,
32.18,9.30,99.25,124.00,
33.26,9.52,97.00,118.75,
33.54,9.64,99.00,122.50,
32.26,9.56,100.13,126.00,
32.38,9.20,97.00,120.10,
32.69,9.24,94.00,116.63,
31.65,9.31,94.75,118.50,
34.39,9.68,102.88,126.38,
33.73,9.22,98.00,119.38,
33.13,9.49,100.38,121.00,
34.04,9.56,96.63,122.25,
32.68,9.07,108.00,125.00,
31.65,9.47,84.63,115.75,
32.31,9.15,90.13,116.00,
33.32,9.27,103.75,122.75,
32.27,9.35,98.25,121.88,
33.87,9.26,97.38,119.50,
33.25,9.53,101.25,123.50,
32.40,9.20,98.88,121.50,
33.38,9.48,93.00,122.13,
32.44,9.86,99.00,122.25,
33.58,9.90,96.88,125.00,
33.72,9.25,103.75,127.00,
32.31,9.05,101.13,123.13,
32.85,9.45,100.25,123.88,
31.48,9.02,101.25,121.50,
33.40,9.48,95.88,120.25),byrow=T,nrow=84)

cov = covF
mean.cov = apply(cov,2,mean)
mean.mat = matrix(rep(mean.cov,84),byrow=T,ncol=ncol(cov))
sd.cov = apply(cov,2,sd)
sd.mat = matrix(rep(sd.cov,84),byrow=T,ncol=ncol(cov))
cov = (cov-mean.mat)/sd.mat

## encounter histories
h <- read.table("SEMfemales.dat")

# number of individuals 
N <- dim(h)[[1]] 

# number of years
Years <- dim(h)[[2]]

# compute the date of first capture
First <- NULL
for (i in 1:N)
{
	temp <- 1:Years
	First <- c(First,min(temp[h[i,]==1]))
}

# init for the states
Xinit <- matrix(NA,nrow=N,ncol=Years)
for (i in 1:N)
{
	for (j in 1:(Years))
	{
		if (j > First[i]) Xinit[i,j] <- 1
	}
}

# first list of inits
init1 <- list(p = 0.2,alive=as.matrix(Xinit),sig12x = .1,sig34x = .1,sigksi1 = .1,gamma12 = 1,gamma13 = .5,gamma14 = .5,beta21 = -1,beta212 = -1)
# second list of inits
init2 <- list(p = 0.8,alive=as.matrix(Xinit),sig12x = 1,sig34x = 1,sigksi1 = 1,gamma12 = 1,gamma13 = .5,gamma14 = .5,beta21 = 1,beta212 = 1)
# concatenate list of initial values
inits <- list(init1,init2)

#- tarsus length (x1)
#- phalanx length (x2)
#- wing length (x3)
#- rectrice length (x4)

# Load rjags package
library(rjags)

# data
datax <- list(N=N,Years=Years,dat=as.matrix(h),First=First,cov=cov)

# store the starting point
deb = Sys.time()

# MCMC simulations 
model <- jags.model(file = "SEMMR.bug", data = datax, inits = inits, n.chains = 2) 
mcmc <- coda.samples(model, c("p","gamma12","gamma13","gamma14","beta21","sig12x","sig34x","sigksi1","beta212"), n.iter = 2500) 

# store the ending point
fin = Sys.time()

# duration of the run 
duration=fin - deb # approx. 5 seconds
duration

plot(mcmc, trace = TRUE, density = FALSE,ask = dev.interactive())
gelman.diag(mcmc)
summary(mcmc)

plot(mcmc, trace = FALSE, density = TRUE,ask = dev.interactive()) 

#------ same thing w/ model selection

# first list of inits
init1 <- list(p = 0.2,alive=as.matrix(Xinit),sig12x = .1,sig34x = .1,sigksi1 = .1,gamma12 = 1,gamma13 = .5,gamma14 = .5,beta21 = -1,beta212 = -1,w=rep(0,5))
# second list of inits
init2 <- list(p = 0.8,alive=as.matrix(Xinit),sig12x = 1,sig34x = 1,sigksi1 = 1,gamma12 = 1,gamma13 = .5,gamma14 = .5,beta21 = 1,beta212 = 1,w=rep(1,5))
# concatenate list of initial values
inits <- list(init1,init2)

#- tarsus length (x1)
#- phalanx length (x2)
#- wing length (x3)
#- rectrice length (x4)

# store the starting point
deb = Sys.time()

# MCMC simulations 
model <- jags.model(file = "MSSEMMR.bug", data = datax, inits = inits, n.chains = 2) 
mcmc <- coda.samples(model, c("p","gamma12","gamma13","gamma14","beta21","beta212","sig12x","sig34x","sigksi1","w"), n.iter = 2500) 

# store the ending point
fin = Sys.time()

# duration of the run 
duration=fin - deb # approx. 10 seconds
duration


plot(mcmc, trace = TRUE, density = FALSE,ask = dev.interactive())
gelman.diag(mcmc)
summary(mcmc)

plot(mcmc, trace = FALSE, density = TRUE,ask = dev.interactive()) 

# param b
w1 = c(mcmc[[1]][,10],mcmc[[2]][,10])
w2 = c(mcmc[[1]][,11],mcmc[[2]][,11])
mean(w1 & w2)
mean(w2)

# param lambda1, 2 et 3
mean(c(mcmc[[1]][,12],mcmc[[2]][,12]))
mean(c(mcmc[[1]][,13],mcmc[[2]][,13]))
mean(c(mcmc[[1]][,14],mcmc[[2]][,14]))

#------ visualisation

densplot2 <- function (x, show.obs = TRUE, bwf, main = "", ylim, ...) 
{
    xx <- as.matrix(x)
    for (i in 1:nvar(x)) {
        y <- xx[, i, drop = TRUE]
        if (missing(bwf)) 
            bwf <- function(x) {
                x <- x[!is.na(as.vector(x))]
                return(1.06 * min(sd(x), IQR(x)/1.34) * length(x)^-0.2)
            }
        bw <- bwf(y)
        width <- 4 * bw
        if (max(abs(y - floor(y))) == 0 || bw == 0) 
            hist(y, prob = TRUE, main = main, ...)
        else {
            scale <- "open"
            if (max(y) <= 1 && 1 - max(y) < 2 * bw) {
                if (min(y) >= 0 && min(y) < 2 * bw) {
                  scale <- "proportion"
                  y <- c(y, -y, 2 - y)
                }
            }
            else if (min(y) >= 0 && min(y) < 2 * bw) {
                scale <- "positive"
                y <- c(y, -y)
            }
            else scale <- "open"
            dens <- density(y, width = width)
            if (scale == "proportion") {
                dens$y <- 3 * dens$y[dens$x >= 0 & dens$x <= 
                  1]
                dens$x <- dens$x[dens$x >= 0 & dens$x <= 1]
            }
            else if (scale == "positive") {
                dens$y <- 2 * dens$y[dens$x >= 0]
                dens$x <- dens$x[dens$x >= 0]
            }
            if (missing(ylim)) 
                ylim <- c(0, max(dens$y))
            if (is.R()) {
                plot(dens, ylab = "", main = main, type = "l", 
                  xlab = '')
            }
            else {
                plot(dens, ylab = "", main = main, type = "l", 
                  xlab = '')
            }
            if (show.obs) 
                lines(y[1:niter(x)], rep(max(dens$y)/100, niter(x)), 
                  type = "h")
        }
        if (!is.null(varnames(x)) && is.null(list(...)$main)) 
            title(paste("Density of", varnames(x)[i]))
    }
    return(invisible(x))
}

par(mfrow=c(3,3))
densplot2(mcmc[[2]][,'beta21'], show.obs=F,main=expression(gamma[1]))
densplot2(mcmc[[2]][,'beta212'],show.obs=F,main=expression(gamma[2]))
densplot2(mcmc[[1]][,'gamma12'], show.obs=F,main=expression(theta[2]))
densplot2(mcmc[[1]][,'gamma13'], show.obs=F,main=expression(theta[3]))
densplot2(mcmc[[1]][,'gamma14'], show.obs=F,main=expression(theta[4]))
densplot2(mcmc[[1]][,'sig12x'], show.obs=F,main=expression(sigma[1]))
densplot2(mcmc[[1]][,'sig34x'], show.obs=F,main=expression(sigma[2]))
densplot2(mcmc[[1]][,'sigksi1'], show.obs=F,main=expression(sigma[3]))
densplot2(mcmc[[1]][,'p'], show.obs=F,main=expression(p))
