# Code built by Jordan Trinka and Mehdi Maadooliat of Marquette University and Hossein Haghbin of Persian Gulf University
# This script acts as the driver for the MFSSA over different domains algorithm developed by Jordan Trinka, Mehdi Maadooliat, and Hossein Haghbin

# Load Libraries

library("fda") # Functional data analysis package

source("MFSSA.R") # Contains the decomposition stage of MFSSA

source("mfssaplot.R") # Plotting function

source("mfreconstruct.R") # Reconstruction Stage

source("utils.R") # Contains variance/covariance function to be used to measure separability of the multivariate functional time series reconstructions

source("msimulation.R") # Simulation function

# Simulate data

omega = c(1/10,1/5) # frequencies

d = 20 # number of basis functions

basetype = "spline" # basis type

sd = 1 # standard deviation

domlength = 100 # the number of points to evaluate along [0,1]

N = 50 # Length of the multivariate functional time series

noisetype = "FWN" # Type of noise. Currently supports Gaussian white noise (GWN), Functional White Noise (FWN), and Functional Auto Regressive noise (FAR1)

p = 0.1 # Parameter for norm of FAR1 Kernel

l = 0 # horizontal shift parameter

Y = mftssim(omega = omega, N = N, domlength = domlength, l=l, d = d, basetype = basetype, noisetype = noisetype, sd = sd, p = p) # simulate data.

# output of mftssim is a list of length 2 containing the bivariate functional time series

L <- 30 # Lag parameter

#mfssa
  
out <- mfssa(Y,L) # output of this is a list that has the length of number of eigenfunctions. Each list element is another list of length p where p is the number of variables.

plot.mfssa(out,d = 10, type = "values") # plotting. Supports eigenvalues ("vals"), eigenvectors ("vecs"), and paired plots ("paired").
# compplot is used to decide which components are plotted if the user specifies "vecs" or "paired"
plot.mfssa(out,d = 4, type = "vectors",varj = 1) # plotting. Supports eigenvalues ("vals"), eigenvectors ("vecs"), and paired plots ("paired").
# 
plot.mfssa(out,d = 5, type = "pairedu") # plotting. Supports eigenvalues ("vals"), eigenvectors ("vecs"), and paired plots ("paired").
#
plot.mfssa(out,d = 5, type = "pairedv") # plotting. Supports eigenvalues ("vals"), eigenvectors ("vecs"), and paired plots ("paired").
#
plot.mfssa(out,d = 5, type ="meanpaired")
# 
plot.mfssa(out,d = 15, type = "wcor") # plotting. Supports eigenvalues ("vals"), eigenvectors ("vecs"), and paired plots ("paired").
# 
# grouping and reconsjtruction

group=list(c(1,2),c(3,4)) # specifying the reconstructions

reconstructions = mfreconstruct(out = out,group = group) # Reconstruction phase returns a list which has the length of the group list a.k.a. number of reconstructions.
# Each list element contains a multivariate functional time series of length p.


## The rest of this stuff is just me plotting junk
reconstruction_1 <- reconstructions[[1]]
reconstruction_2 <- reconstructions[[2]]

par(mfrow=c(3,2))

plot(Y[[1]], main="Observed Variable 1")

plot(Y[[2]], main="Observed Variable 2")

plot(reconstruction_1[[1]], main="Recon 1 Variable 1")

plot(reconstruction_1[[2]], main="Recon 1 Variable 2")

plot(reconstruction_2[[1]], main ="Recon 2 Variable 1")

plot(reconstruction_2[[2]], main ="Recon 2 Variable 2")

##




