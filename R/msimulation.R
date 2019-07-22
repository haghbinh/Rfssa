## Simulation



mftssim<-function(omega = c(1/10, 1/5), N = 50, domlength = 100, l = 0, d = 20, basetype = "spline", noisetype = "GWN", sd = 0.1, p = (1/2)){
  
  library("fda")
  
  s <- seq(0, 1, length.out = domlength)
  
  y_1 <- matrix(nrow = domlength, ncol = N, data = NA)
  
  y_2 <- matrix(nrow = domlength, ncol = N, data = NA)
  
  
  for(t in 1:N){
    
    y_1[,t] <- exp(s^(2))*cos(2*pi*omega[1]*(t+l))-cos(4*pi*s)*sin(2*pi*omega[2]*t)+exp(1-s^(2))*cos(2*pi*omega[1]*t)+sin(pi*s)*sin(2*pi*omega[2]*t)
    
    y_2[,t] <- exp(s^(2))*sin(2*pi*omega[1]*t)+cos(4*pi*s)*cos(2*pi*omega[2]*(t-l))
    
    
  }
  
  y=list(y_1,y_2)
  
  Y=list()
  
  if(noisetype == "GWN"){
    
    for(j in 1:2){
      
      if(basetype == "spline"){
        
        Y[[j]] = smooth.basis(s, y[[j]]+matrix(nrow = domlength, ncol = N, data = rnorm(n = domlength*N, mean = 0, sd = sd)), create.bspline.basis(rangeval=c(0,1),nbasis=d))$fd
      
      }else if(basetype == "fourier"){
        
        Y[[j]] = smooth.basis(s, y[[j]]+matrix(nrow = domlength, ncol = N, data = rnorm(n = domlength*N, mean = 0, sd = sd)), create.fourier.basis(rangeval=c(0,1),nbasis=d))$fd
        
      }
      
    }
    
    
  }else if(noisetype == "FWN"){
    
    for(j in 1:2){
      
      if(basetype == "spline"){
        
        Y[[j]] = smooth.basis(s, y[[j]]+eval.fd(s,fd(matrix(rnorm(N*d, mean = 0, sd = sd),ncol = N), create.bspline.basis(c(0,1), d))), create.bspline.basis(rangeval=c(0,1),nbasis=d))$fd
      
      }else if(basetype == "fourier"){
        
        Y[[j]] = smooth.basis(s, y[[j]]+eval.fd(s,fd(matrix(rnorm(N*d, mean = 0, sd = sd),ncol = N), create.fourier.basis(c(0,1), d))), create.fourier.basis(rangeval=c(0,1),nbasis=d))$fd
        
      }
      
    }
    
    
  }else if(noisetype == "FAR1"){
    
    
    rfar <- function(N,norm,psi,Eps,basis){
      # Create Corresponding matrix of an integral (kernel) operator
      # Get an kernel function k(s,t) corresponds to an integral operator and 
      # a basis system incuding d basis functions, return corresponding d*d matrix of the operator with respect
      # to basis system.
      OpsMat <- function(kernel,basis){
        u <- seq(0,1,by=0.01)
        n <- length(u)
        K_mat <- outer(u,u,FUN = kernel)
        K_t <- smooth.basis(u,K_mat,basis)$fd # the kernel function convert
        # to fd object w.r.t the first argument. K_t is an objcect of n fd.
        A <- inprod(K_t,basis) #An n*d matrix.
        K <- smooth.basis(u,A,basis)$fd # d fd object.
        B <- inprod(K,basis) # An d*d matrix.
        return(B) # return to OpsMat.
      }
      Psi_mat0 <- OpsMat(psi,basis)
      Gram <- inprod(basis,basis)
      Psi_mat <- solve(Gram)%*%Psi_mat0
      E <- Eps$coefs
      X <- E;
      for(i in 2:N) X[,i] <- Psi_mat%*%X[,i-1]+E[,i]
      X_fd <- fd(X,basis)
      return(X_fd) # return to rfar function.
    }
    gamma0 <- function(norm){
      f <- function(x){
        g <- function(y) psi0(x,y)^2
        return(integrate(g,0,1)$value) #return into f.
      }
      f <- Vectorize(f)
      A <- integrate(f,0,1)$value
      return(norm/A) #return into gamma.
    }
    psi0 <- function(x,y) 2-(2*x-1)^2-(2*y-1)^2
    
    set.seed(domlength*N*(1/10)*sd);
    Z <- matrix(rnorm(N*domlength, 0, sd), nrow=domlength);
    k0 <- gamma0(p)
    psi <- function(x,y) k0*psi0(x,y)
    Z[1, ] <- 0; noise_1 <- apply(Z, 2, cumsum)
    for(j in 1:2){
      if (basetype == "spline"){ 
        
        basis.Z <- create.bspline.basis(c(0, 1), d)
        Eps <- smooth.basis(s, noise_1, basis.Z)$fd
        basis.noise <- rfar(N,p,psi,Eps,basis.Z)
        Y[[j]] <- smooth.basis(s, y[[j]]+eval.fd(s, basis.noise), create.bspline.basis(rangeval=c(0,1),nbasis=d))$fd
        
      }
      else if(basetype == "fourier"){
        
        basis.Z <- create.fourier.basis(c(0, 1), d)
        Eps <- smooth.basis(s, noise, basis.Z)$fd
        basis.noise <- rfar(N,p,psi,Eps,basis.Z)
        Y[[j]] <- smooth.basis(s, y[[j]]+eval.fd(s, basis.noise), create.fourier.basis(rangeval=c(0,1),nbasis=d))$fd
        
      }
    
    }
    
    
    
  }
  
  par(mfrow=c(1,2))
  
  plot(Y[[1]], main = "Simulated Variable 1")
  
  plot(Y[[2]], main = "Simulated Variable 2")
  
  par(mfrow=c(1,1))
  
  return(Y)
  
}


