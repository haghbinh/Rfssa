# Code built by Jordan Trinka and Mehdi Maadooliat of Marquette University and Hossein Haghbin of Persian Gulf University

plot.mfssa <-function(out = NULL, d = length(out$values), type = "vecs", varj = 1){
  
  p <- length(out$U[[1]])
  
  finddim <- function(d){
    
    if(d==1){
      
      return(c(1,1))
      
    }else if(d==2){
    
      return(c(1,2))
    
    }else if(d==3){
      
      return(c(1,3))
      
    }else if(d>=4 && d<=8){
    
      return(c(ceiling(d/2),2))
    
    
    }else if(d>=9 && d<=12){
    
      return(c(ceiling(d/3),3))  
      
    }else if(d>=13 && d<=16){
      
      return(c(ceiling(d/4),4))
      
    }else if(d>=17 && d<=20)
      
      return(c(ceiling(d/5),5))
  
  }
  
  # plot eigenvalues
  if(type == "values"){
    
    plot(out$values[1:d],log = "y", type = "o", lwd = 2L, col = "dodgerblue3", pch = 19L, cex = 0.8, main = "Eigenvalues", ylab = "Log of Value")
    #plot eigenvejtores
  }else if(type == "pairedu"){
    
    
    t_v <- sum(out$values)
    
    if(d==1){
      
      d=2
      
    }
    
    dimensions <- finddim(d-1)
    
    par(mfrow=c(dimensions[1],dimensions[2]))
      
    for(i in 1:(d-1)){
        
      plot(out$coef[,i],out$coef[,i+1],type = "l", main = paste(as.character(i),"(",as.character(round(100*out$values[i]/t_v,digits = 2)),"%)","vs.",as.character(i+1),"(",as.character(round(100*out$values[i+1]/t_v,digits = 2)),"%)"),xlab = "", ylab = "",xaxt = "n", yaxt = "n", col = "dodgerblue3")
        
    }
    
    par(mfrow=c(1,1))
      
  
    
  }else if(type == "vectors"){
    
     t_v <- sum(out$values)  
    
     dimensions <- finddim(d)
     
     par(mfrow=c(dimensions[1],dimensions[2]))
      
       for(i in 1:d){
        
         mult_pc <- out$U[[i]]
         plot(mult_pc[[varj]],main = paste(as.character(i),"(",as.character(round(100*out$values[i]/t_v,digits = 2)),"%)"),xlab = "", ylab = "",xaxt = "n", yaxt = "n")
          
      }
       
      par(mfrow=c(1,1))
      
    
  }else if (type == "wcor"){
    library(Rfssa)
    source("mfwcor.R")
    cor_mat <- mfwcor(out,d)
    wplot(cor_mat)
    
  }else if(type == "meanpaired"){
    
    library("fda")
    
    t_v <- sum(out$values)
    
    if(d==1){
      
      d=2
      
    }
    
    dimensions <- finddim(d-1)
    
    par(mfrow=c(dimensions[1],dimensions[2]))
    
    for(i in 1:(d-1)){
      
      p_c_i = out$U[[i]]
      
      p_c_i_p = out$U[[i+1]]
      
      mat_i = matrix(ncol=out$L)
      
      mat_i_p = matrix(ncol=out$L)
      
      for(j in 1:p){
      
        p_c_i_j = p_c_i[[j]]
      
        p_c_i_p_j = p_c_i_p[[j]]
      
        e_i_j = eval.fd(evalarg = seq(p_c_i_j$basis$rangeval[1],p_c_i_j$basis$rangeval[2],length.out = 100) ,fdobj = p_c_i_j)
      
        e_i_p_j = eval.fd(evalarg = seq(p_c_i_p_j$basis$rangeval[1],p_c_i_p_j$basis$rangeval[2],length.out = 100) ,fdobj = p_c_i_p_j)
        
        mat_i=rbind(mat_i,e_i_j)
        
        mat_i_p=rbind(mat_i_p,e_i_p_j)
      
      }
      
      mat_i=mat_i[2:nrow(mat_i),]
      
      mat_i_p=mat_i_p[2:nrow(mat_i_p),]
      
      plot(colMeans(mat_i),colMeans(mat_i_p),type = "l",lwd = 3, main = paste(as.character(i),"(",as.character(round(100*out$values[i]/t_v,digits = 2)),"%)","vs.",as.character(i+1),"(",as.character(round(100*out$values[i+1]/t_v,digits = 2)),"%)"),xlab = "", ylab = "",xaxt = "n", yaxt = "n", col = "dodgerblue3")
      
    }
    
    par(mfrow=c(1,1))
    
    
  }else if(type == "pairedv"){
    
    library("fda")
    
    t_v <- sum(out$values)
    
    if(d==1){
      
      d=2
      
    }
    
    dimensions <- finddim(d-1)
    
    par(mfrow=c(dimensions[1],dimensions[2]))
    
    for(i in 1:(d-1)){
      
      
      plot(out$V[,i],out$V[,i+1],type = "l",lwd = 3, main = paste(as.character(i),"(",as.character(round(100*out$values[i]/t_v,digits = 2)),"%)","vs.",as.character(i+1),"(",as.character(round(100*out$values[i+1]/t_v,digits = 2)),"%)"),xlab = "", ylab = "",xaxt = "n", yaxt = "n", col = "dodgerblue3")
      
    }
    
    par(mfrow=c(1,1))
    
    
  }
}
