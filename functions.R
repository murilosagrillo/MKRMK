library(extraDistr)
library(AdequacyModel)
library(GenSA)
library(extraDistr)

#------------------------------------------------------------------
#Modified Kumaraswamy (MK) distribution functions#
#------------------------------------------------------------------

#------------------------------------------------------------------
# PDF MK
#------------------------------------------------------------------

dmk <- Vectorize(function(y,alpha,beta,log = FALSE){
  logden <- log((alpha*beta*exp(alpha-alpha/y)*(1-exp(alpha-alpha/y))^(beta-1))/y^2)
  val <- ifelse(log, logden, exp(logden)) 
  return(val)
})

# dmk(0.3,0.3,0.4)

#------------------------------------------------------------------
# Quantile function MK
#------------------------------------------------------------------

qmk <- Vectorize(function(u,alpha,beta){
  val <-  alpha/(alpha-log(1-(1-u)^(1/beta)))
  return(val)
})

# qmk(0.5,1.5,1.3)

#------------------------------------------------------------------
# CDF MK
#------------------------------------------------------------------

pmk <- Vectorize(function(q,alpha,beta,log.p = FALSE){
  cdf <-  1-(1-exp(alpha-alpha/q))^beta
  val <- ifelse(log.p, log(cdf), cdf)
  return(val)
})

# pmk(0.6292858,1.5,1.3)

#------------------------------------------------------------------
# Random number generation MK
#------------------------------------------------------------------

rmk <- function(n,alpha,beta){
  u <- runif(n)
  val <- alpha/(alpha-log(1-(1-u)^(1/beta)))
  return(val)
}

# rmk(100,0.3,0.9)

#------------------------------------------------------------------
# MK Estimation MK
#------------------------------------------------------------------

"mkfit"<-function(y){
  
  y <- as.vector(y)
  n<-length(y)
  set.seed(137)
  
  fdpMK<-function(par){
    alpha=par[1]
    beta=par[2]
    alpha*beta*exp(alpha-alpha/y)*(1-exp(alpha-alpha/y))^(beta-1)*y^(-2)  
  }
  
  #loglike_theta
  l_theta<- function(par){
    sum(log(fdpMK(par)))
  }
  
  gr<-function(par){
    alpha<-par[1] 
    beta<-par[2] 
    c( -(n+(n/alpha)+(beta-1)*exp(alpha)*sum((y-1)/(y*(exp(alpha)-exp(alpha/y))))-sum(1/y)),
       -(n/beta + sum(log(1-exp(alpha-alpha/y))))
    )
  }
  
  ##INICIO genSA para obter os valores de chutes
  fit.sa <- function(y,fdpMK) {
    minusllike <- function(y) -sum(log(fdpMK(c(y[1],y[2]))))
    lower <- c(0.001,0.001) #may need some changes here
    upper <- c(100,100)
    out <- GenSA(lower = lower, upper = upper,
                 fn = minusllike, control=list(verbose=TRUE,max.time=2))
    return(out[c("value","par","counts")])
  }
  
  chute<-fit.sa(y,fdpMK)$par
  ##fim do genSA
  
  res<-optim(chute, l_theta, method="BFGS",#x=x,
             control=list(fnscale=-1, pgtol=1e-20, maxit=200),gr=gr)#, silent=T)
  # if(class(res)=="try-error" || res$conv != 0 ||res$par[2]>.999||res$par[2]<.07) # a classe dos objetos que cont?m o erro,
  
  alpha_emv<-res$par[1]
  beta_semi<--(n/sum(log(1-exp(res$par[1]-res$par[1]/y))))
  
  
  fdpMK_densidade<-function(y){
    alpha_emv*beta_semi*exp(alpha_emv-alpha_emv/y)*(1-exp(alpha_emv-alpha_emv/y))^(beta_semi-1)*y^(-2)  
  }
  
  p<-2
  AIC<-(-2*sum(log( alpha_emv*beta_semi*exp(alpha_emv-alpha_emv/y)*(1-exp(alpha_emv-alpha_emv/y))^(beta_semi-1)*y^(-2)     ))+2*p)
  # AICc<-(-2*sum(log(  alpha_emv*beta_semi*exp(alpha_emv-alpha_emv/y)*(1-exp(alpha_emv-alpha_emv/y))^(beta_semi-1)*y^(-2)    ))+2*p+2*((p*(p+1))/(n-p-1)))
  BIC<-(-2*sum(log(  alpha_emv*beta_semi*exp(alpha_emv-alpha_emv/y)*(1-exp(alpha_emv-alpha_emv/y))^(beta_semi-1)*y^(-2)   ))+p*log(length(y)))
  # HQIC<-(-2*sum(log(   alpha_emv*beta_semi*exp(alpha_emv-alpha_emv/y)*(1-exp(alpha_emv-alpha_emv/y))^(beta_semi-1)*y^(-2)    ))+2*p*log(log(length(y))))
  
  emvs<-c(alpha_emv,beta_semi)
  dados_ordenados = sort(y)
  cdfMK<-function(emvs,dados_ordenados){
    1-(1-exp(alpha_emv-alpha_emv/dados_ordenados))^beta_semi
  }
  
  cdfMK_ks_test<-function(y,par){
    1-(1-exp(alpha_emv-alpha_emv/y))^beta_semi
  }
  
  v = cdfMK(emvs, dados_ordenados)
  s = qnorm(v)
  u = pnorm((s - mean(s))/sqrt(var(s)))
  W_temp <- vector()
  A_temp <- vector()
  for (i in 1:n) {
    W_temp[i] = (u[i] - (2 * i - 1)/(2 * n))^2
    A_temp[i] = (2 * i - 1) * log(u[i]) + (2 * n + 1 - 2 * i) * log(1 - u[i])
  }
  A_2 = -n - mean(A_temp)
  W_2 = sum(W_temp) + 1/(12 * n)
  W_star = W_2 * (1 + 0.5/n)
  A_star = A_2 * (1 + 0.75/n + 2.25/n^2)
  KS = ks.test(x = y, y="cdfMK_ks_test", par = emvs)
  
  g<-c()
  estim <- emvs
  g$estimativas <- estim
  g$convergencia <- res$conv
  # g$valor<-res$value
  data<-y
  hist(data,nclass = 10, freq = F,xlim=c(0.001,1),ylim=c(0.001,5))
  # plot(density(y),xlim=c(0.001,max(y)))
  curve(fdpMK_densidade,xlim=c(0.001,1),ylim=c(0.001,.5),lty=3,add=T)
  # g$chute<-chute
  g$AIC<-AIC
  # g$AICc<-AICc
  g$BIC<-BIC
  # g$HQIC<-HQIC
  g$W_star<- W_star
  g$A_star<- A_star
  g$KS<- KS
  return(g)
}

# mkfit((1-log(rkumar(210,0.1,.3)))^(-1))

#------------------------------------------------------------------
#Reflected Modified Kumaraswamy (RMK) distribution functions#
#------------------------------------------------------------------

#------------------------------------------------------------------
# PDF RMK
#------------------------------------------------------------------

drmk <- Vectorize(function(z,alpha,beta,log = FALSE){
  logden <-log((alpha*beta*exp(z*alpha/(z-1))*(1-exp(z*alpha/(z-1)))^(beta-1))/(z-1)^2)
    #log(alpha)+log(beta)+(z*alpha/(z-1))+(beta-1)*log(1-exp(z*alpha/(z-1)))-log((z-1)^2)
  val <- ifelse(log, logden, exp(logden)) 
  return(val)
})

 # drmk(0.7,0.3,0.4)

#------------------------------------------------------------------
# Quantile function RMK
#------------------------------------------------------------------

qrmk <- Vectorize(function(u,alpha,beta){
  val <-  log(1-u^(1/beta))/(log(1-u^(1/beta))-alpha)
  return(val)
})

 # qrmk(0.5,1.5,1.3)

#------------------------------------------------------------------
# CDF RMK
#------------------------------------------------------------------

prmk <- Vectorize(function(q,alpha,beta,log.p = FALSE){
  cdf <-  (1-exp(q*alpha/(q-1)))^beta
  val <- ifelse(log.p, log(cdf), cdf)
  return(val)
})

 # prmk(0.3707142,1.5,1.3)

#------------------------------------------------------------------
# Random number generation RMK
#------------------------------------------------------------------

rrmk <- function(n,alpha,beta){
  u <- runif(n)
  val <- log(1-u^(1/beta))/(log(1-u^(1/beta))-alpha)
  return(val)
}

 rrmk(100,0.3,0.9)

#------------------------------------------------------------------
# RMK Estimation MK
#------------------------------------------------------------------

 "rmkfit"<-function(z){
   
   z <- as.vector(z)
   n<-length(z)
   set.seed(137)
   
   fdpRMK<-function(par){
     alpha=par[1]
     beta=par[2]
     alpha*beta*exp(z*alpha/(z-1))*(1-exp(z*alpha/(z-1)))^(beta-1)*(z-1)^(-2)
   }
   
   #loglike_theta
   l_theta<- function(par){
     sum(log(fdpRMK(par)))
   }
   
   gr<-function(par){
     alpha<-par[1]
     beta<-par[2]
     y=1-z
     c( -(n+(n/alpha)+(beta-1)*exp(alpha)*sum((y-1)/(y*(exp(alpha)-exp(alpha/y))))-sum(1/y)),
        -(n/beta + sum(log(1-exp(alpha-alpha/y))))
     )
   }
   
   ##INICIO genSA para obter os valores de chutes
   fit.sa <- function(z,fdpRMK) {
     minusllike <- function(z) -sum(log(fdpRMK(c(z[1],z[2]))))
     lower <- c(0.001,0.001) #may need some changes here
     upper <- c(100,100)
     out <- GenSA(lower = lower, upper = upper,
                  fn = minusllike, control=list(verbose=TRUE,max.time=2))
     return(out[c("value","par","counts")])
   }
   
   chute<-fit.sa(z,fdpRMK)$par
   ##fim do genSA
   
   res<-optim(chute, l_theta, method="BFGS",#x=x,
              control=list(fnscale=-1, pgtol=1e-20, maxit=200),gr=gr)#, silent=T)
   # if(class(res)=="try-error" || res$conv != 0 ||res$par[2]>.999||res$par[2]<.07) # a classe dos objetos que cont?m o erro,
   
   alpha_emv<-res$par[1]
   beta_semi<--(n/sum(log(1-exp(res$par[1]-res$par[1]/(1-z)))))
   
   
   fdpRMK_densidade<-function(z){
     alpha_emv*beta_semi*exp(z* alpha_emv/(z-1))*(1-exp(z* alpha_emv/(z-1)))^(beta_semi-1)*(z-1)^(-2)
   }
   
   p<-2
   AIC<-(-2*sum(log(  alpha_emv*beta_semi*exp(z* alpha_emv/(z-1))*(1-exp(z* alpha_emv/(z-1)))^(beta_semi-1)*(z-1)^(-2)     ))+2*p)
   # AICc<-(-2*sum(log(   alpha_emv*beta_semi*exp(z* alpha_emv/(z-1))*(1-exp(z* alpha_emv/(z-1)))^(beta_semi-1)*(z-1)^(-2)   ))+2*p+2*((p*(p+1))/(n-p-1)))
   BIC<-(-2*sum(log(   alpha_emv*beta_semi*exp(z* alpha_emv/(z-1))*(1-exp(z* alpha_emv/(z-1)))^(beta_semi-1)*(z-1)^(-2)   ))+p*log(length(z)))
   # HQIC<-(-2*sum(log(    alpha_emv*beta_semi*exp(z* alpha_emv/(z-1))*(1-exp(z* alpha_emv/(z-1)))^(beta_semi-1)*(z-1)^(-2)    ))+2*p*log(log(length(z))))
   
   emvs<-c(alpha_emv,beta_semi)
   dados_ordenados = sort(z)
   
   cdfRMK<-function(emvs,dados_ordenados){
     (1-exp(dados_ordenados*alpha_emv/(dados_ordenados-1)))^beta_semi
   }
   
   cdfRMK_ks_test<-function(z,par){
     (1-exp(z*alpha_emv/(z-1)))^beta_semi
   }
   
   v = cdfRMK(emvs, dados_ordenados)
   s = qnorm(v)
   u = pnorm((s - mean(s))/sqrt(var(s)))
   W_temp <- vector()
   A_temp <- vector()
   for (i in 1:n) {
     W_temp[i] = (u[i] - (2 * i - 1)/(2 * n))^2
     A_temp[i] = (2 * i - 1) * log(u[i]) + (2 * n + 1 - 2 * i) * log(1 - u[i])
   }
   A_2 = -n - mean(A_temp)
   W_2 = sum(W_temp) + 1/(12 * n)
   W_star = W_2 * (1 + 0.5/n)
   A_star = A_2 * (1 + 0.75/n + 2.25/n^2)
   KS = ks.test(x = z, y="cdfRMK_ks_test", par = emvs)
   
   g<-c()
   estim <- emvs
   g$estimativas <- estim
   g$convergencia <- res$conv
   # g$valor<-res$value
   data<-z
   hist(data,nclass = 10, freq = F,xlim=c(0.001,1),ylim=c(0.001,5))
   # plot(density(y),xlim=c(0.001,max(y)))
   curve(fdpRMK_densidade,xlim=c(0.001,1),ylim=c(0.001,.5),lty=3,add=T)
   # g$chute<-chute
   g$AIC<-AIC
   # g$AICc<-AICc
   g$BIC<-BIC
   # g$HQIC<-HQIC
   g$W_star<- W_star
   g$A_star<- A_star
   g$KS<- KS
   return(g)
 }

# rmkfit(1-(1-log(rkumar(210,0.1,.3)))^(-1))

 