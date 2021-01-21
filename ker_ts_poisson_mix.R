rm(list = ls())
setwd("/Users/yangwang/Documents/paper_code")
library(data.table)
library(intervals)
f.mu<-function(x){8*sqrt(x)+4}
#f.mu<-function(x){x^2+4}
f.phi<-function(x){x^2/9}
#f.phi<-function(x){sin(pi*x/3)}

start<-Sys.time()
#create a dataset with mixed Poisson process
data.gen.mix<-function(
  nsub,                       #Sample size
  max.obs,                    #Maximal number of observation for each individual
  irate,                      #Rate for interarrival time
  len                         #length of study
)
{
  Z1<-round(runif(nsub,0,3),4)         #define covariate
  
  drnd <- function(x, p, n) {
    z <- NULL
    ps <- cumsum(p)
    r <- runif(n)
    for (i in 1:n) z <- c(z, x[which(r[i] <= ps)[1]])
    return(z)
  }
  alpha<-c(-0.4,0,0.4)
  prob<-c(0.25,0.5,0.25)
  reff<-drnd(alpha,prob,nsub)  #####generate frailty or random effect on the intensity
  
  NumObs<-rep(nsub,0)
  Mobs<-rep(nsub,0)
  rexpf<-function(n,irate){
    u<-runif(n,0.001,0.999)     
    x<--log(u)/irate          
    return(x)
  }                        #####construct the expoentional distribution from uniform distribution
  
  x<-round(rexpf(max.obs,irate),4)     #generate time for the first subject
  while(x[1]>=len)                     #if the first time larger than the length of study,generate time again
    x<-round(rexpf(max.obs,irate),4)
  for(j in 1:max.obs){
    if(sum(x[1:j])<len)
      k<-j
  }                           #choose the cumulative time less than the length of study
  indtimeObs<-cumsum(x[1:k])  #interarrival time
  NumObs[1]<-length(indtimeObs) #the number of observations for the first subject
  Mobs[1]<-indtimeObs[length(indtimeObs)] #the maximum observation times for the first subject
  
  subj<-rep(1,NumObs[1])      #identification for the first subject
  indcount<-rep(0,NumObs[1])  #cumulative counts
  lambda<-rep(0,NumObs[1])
  lambda[1]<-(f.mu(indtimeObs[1])+reff[1])*exp(f.phi(Z1[1]))
  indcount[1]<-rpois(1,lambda[1])  #cumulative count when observation only one
  
  if(NumObs[1]>1){
    for(j in 2:NumObs[1]){
      lambda[j]<-(f.mu(indtimeObs[j])+reff[1])*exp(f.phi(Z1[1]))
      indcount[j]<-indcount[j-1]+rpois(1,lambda[j]-lambda[j-1])
    }     #end of j
  }
  
  lambda.t<-lambda
  timeObs<-indtimeObs          #define observation time
  count<-indcount              #define observation count
  covar1<-rep(Z1[1],NumObs[1]) #covariate for subject 1
  reffect<-rep(reff[1],NumObs[1])
  Mobstime<-rep(Mobs[1],NumObs[1])     
  
  #To expand the database by adding the subsequence observations
  for(i in 2:nsub){
    x<-round(rexpf(max.obs,irate),4)  #generate follow-up time for the ith subject
    while(x[1]>=len)
      x<-round(rexpf(max.obs,irate),4)#if the first time larger than the length of study, generate time again
    for(j in 1:max.obs){
      if(sum(x[1:j])<len)
        k<-j
    }                                 #choose the cumulative time less than the length of study
    indtimeObs<-cumsum(x[1:k])        #interarrival time
    NumObs[i]<-length(indtimeObs)     #the number of observations for the ith subject
    Mobs[i]<-indtimeObs[length(indtimeObs)] #the maximum observation times for the ith subject
    
    subj<-c(subj,rep(i,NumObs[i]))    #generate identification from i to nsub
    indcount<-rep(0,NumObs[i])
    lambda<-rep(0,NumObs[i])
    lambda[1]<-(f.mu(indtimeObs[1])+reff[i])*exp(f.phi(Z1[i]))
    indcount[1]<-rpois(1,lambda[1])
    if(NumObs[i]>1)
      for(j in 2:NumObs[i]){
        lambda[j]<-(f.mu(indtimeObs[j])+reff[i])*exp(f.phi(Z1[i]))
        indcount[j]<-indcount[j-1]+rpois(1,lambda[j]-lambda[j-1])
      }
    lambda.t<-c(lambda.t,lambda)
    timeObs<-c(timeObs,indtimeObs)      #combine all the observation
    count<-c(count,indcount)            #combine all the count
    covar1<-c(covar1,rep(Z1[i],NumObs[i]))#combine all the covariate
    reffect<-c(reffect,rep(reff[i],NumObs[i]))
    Mobstime<-c(Mobstime,rep(Mobs[i],NumObs[i]))#combine all the maximum observation times
  }
  cbind(subj,timeObs,count,covar1,Mobstime,lambda.t,reffect)
}   #end of generate data

#create a dataset with Poisson process
data.gen<-function(
  nsub,                       #Sample size
  max.obs,                    #Maximal number of observation for each individual
  irate,                      #Rate for interarrival time
  len                         #length of study
)
{
  Z1<-round(runif(nsub,0,3),4)         #define covariate 
  #To creat a database starting from the first subject
  NumObs<-rep(nsub,0)
  Mobs<-rep(nsub,0)
  Mcoun<-rep(nsub,0)
  rexpf<-function(n,irate){
    u<-runif(n,0.001,0.999)     
    x<--log(u)/irate          
    return(x)
  }    #####construct the expoentional distribution from uniform distribution
  
  x<-round(rexpf(max.obs,irate),4)   #generate time for the first subject
  while(x[1]>=len)            #if the first time larger than the length of study, generate time again
    x<-round(rexpf(max.obs,irate),4)
  for(j in 1:max.obs){
    if(sum(x[1:j])<len)
      k<-j
  }                           #choose the cumulative time less than the length of study
  indtimeObs<-cumsum(x[1:k])  #interarrival time
  NumObs[1]<-length(indtimeObs) #the number of observations for the first subject
  Mobs[1]<-indtimeObs[length(indtimeObs)] #the maximum observation times for the first subject
    
  subj<-rep(1,NumObs[1])      #identification for the first subject
  indcount<-rep(0,NumObs[1])  #cumulative counts
  lambda<-rep(0,NumObs[1])
  lambda[1]<-f.mu(indtimeObs[1])*exp(f.phi(Z1[1]))
  indcount[1]<-rpois(1,lambda[1])  #cumulative count when observation only one
  
  if(NumObs[1]>1){
    for(j in 2:NumObs[1]){
      lambda[j]<-f.mu(indtimeObs[j])*exp(f.phi(Z1[1]))
      indcount[j]<-indcount[j-1]+rpois(1,lambda[j]-lambda[j-1])
    }     #end of j
  }
  Mcoun[1]<-indcount[length(indtimeObs)]
  lambda.t<-lambda
  timeObs<-indtimeObs          #define observation time
  count<-indcount              #define observation count
  covar1<-rep(Z1[1],NumObs[1]) #covariate1 for subject 1
  Mobstime<-rep(Mobs[1],NumObs[1])     #the maximum observation time for subject1
  Mcount<-rep(Mcoun[1],NumObs[1])
  #To expand the database by adding the subsequence observations
  for(i in 2:nsub){
    x<-round(rexpf(max.obs,irate),4)  #generate time for the ith subject
    while(x[1]>=len)
      x<-round(rexpf(max.obs,irate),4)#if the first time larger than the length of study,generate time again
    for(j in 1:max.obs){
      if(sum(x[1:j])<len)
        k<-j
    }                                 #choose the cumulative time less than the length of study
    indtimeObs<-cumsum(x[1:k])        #interarrival time
    NumObs[i]<-length(indtimeObs)     #the number of observations for the ith subject
    Mobs[i]<-indtimeObs[length(indtimeObs)] #the maximum observation times
    
    subj<-c(subj,rep(i,NumObs[i]))    #generate identification from i to nsub
    indcount<-rep(0,NumObs[i])
    lambda<-rep(0,NumObs[i])
    lambda[1]<-f.mu(indtimeObs[1])*exp(f.phi(Z1[i]))
    indcount[1]<-rpois(1,lambda[1])
    if(NumObs[i]>1)
      for(j in 2:NumObs[i]){
        lambda[j]<-f.mu(indtimeObs[j])*exp(f.phi(Z1[i]))
        indcount[j]<-indcount[j-1]+rpois(1,lambda[j]-lambda[j-1])
      }
    Mcoun[i]<-indcount[length(indtimeObs)]
    lambda.t<-c(lambda.t,lambda)
    timeObs<-c(timeObs,indtimeObs)      #combine all the observation
    count<-c(count,indcount)            #combine all the count
    covar1<-c(covar1,rep(Z1[i],NumObs[i]))#combine all the covariate1
    Mobstime<-c(Mobstime,rep(Mobs[i],NumObs[i]))#combine all the maximum observation times
    Mcount<-c(Mcount,rep(Mcoun[i],NumObs[i])) #combine all the last event count
  }
  cbind(subj,timeObs,count,covar1,Mobstime,Mcount,lambda.t)
}   #end of generate data

###########simulation ############

ker<-function(x,h){ re<-0.75*(1-(x/h)^2)*(abs(x/h)<1)/h }  #epanechnikov kernel function

beta.fun<-function(data,beta,z,h){
  obstime<-data[,2]
  count<-data[,3]
  covar<-data[,4]
  mtime<-data[,5]
  fun<-function(beta){
    S<-sapply(obstime,function(x){
      a0<-sum(exp(beta[1]*(covar-z)+beta[2]*(covar-z)^2)*ker(covar-z,h)*(mtime>=x))
      a1<-sum(exp(beta[1]*(covar-z)+beta[2]*(covar-z)^2)*ker(covar-z,h)*(covar-z)*(mtime>=x))
      a2<-sum(exp(beta[1]*(covar-z)+beta[2]*(covar-z)^2)*ker(covar-z,h)*(covar-z)^2*(mtime>=x))
      a3<-sum(exp(beta[1]*(covar-z)+beta[2]*(covar-z)^2)*ker(covar-z,h)*(covar-z)^3*(mtime>=x))
      a4<-sum(exp(beta[1]*(covar-z)+beta[2]*(covar-z)^2)*ker(covar-z,h)*(covar-z)^4*(mtime>=x))
      if(a0==0){re<-c(0,0,0,0,0)}else{re<-c(a1/a0,a2/a0,a1^2/a0^2-a2/a0,a1*a2/a0^2-a3/a0,a2^2/a0^2-a4/a0)}
    })
    #for each observed number T(il),calculate the at risk number of score function and hessian matrix
    sum11<-sum(ker(covar-z,h)*count*(covar-z-S[1,]))
    sum12<-sum(ker(covar-z,h)*count*((covar-z)^2-S[2,]))
    df<-c(sum11,sum12)
    sum21<-sum(ker(covar-z,h)*count*S[3,])
    sum22<-sum(ker(covar-z,h)*count*S[4,])
    sum23<-sum(ker(covar-z,h)*count*S[5,])
    jac<-matrix(c(sum21,sum22,sum22,sum23),nr=2)
    list(df=df,jac=jac)
  }            #score function df and hessian matrix jac
  
  Newton<-function(fun,beta,eps=1e-5){
    k<-0
    repeat{
      k<-k+1
      beta1<-beta
      obj<-fun(beta)
      beta<-beta-solve(obj$jac,obj$df)
      #print(obj)
      if((beta-beta1)%*%(beta-beta1)<eps){
        print(k)
        return(beta)
        break
      }
    }
  }           
  return(Newton(fun,beta))
}     #newton-raphson iteration algorithm of the proposed method

beta.fun.ts<-function(data,beta,z,h,Ft){
  obstime<-data[,2]
  count<-data[,3]
  covar<-data[,4]
  
  fun<-function(beta){
    sum11<-sum(ker(covar-z,h)*(count/Ft-exp(beta[1]+beta[2]*(covar-z))))
    sum12<-sum(ker(covar-z,h)*(covar-z)*(count/Ft-exp(beta[1]+beta[2]*(covar-z))))
    df<-c(sum11,sum12)
    sum21<--sum(ker(covar-z,h)*exp(beta[1]+beta[2]*(covar-z)))
    sum22<--sum(ker(covar-z,h)*(covar-z)*exp(beta[1]+beta[2]*(covar-z)))
    sum23<--sum(ker(covar-z,h)*(covar-z)^2*exp(beta[1]+beta[2]*(covar-z)))
    jac<-matrix(c(sum21,sum22,sum22,sum23),nr=2)
    list(df=df,jac=jac)
  }            #score function df and hessian matrix jac
  
  Newton<-function(fun,beta,eps=1e-4){
    k<-0
    repeat{
      k<-k+1
      beta1<-beta
      obj<-fun(beta)
      beta<-beta-solve(obj$jac,obj$df)
      #print(obj)
      if((beta-beta1)%*%(beta-beta1)<eps){
        print(k)
        return(beta)
        break
      }
    }
  }           
  return(Newton(fun,beta))
}#newton-raphson iteration algorithm of the two-stage approach
shape.fun<-function(data,tim){
  Insub<-data[,1]
  obstime<-data[,2]
  count<-data[,3]
  covar<-data[,4]
  mtime<-data[,5]
  mcount<-data[,6]
  unmtime<-unname(tapply(mtime,Insub,unique))
  unmcount<-unname(tapply(mcount,Insub,unique))
  ordtime<-unique(sort(obstime))  
  nordtime<-c(ordtime,len)
  Matime<-as.matrix(cbind(c(0,ordtime),c(ordtime,len)))
  Intime<-Intervals(Matime,closed=c(TRUE,TRUE),type="R")
  
  count0<-unname(unlist(tapply(count,Insub,function(x){c(0,x[-length(x)])})))
  decount<-count-count0
  obstime0<-unname(unlist(tapply(obstime,Insub,function(x){c(0,x[-length(x)])})))
  Maobstime<-as.matrix(cbind(obstime0,obstime))
  Inobstime<-Intervals(Maobstime,closed=c(TRUE,TRUE),type="R")
  
  Alist<-interval_included(Inobstime,Intime)
  
  Amatrix<-matrix(0,nrow = dim(Inobstime)[1],ncol = dim(Intime)[1])
  Bmatrix<-matrix(0,nrow = nsub,ncol = dim(Intime)[1])
  for (i in 1:dim(Inobstime)[1]) {
    Amatrix[i,Alist[[i]]]<-1
  }
  for (i in 1:dim(Intime)[1]) {
    Bmatrix[,i]<-1*(nordtime[i]<=unmtime)
  }
  p<-rep(1/dim(Intime)[1],dim(Intime)[1])
  for (step in 1:10000) {
    ##E-step
    pa<-Amatrix%*%p
    pb<-Bmatrix%*%p
    D<-sapply(1:dim(Intime)[1],function(x){
      re<-sum(decount*(Amatrix[,x]*p[x]/pa))+sum(unmcount*c((1-Bmatrix[,x])*p[x]/pb))
    })
    oldp<-p
    ##M-step
    p<-D/sum(D)
    epslio<-1e-4
    if(sum(abs(oldp-p))<epslio) break
    #print(step)
  }  
  phat<-p
  Fhat<-sapply(tim,function(x){
    sum(phat[which(nordtime<=x)])
  })
  Flast<-c(Bmatrix%*%phat)
  return(list(Flast,Fhat))
}#the estimates for shape function
cvscore<-function(data,beta1,beta2,z,h){
  obstime<-data[,2]
  count<-data[,3]
  covar<-data[,4]
  mtime<-data[,5]
  S<-sapply(obstime,function(x){
    a0<-sum(exp(beta1*(covar-z)+beta2*(covar-z)^2)*ker(covar-z,h)*(mtime>=x))
    a1<-sum(exp(beta1*(covar-z)+beta2*(covar-z)^2)*ker(covar-z,h)*(covar-z)*(mtime>=x))
    a2<-sum(exp(beta1*(covar-z)+beta2*(covar-z)^2)*ker(covar-z,h)*(covar-z)^2*(mtime>=x))
    a3<-sum(exp(beta1*(covar-z)+beta2*(covar-z)^2)*ker(covar-z,h)*(covar-z)^3*(mtime>=x))
    a4<-sum(exp(beta1*(covar-z)+beta2*(covar-z)^2)*ker(covar-z,h)*(covar-z)^4*(mtime>=x))
    if(a0==0){re<-c(0,0,0,0,0,0)}else{re<-c(a1/a0,a2/a0,a1^2/a0^2-a2/a0,a1*a2/a0^2-a3/a0,a2^2/a0^2-a4/a0,log(a0/nsub))}
  })
 
  sum11<-sum(ker(covar-z,h)^2*count^2*(covar-z-S[1,])^2)
  sum12<-sum(ker(covar-z,h)^2*count^2*((covar-z)^2-S[2,])*(covar-z-S[1,]))
  sum13<-sum(ker(covar-z,h)^2*count^2*((covar-z)^2-S[2,])^2)
  dfm<-matrix(c(sum11,sum12,sum12,sum13),nr=2)
  sum21<-sum(ker(covar-z,h)*count*S[3,])
  sum22<-sum(ker(covar-z,h)*count*S[4,])
  sum23<-sum(ker(covar-z,h)*count*S[5,])
  jac<-matrix(c(sum21,sum22,sum22,sum23),nr=2)
  
  cvl1<-sum(diag(solve(jac)%*%dfm))
  cvl2<-sum(ker(covar-z,h)*count*(beta1*(covar-z)+beta2*(covar-z)^2-S[6,]))
  return(cvl1+cvl2)
}   #the cross-validation score function of the proposed method
cvscore.ts<-function(data,beta1,beta2,z,h,Ft){
  obstime<-data[,2]
  count<-data[,3]
  covar<-data[,4]
  sum11<-sum(ker(covar-z,h)^2*(count/Ft-exp(beta1+beta2*(covar-z)))^2)
  sum12<-sum(ker(covar-z,h)^2*(covar-z)*(count/Ft-exp(beta1+beta2*(covar-z)))^2)
  sum13<-sum(ker(covar-z,h)^2*(covar-z)^2*(count/Ft-exp(beta1+beta2*(covar-z)))^2)
  dfm<-matrix(c(sum11,sum12,sum12,sum13),nr=2)
  sum21<--sum(ker(covar-z,h)*exp(beta1+beta2*(covar-z)))
  sum22<--sum(ker(covar-z,h)*(covar-z)*exp(beta1+beta2*(covar-z)))
  sum23<--sum(ker(covar-z,h)*(covar-z)^2*exp(beta1+beta2*(covar-z)))
  jac<-matrix(c(sum21,sum22,sum22,sum23),nr=2)
  
  cvl1<-sum(diag(solve(jac)%*%dfm))
  cvl2<-sum(ker(covar-z,h)*((beta1+beta2*(covar-z))*(count/Ft)-exp(beta1+beta2*(covar-z))))
  return(cvl1+cvl2)
}#the cross-validation score function of the two-stage approach
choose.hb<-function(data,z){
  hvalue<-seq(0.5,0.9,by=0.2)
  cvh<-function(z,h){
    beta0<-beta.fun(data,beta,z,h)
    cv<-cvscore(data,beta0[1],beta0[2],z,h)
    return(c(beta0,cv))
  }
  cv.h<-sapply(hvalue,function(x){
    re<-cvh(z,h=x)
  })
  opt.h<-which.max(cv.h[3,])
  hopt<-hvalue[opt.h]
  betaopt0<-cv.h[1,][opt.h]
  betaopt1<-cv.h[2,][opt.h]
  return(c(hopt,betaopt0,betaopt1))
} #the cross-validation approach for bandwidth of the proposed method
choose.hb.ts<-function(data,z,Ft){
  hvalue<-seq(0.5,0.9,by=0.2)
  cvh<-function(z,h){
    beta0<-beta.fun.ts(data,beta,z,h,Ft)
    cv<-cvscore.ts(data,beta0[1],beta0[2],z,h,Ft)
    return(c(beta0,cv))
  }
  cv.h<-sapply(hvalue,function(x){
    re<-cvh(z,h=x)
  })
  opt.h<-which.max(cv.h[3,])
  hopt<-hvalue[opt.h]
  betaopt0<-cv.h[1,][opt.h]
  betaopt1<-cv.h[2,][opt.h]
  return(c(hopt,betaopt0,betaopt1))
}#the cross-validation approach for bandwidth of the two-stage approach
############simulation for 500 replications#############
simulation.est<-function(nr,data){
  All<-list()
  for(i in 1:nr){
    set.seed(100*i+i)
    print(100*i+i)
    obstime<-data[,2]
    count<-data[,3]
    covar<-data[,4]
    mtime<-data[,5]
    zp<-seq(0.04,2.98,by=0.06)
    bhopt<-sapply(zp,function(x){
      re<-choose.hb(data,z=x)
    })
    hopt<-bhopt[1,]
    hbeta1<-bhopt[2,]
    hbeta2<-bhopt[3,]
    Var<-rep()
    hbeta0<-rep(0,50) 
    hbeta0[1]<-hbeta1[1]*zp[1]  
    for(j in 2:50){
      hbeta0[j]<-(hbeta1[j]+hbeta1[j-1])*(zp[j]-zp[j-1])/2
    }
    Hbeta<-cumsum(hbeta0)  #####caculate regression function at zp
    beta_covar<-sapply(covar,function(x){
      ind<-order(abs(zp-x))[c(1,2)]
      betau<-mean(Hbeta[ind])
      return(betau)
    })#####calculate the regression function at covariate points
    
    Mut<-sapply(obstime,function(x){
      sum1<-sum(count*(obstime==x))
      sum2<-sum(exp(beta_covar)*(obstime==x))
      if(sum2==0){re<-0}else{re<-sum1/sum2}
    })######calculate the baseline function at observation times
    Tim<-seq(0.1,9.9,by=0.2)
    mu.t<-sapply(Tim,function(x){
      index<-which(between(obstime,x-0.1,x+0.1))
      mut<-mean(Mut[index])
      return(mut)
    })#####calculate the baseline function at Tim points
    for(k in 1:50){
      V<-sapply(obstime,function(x){
        sum0<-sum(ker(covar-zp[k],hopt[k])*exp(hbeta1[k]*(covar-zp[k])+hbeta2[k]*(covar-zp[k])^2)*(mtime>=x))
        sum1<-sum(ker(covar-zp[k],hopt[k])*exp(hbeta1[k]*(covar-zp[k])+hbeta2[k]*(covar-zp[k])^2)*((covar-zp[k])/hopt[k])*(mtime>=x))
        sum2<-sum(ker(covar-zp[k],hopt[k])*exp(hbeta1[k]*(covar-zp[k])+hbeta2[k]*(covar-zp[k])^2)*((covar-zp[k])/hopt[k])^2*(mtime>=x))
        sum3<-sum(ker(covar-zp[k],hopt[k])*exp(hbeta1[k]*(covar-zp[k])+hbeta2[k]*(covar-zp[k])^2)*((covar-zp[k])/hopt[k])^3*(mtime>=x))
        sum4<-sum(ker(covar-zp[k],hopt[k])*exp(hbeta1[k]*(covar-zp[k])+hbeta2[k]*(covar-zp[k])^2)*((covar-zp[k])/hopt[k])^4*(mtime>=x))
        if(sum0==0){re<-c(0,0,0,0,0)} else{re<-c(sum1/sum0,sum2/sum0,sum2/sum0-(sum1/sum0)^2,(sum0*sum3-sum1*sum2)/sum0^2,sum4/sum0-(sum2/sum0)^2)}
      })
      a11<-sum(hopt[k]*ker(covar-zp[k],hopt[k])^2*Mut^2*exp(2*beta_covar)*((covar-zp[k])/hopt[k]-V[1,])^2)/nsub 
      a12<-sum(hopt[k]*ker(covar-zp[k],hopt[k])^2*Mut^2*exp(2*beta_covar)*((covar-zp[k])/hopt[k]-V[1,])*(((covar-zp[k])/hopt[k])^2-V[2,]))/nsub 
      a22<-sum(hopt[k]*ker(covar-zp[k],hopt[k])^2*Mut^2*exp(2*beta_covar)*(((covar-zp[k])/hopt[k])^2-V[2,])^2)/nsub 
      b11<-sum(ker(covar-zp[k],hopt[k])*count*V[3,])/nsub
      b12<-sum(ker(covar-zp[k],hopt[k])*count*V[4,])/nsub
      b22<-sum(ker(covar-zp[k],hopt[k])*count*V[5,])/nsub
      A<-matrix(c(a11,a12,a12,a22),nr=2)
      B<-matrix(c(b11,b12,b12,b22),nr=2)    
      Sig1<-solve(B)%*%A%*%solve(B)
      Var[k]<-Sig1[1,1]/(nsub*hopt[k]^3)
    }#####calculate the estimated covariance
    All[[i]]<-list(Hbeta,hbeta1,Var,mu.t)
  }
  return(All)
} #the simulation study for 500 replications of the proposed method
simulation.ts<-function(nr,data){
  All<-list()
  tim<-seq(0.05,9.95,by=0.1)
  zp<-seq(0.04,2.98,by=0.06)
  for(i in 1:nr){
    set.seed(100*i+i)
    print(100*i+i)
    Fvalue<-shape.fun(data,tim)
    Fk<-Fvalue[[1]]
    Fall<-Fvalue[[2]]
    obstime<-data[,2]
    Ft<-sapply(obstime,function(x){
      mean(Fall[order(abs(tim-x))[1:3]])
    })
    
    bhopt1<-sapply(zp,function(x){
      re<-choose.hb.ts(data,z=x,Ft)
    })
    hopt1<-bhopt1[1,]
    gamma10<-bhopt1[2,]
    gamma11<-bhopt1[3,]
    
    hbeta10<-rep(0,50) 
    hbeta10[1]<-gamma11[1]*zp[1]  
    for(j in 2:50){
      hbeta10[j]<-(gamma11[j]+gamma11[j-1])*(zp[j]-zp[j-1])/2
    }
    Hbeta1<-cumsum(hbeta10) 
    Heta1<-mean(gamma10-Hbeta1)
    
    Muall1<-Fall*exp(Heta1)
    
    All[[i]]<-list(Hbeta1,gamma11,Muall1)
  }
  return(All)
}#the simulation study for 500 replications of the two-stage approach
nsub<-300
max.obs<-6
irate<-0.3
len<-10
beta<-c(0.01,0.01)
nr<-500
rdata<-data.gen(nsub,max.obs,irate,len)
mixdata<-data.gen.mix(nsub,max.obs,irate,len)
mixout<-simulation.est(nr,data = mixdata)
kerout<-simulation.est(nr,data=rdata)
tsout<-simulation.ts(nr,data=rdata)
save(kerout,file=paste(getwd(),"/kerpoisson.rda",sep=""))####save the simulation outcome
save(mixout,file=paste(getwd(),"/kermixpoisson.rda",sep=""))
save(tsout,file=paste(getwd(),"/tspoisson.rda",sep=""))
end<-Sys.time()
end-start
######calculate the regression function estimators, the derivative function estimators, 
######the estimated and empirical standard errors, the coverage probility, and the baseline function estimators
  tim<-seq(0.1,9.9,by=0.2)
  tim.ts<-seq(0.05,9.95,by=0.1)
  zk<-seq(0.04,2.98,by=0.06)
  trub<-2*zk/9
  Beta0.ts<-sapply(tsout,function(x)x[[1]])
  Beta1.ts<-sapply(tsout,function(x)x[[2]])
  Mu.ts<-sapply(tsout,function(x)x[[3]])
  beta0.ts<-apply(Beta0.ts, 1, mean)
  beta1.ts<-apply(Beta1.ts, 1, mean)
  esd.ts<-apply(Beta1.ts, 1, sd)
  mu.ts<-apply(Mu.ts, 1, mean)
  Beta0.ker<-sapply(kerout,function(x)x[[1]])  
  Beta1.ker<-sapply(kerout,function(x)x[[2]]) 
  Var.ker<-sapply(kerout,function(x)x[[3]]) 
  Mu.ker<-sapply(kerout,function(x)x[[4]])
  beta0.ker<-apply(Beta0.ker,1,mean)             
  beta1.ker<-apply(Beta1.ker,1,mean)
  esd.ker<-apply(Beta1.ker,1,sd)
  msd.ker<-apply(sqrt(Var.ker),1,mean)
  mu.ker<-apply(Mu.ker,1,mean)
  CIK<-lapply(1:50,function(k){          
    mapply(function(x,y){               
      lc<-mean(x)-1.96*sqrt(y)           
      rc<-mean(x)+1.96*sqrt(y)
      return(c(lc,rc))
    },Beta1.ker[k,],Var.ker[k,])
  })
  KCount<-sapply(1:50,function(x){      
    data.table::between(trub[x],CIK[[x]][1,],CIK[[x]][2,],incbounds = FALSE)   
  })                                                             
  KCovpro<-apply(KCount,2,function(x)sum(x)/500)
  
  Beta0.mix<-sapply(mixout,function(x)x[[1]])  
  Beta1.mix<-sapply(mixout,function(x)x[[2]]) 
  Var.mix<-sapply(mixout,function(x)x[[3]]) 
  Mu.mix<-sapply(mixout,function(x)x[[4]])
  beta0.mix<-apply(Beta0.mix,1,mean)             
  beta1.mix<-apply(Beta1.mix,1,mean)
  esd.mix<-apply(Beta1.mix,1,sd)
  msd.mix<-apply(sqrt(Var.mix),1,mean)
  mu.mix<-apply(Mu.mix,1,mean)
  CIM<-lapply(1:50,function(k){          
    mapply(function(x,y){               
      lc<-mean(x)-1.96*sqrt(y)           
      rc<-mean(x)+1.96*sqrt(y)
      return(c(lc,rc))
    },Beta1.mix[k,],Var.mix[k,])
  })
  MCount<-sapply(1:50,function(x){      
    data.table::between(trub[x],CIM[[x]][1,],CIM[[x]][2,],incbounds = FALSE)   
  })                                                             
  MCovpro<-apply(MCount,2,function(x)sum(x)/500)
  
  ind<-seq(3,48,by=3)
  zp<-zk[ind]
  MCsd.ts<-round(esd.ts[ind],3)
  MCsd.ker<-round(esd.ker[ind],3)
  MCsd.mix<-round(esd.mix[ind],3)
  ARE.ts<-MCsd.ts/MCsd.ker
  ARE.mix<-MCsd.ker/MCsd.mix
  
  ########plot
  postscript(file = paste(getwd(),"/cprandom.eps",sep = ""),onefile = FALSE,horizontal = FALSE)
  par(mfrow=c(3,2),mar= c(4, 4, 2, 2))
  plot(zk,beta1.ker,pch=".",cex=5,ann = F,xaxt="n",col=1,ylim = c(-0.2,1.2))
  axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
  title(xlab = "z", ylab = expression(paste(phi^(1),(z))),line = 2)
  lines(zk,beta1.ts,lty=2,lwd=2,col=1)
  lines(zk,2*zk/9,lwd=2,lty=1,col=1)
  legend("topleft",c("True curve","Proposed method","Two-satge method"),col = c(1,1,1),lty = c(1,3,2),lwd=c(1.5,1.8,1.5),bty = "o",cex = 0.9)
  title(sub = list("(c1)",cex=1.2,font=1),mgp=c(2,1,0))
  
  plot(zk,beta0.ker,pch=".",cex=5,ann = F,xaxt="n",col=1,ylim = c(0,1))
  axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
  title(xlab = "z", ylab = expression(paste(phi,(z))),line = 2)
  lines(zk,beta0.ts,lty=2,lwd=2,col=1)
  lines(zk,zk^2/9,lwd=2,lty=1,col=1)
  legend("topleft",c("True curve","Proposed method","Two-satge method"),col = c(1,1,1),lty = c(1,3,2),lwd=c(1.5,1.8,1.5),bty = "o",cex = 0.9)
  title(sub = list("(c2)",cex=1.2,font=1),mgp=c(2,1,0))
  
  plot(zk,esd.ker,pch=".",cex=5,col=1,ann = F,xaxt="n",ylim = c(0,1.5))
  axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
  title(xlab = "z",ylab = "S.E.",line = 2)
  lines(zk,esd.ts,lty=2.5,lwd=2,col=1)
  legend("top",c("Proposed method","Two-satge method"),col = c(1,1),lty = c(3,2),lwd=c(2,1.5),bty = "o",cex=0.9)
  title(sub = list("(c3)",cex=1.2,font=1),mgp=c(2,1,0))
  
  plot(zp,ARE.ts,typ="b",lty=2,lwd=2,col=1,ann = F,xaxt="n",ylim = c(0,2))
  axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
  abline(h=1,col=1,lwd=2)
  title(xlab = expression(z),ylab = expression("ARE"),line = 2)
  title(sub = list(expression("(c4)"),cex=1.2,font=1),mgp=c(2,1,0))
  
  plot(tim,mu.ker,pch=".",col=1,cex=5,ann = F,xaxt="n",ylim = c(2,40))
  axis(1,0:10,0:10)
  title(xlab = expression(t),ylab = expression(mu(t)),line = 2)
  lines(tim.ts,mu.ts,lty=2,lwd=2,col=1)
  lines(tim,8*sqrt(tim)+4,lwd=2,lty=1,col=1)
  legend("topleft",c("True curve","Proposed method","Two-satge method"),lwd = c(1.5,2,1.5),lty = c(1,3,2),col=c(1,1,1),bty = "o",cex = 0.9)
  title(sub = list(expression("(c5)"),cex=1.2,font=1),mgp=c(2,1,0))
  dev.off()
  
    postscript(file = paste(getwd(),"/cpmixpoisson.eps",sep = ""),onefile = FALSE,horizontal = FALSE)
    par(mfrow=c(3,2),mar= c(4, 4, 2, 2) )
    plot(zk,beta1.mix,pch=".",cex=5,ann = F,xaxt="n",col=1,ylim = c(-0.2,1))
    axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
    title(xlab = "z", ylab = expression(paste(phi^(1),(z))),line = 2)
    lines(zk,beta1.ker,lwd=2.5,lty=2,col=1)
    lines(zk,2*zk/9,lwd=2,lty=1,col=1)
    legend("topleft",c("True curve","Poisson process","MixedPoisson process"),col = c(1,1,1),lty = c(1,2,3),lwd=c(1.5,1.5,2),bty = "o",cex = 0.9)
    title(sub = list("(a1)",cex=1.2,font=1),mgp=c(2,1,0))
    
    plot(zk,beta0.mix,pch=".",cex=5,ann = F,xaxt="n",col=1,ylim = c(0,1))
    axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
    title(xlab = "z", ylab = expression(paste(phi,(z))),line = 2)
    lines(zk,beta0.ker,lty=2,lwd=2.5,col=1)
    lines(zk,zk^2/9,lwd=1.5,lty=1,col=1)
    legend("topleft",c("True curve","Poisson process","MixedPoisson process"),col = c(1,1,1),lty = c(1,2,3),lwd=c(1.5,1.5,2),bty = "o",cex = 0.9)
    title(sub = list("(a2)",cex=1.2,font=1),mgp=c(2,1,0))
    
    plot(zk,esd.ker,type = "l",lty=2,lwd=2,col=1,ann = F,xaxt="n",ylim = c(0,3))
    axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
    title(xlab = "z",ylab = "S.E.",line = 2)
    lines(zk,msd.ker,lty=3,lwd=3,col=1)
    lines(zk,esd.mix,lty=4,lwd=2,col=1)
    lines(zk,msd.mix,lty=5,lwd=2,col=1)
    legend("topleft",c("ESE_Poisson","MSE_Poisson","ESE_MixedPoisson","MSE_MixedPoisson"),col = c(1,1,1,1),lty = c(2,3,4,5),lwd=c(1.5,2,1.5,1.5),bty = "o",cex=0.9)
    title(sub = list("(a3)",cex=1.2,font=1),mgp=c(2,1,0))
    
    plot(zp,ARE.mix,typ="b",lty=2,lwd=2,col=1,ann = F,xaxt="n",ylim = c(0,2))
    axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
    abline(h=1,col=1,lwd=2)
    title(xlab = expression(z),ylab = expression("ARE"),line = 2)
    title(sub = list("(a4)",cex=1.2,font=1),mgp=c(2,1,0))
    
    plot(zk,MCovpro,pch=".",cex=4.5,ann = F,xaxt="n",col=1,ylim = c(0.5,1))
    axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
    title(xlab = "z", ylab = expression("cov. prob."),line = 2)
    lines(zk,KCovpro,lty=2,lwd=2,col=1)
    abline(h=0.95,col=1,lwd=2)
    legend("bottomleft",c("Poisson process","MixedPoisson process"),col = c(1,1),lty = c(2,3),lwd=c(1.5,2),bty = "o",cex = 0.9)
    title(sub = list("(a5)",cex=1.2,font=1),mgp=c(2,1,0))
    
    plot(tim,mu.mix,pch=".",col=1,cex=5,ann = F,xaxt="n",ylim = c(2,40))
    axis(1,0:10,0:10)
    title(xlab = expression(t),ylab = expression(mu(t)),line = 2)
    lines(tim,mu.ker,lty=2,lwd=2,col=1)
    lines(tim,8*sqrt(tim)+4,lwd=2,lty=1,col=1)
    legend("topleft",c("True curve","Poisson process","MixedPoisson process"),col = c(1,1,1),lty = c(1,2,3),lwd=c(1.5,1.5,2),bty = "o",cex = 0.9)
    title(sub = list("(a6)",cex=1.2,font=1),mgp=c(2,1,0))
    dev.off()
  
  
  
  
 
