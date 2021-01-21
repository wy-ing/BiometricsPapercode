rm(list = ls())
setwd("/Users/yangwang/Documents/paper_code")
library(data.table)
library(intervals)
f.mu<-function(x){8*sqrt(x)+4}
f.phi<-function(x){x^2/9}
start<-Sys.time()
##########################generate a data set##############
#create a dataset
data.gen<-function(
  nsub,                       #Sample size
  max.obs,                    #Maximal number of observation for each individual
  len                         #length of study
)
{
  Z1<-round(runif(nsub,0,3),4)         #define covariate
  NumObs<-rep(nsub,0)
  Mobs<-rep(nsub,0)
  Mcoun<-rep(nsub,0)
  
  indtimeObs<-seq(0.01,len-0.01,length=max.obs)  #interarrival times
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
  
  Mcoun[1]<-indcount[length(indtimeObs)]  #the event count at the last observation time
  lambda.t<-lambda
  timeObs<-indtimeObs          #define observation time
  count<-indcount              #define observation count
  covar1<-rep(Z1[1],NumObs[1]) #covariate for subject 1
  Mobstime<-rep(Mobs[1],NumObs[1]) #the maximum observation time   
  Mcount<-rep(Mcoun[1],NumObs[1])  #the last event count 
  
  #To expand the database by adding the subsequence observations
  for(i in 2:nsub){
    
    indtimeObs<-seq(0.01,len-0.01,length=max.obs)      #interarrival times
    NumObs[i]<-length(indtimeObs)     #the number of observations for the i subject
    Mobs[i]<-indtimeObs[length(indtimeObs)]
    
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
    covar1<-c(covar1,rep(Z1[i],NumObs[i]))#combine all the covariate
    Mobstime<-c(Mobstime,rep(Mobs[i],NumObs[i]))#combine all the maximum observation times
    Mcount<-c(Mcount,rep(Mcoun[i],NumObs[i]))#combine all the last event count
  }
  cbind(subj,timeObs,count,covar1,Mobstime,Mcount,lambda.t)
}   #end of generate data


ker<-function(x,h){ re<-0.75*(1-(x/h)^2)*(abs(x/h)<1)/h }   #epanechnikov kernel function

beta.fun.ker<-function(data,beta,z,h){
  obstime<-data[,2]
  count<-data[,3]
  covar<-data[,4]
  fun<-function(beta){
    a0<-sum(exp(beta[1]*(covar-z)+beta[2]*(covar-z)^2)*ker(covar-z,h))
    a1<-sum(exp(beta[1]*(covar-z)+beta[2]*(covar-z)^2)*ker(covar-z,h)*(covar-z))
    a2<-sum(exp(beta[1]*(covar-z)+beta[2]*(covar-z)^2)*ker(covar-z,h)*(covar-z)^2)
    a3<-sum(exp(beta[1]*(covar-z)+beta[2]*(covar-z)^2)*ker(covar-z,h)*(covar-z)^3)
    a4<-sum(exp(beta[1]*(covar-z)+beta[2]*(covar-z)^2)*ker(covar-z,h)*(covar-z)^4)
    if(a0==0){S<-c(0,0,0,0,0)}else{S<-c(a1/a0,a2/a0,a1^2/a0^2-a2/a0,a1*a2/a0^2-a3/a0,a2^2/a0^2-a4/a0)}
    sum11<-sum(ker(covar-z,h)*count*(covar-z-S[1]))
    sum12<-sum(ker(covar-z,h)*count*((covar-z)^2-S[2]))
    df<-c(sum11,sum12)
    sum21<-sum(ker(covar-z,h)*count*S[3])
    sum22<-sum(ker(covar-z,h)*count*S[4])
    sum23<-sum(ker(covar-z,h)*count*S[5])
    jac<-matrix(c(sum21,sum22,sum22,sum23),nr=2)
    list(df=df,jac=jac)
  }            
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
}     #newton-raphson iteration algorithm for the proposed method

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
      if(!is.na((beta-beta1)%*%(beta-beta1)<eps)&&(beta-beta1)%*%(beta-beta1)<eps){
        print(k)
        return(beta)
        break
      }
    }
  }          
  return(Newton(fun,beta))
} #newton-raphson iteration algorithm for the two-stage approach

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
    print(step)
  }  
  phat<-p
  Fhat<-sapply(tim,function(x){
    sum(phat[which(nordtime<=x)])
  })
  Flast<-c(Bmatrix%*%phat)
  return(list(Flast,Fhat))
}   #the estimates for shape function

cvscore.ker<-function(data,beta1,beta2,z,h){
  obstime<-data[,2]
  count<-data[,3]
  covar<-data[,4]
    a0<-sum(exp(beta1*(covar-z)+beta2*(covar-z)^2)*ker(covar-z,h))
    a1<-sum(exp(beta1*(covar-z)+beta2*(covar-z)^2)*ker(covar-z,h)*(covar-z))
    a2<-sum(exp(beta1*(covar-z)+beta2*(covar-z)^2)*ker(covar-z,h)*(covar-z)^2)
    a3<-sum(exp(beta1*(covar-z)+beta2*(covar-z)^2)*ker(covar-z,h)*(covar-z)^3)
    a4<-sum(exp(beta1*(covar-z)+beta2*(covar-z)^2)*ker(covar-z,h)*(covar-z)^4)
    if(a0==0){S<-c(0,0,0,0,0,0)}else{S<-c(a1/a0,a2/a0,a1^2/a0^2-a2/a0,a1*a2/a0^2-a3/a0,a2^2/a0^2-a4/a0,log(a0/nsub))}
  sum11<-sum(ker(covar-z,h)^2*count^2*(covar-z-S[1])^2)
  sum12<-sum(ker(covar-z,h)^2*count^2*((covar-z)^2-S[2])*(covar-z-S[1]))
  sum13<-sum(ker(covar-z,h)^2*count^2*((covar-z)^2-S[2])^2)
  dfm<-matrix(c(sum11,sum12,sum12,sum13),nr=2)
  sum21<-sum(ker(covar-z,h)*count*S[3])
  sum22<-sum(ker(covar-z,h)*count*S[4])
  sum23<-sum(ker(covar-z,h)*count*S[5])
  jac<-matrix(c(sum21,sum22,sum22,sum23),nr=2)
  cvl1<-sum(diag(solve(jac)%*%dfm))
  cvl2<-sum(ker(covar-z,h)*count*(beta1*(covar-z)+beta2*(covar-z)^2-S[6]))
  return(cvl1+cvl2)
}    #the cross-validation score for the proposed method
choose.hb.ker<-function(data,z){
  hvalue<-seq(0.5,0.9,by=0.2)
  cvh<-function(z,h){
    beta0<-beta.fun.ker(data,beta,z,h)
    cv<-cvscore.ker(data,beta0[1],beta0[2],z,h)
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
}#the cross-validation approach for bandwidth selection of the proposed method

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
}   #the cross-validation score function
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
}   #the cross-validation approach for bandwidth selection of the two-stage approach
###################simulation for 500 replications################
simulation.ker<-function(nr){
  All<-list()
  for(i in 1:nr){
    set.seed(100*i+i)
    print(100*i+i)
    data<-data.gen(nsub,max.obs,len)
    obstime<-data[,2]
    count<-data[,3]
    covar<-data[,4]
    mtime<-data[,5]
    zp<-seq(0.04,2.98,by=0.06)
    bhopt<-sapply(zp,function(x){
      re<-choose.hb.ker(data,z=x)
    })
    hopt<-bhopt[1,]
    hbeta1<-bhopt[2,]
    hbeta2<-bhopt[3,]
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
    })     #####calculate the regression function at covariate points
    Mut<-sapply(obstime,function(x){
      sum1<-sum(count*(obstime==x))
      sum2<-sum(exp(beta_covar)*(obstime==x))
      if(sum2==0){re<-0}else{re<-sum1/sum2}
    })   ######calculate the baseline function at observation times
    Tim<-seq(0.05,9.95,by=0.1)
    mu.t<-sapply(Tim,function(x){
      index<-order(abs(obstime-x))[1:2]
      mut<-mean(Mut[index])
      return(mut)
    })     #####calculate the baseline function at Tim points
    Var<-rep()
    for(k in 1:50){
      sum0<-sum(ker(covar-zp[k],hopt[k])*exp(hbeta1[k]*(covar-zp[k])+hbeta2[k]*(covar-zp[k])^2))
      sum1<-sum(ker(covar-zp[k],hopt[k])*exp(hbeta1[k]*(covar-zp[k])+hbeta2[k]*(covar-zp[k])^2)*((covar-zp[k])/hopt[k]))
      sum2<-sum(ker(covar-zp[k],hopt[k])*exp(hbeta1[k]*(covar-zp[k])+hbeta2[k]*(covar-zp[k])^2)*((covar-zp[k])/hopt[k])^2)
      sum3<-sum(ker(covar-zp[k],hopt[k])*exp(hbeta1[k]*(covar-zp[k])+hbeta2[k]*(covar-zp[k])^2)*((covar-zp[k])/hopt[k])^3)
      sum4<-sum(ker(covar-zp[k],hopt[k])*exp(hbeta1[k]*(covar-zp[k])+hbeta2[k]*(covar-zp[k])^2)*((covar-zp[k])/hopt[k])^4)
      if(sum0==0){V<-c(0,0,0,0,0)} else{V<-c(sum1/sum0,sum2/sum0,sum2/sum0-(sum1/sum0)^2,(sum0*sum3-sum1*sum2)/sum0^2,sum4/sum0-(sum2/sum0)^2)}
      
      a11<-sum(hopt[k]*ker(covar-zp[k],hopt[k])^2*Mut^2*exp(2*beta_covar)*((covar-zp[k])/hopt[k]-V[1])^2)/nsub 
      a12<-sum(hopt[k]*ker(covar-zp[k],hopt[k])^2*Mut^2*exp(2*beta_covar)*((covar-zp[k])/hopt[k]-V[1])*(((covar-zp[k])/hopt[k])^2-V[2]))/nsub 
      a22<-sum(hopt[k]*ker(covar-zp[k],hopt[k])^2*Mut^2*exp(2*beta_covar)*(((covar-zp[k])/hopt[k])^2-V[2])^2)/nsub 
      
      b11<-sum(ker(covar-zp[k],hopt[k])*count*V[3])/nsub
      b12<-sum(ker(covar-zp[k],hopt[k])*count*V[4])/nsub
      b22<-sum(ker(covar-zp[k],hopt[k])*count*V[5])/nsub
      
      A<-matrix(c(a11,a12,a12,a22),nr=2)
      B<-matrix(c(b11,b12,b12,b22),nr=2)    
      
      Sig1<-solve(B)%*%A%*%solve(B)
      Var[k]<-Sig1[1,1]/(nsub*hopt[k]^3)
    }   #####calculate the estimated covariance
    All[[i]]<-list(Hbeta,hbeta1,Var,mu.t)
  }
  return(All)
} ##the proposed method

simulation.ts<-function(nr){
  All<-list()
  tim<-seq(0.05,9.95,by=0.1)
  zp<-seq(0.04,2.98,by=0.06)
  for(i in 1:nr){
    set.seed(100*i+i)
    print(100*i+i)
    data<-data.gen(nsub,max.obs,len)
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
    gamma0<-bhopt1[2,]
    gamma1<-bhopt1[3,]
    
    hbeta10<-rep(0,50) 
    hbeta10[1]<-gamma1[1]*zp[1]  
    for(j in 2:50){
      hbeta10[j]<-(gamma1[j]+gamma1[j-1])*(zp[j]-zp[j-1])/2
    }
    Hbeta1<-cumsum(hbeta10) 
    Heta1<-mean(gamma0-Hbeta1)
    
    Muall1<-Fall*exp(Heta1)
    
    All[[i]]<-list(Hbeta1,gamma1,Muall1)
  }
  return(All)
} ##the two-stage approach
nsub<-300
max.obs<-10
irate<-0.3
len<-10
beta<-c(0.01,0.01)
nr<-500
kerout<-simulation.ker(nr)
tsout<-simulation.ts(nr)
save(kerout,file=paste(getwd(),"/kerfixed.rda",sep=""))  ####save the simulation outcome 
save(tsout,file=paste(getwd(),"/tsfixed.rda",sep=""))  ####save the simulation outcome 
end<-Sys.time()
end-start

#calculate the regression function estimators, the derivative function estimators, 
#the estimated and empirical standard errors, the coverage probility, and the baseline function estimators
  tim<-seq(0.05,9.95,by=0.1)
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
  CI<-lapply(1:50,function(k){          
    mapply(function(x,y){               
      lc<-mean(x)-1.96*sqrt(y)           
      rc<-mean(x)+1.96*sqrt(y)
      return(c(lc,rc))
    },Beta1.ker[k,],Var.ker[k,])
  })
  Count<-sapply(1:50,function(x){      
    data.table::between(trub[x],CI[[x]][1,],CI[[x]][2,],incbounds = FALSE)   
  })                                                             
  Covpro<-apply(Count,2,function(x)sum(x)/500)
  
  ind<-seq(3,48,by=3)
  zp<-zk[ind]
  MCsd.ts<-round(esd.ts[ind],3)
  MCsd.ker<-round(esd.ker[ind],3)
  ARE<-MCsd.ts/MCsd.ker

  ########plot
    postscript(file = paste(getwd(),"/cpfixed.eps",sep = ""),onefile = FALSE,horizontal = FALSE)
    par(mfrow=c(3,2),mar= c(4, 4, 2, 2))
    plot(zk,beta1.ker,pch=".",cex=5,ann = F,xaxt="n",col=1,ylim = c(-0.2,1.2))
    axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
    title(xlab = "z", ylab = expression(paste(phi^(1),(z))),line = 2)
    lines(zk,beta1.ts,lty=2,lwd=2,col=1)
    lines(zk,2*zk/9,lwd=2,lty=1,col=1)
    legend("topleft",c("True curve","Proposed method","Two-satge method"),col = c(1,1,1),lty = c(1,3,2),lwd=c(1.5,1.8,1.5),bty = "o",cex = 0.9)
    title(sub = list("(d1)",cex=1.2,font=1),mgp=c(2,1,0))
    
    plot(zk,beta0.ker,pch=".",cex=5,ann = F,xaxt="n",col=1,ylim = c(0,1))
    axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
    title(xlab = "z", ylab = expression(paste(phi,(z))),line = 2)
    lines(zk,beta0.ts,lty=2,lwd=2,col=1)
    lines(zk,zk^2/9,lwd=2,lty=1,col=1)
    legend("topleft",c("True curve","Proposed method","Two-satge method"),col = c(1,1,1),lty = c(1,3,2),lwd=c(1.5,1.8,1.5),bty = "o",cex = 0.9)
    title(sub = list("(d2)",cex=1.2,font=1),mgp=c(2,1,0))
    
    plot(zk,esd.ker,pch=".",cex=5,col=1,ann = F,xaxt="n",ylim = c(0,1.5))
    axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
    title(xlab = "z",ylab = "S.E.",line = 2)
    lines(zk,esd.ts,lty=2.5,lwd=2,col=1)
    legend("top",c("Proposed method","Two-satge method"),col = c(1,1),lty = c(3,2),lwd=c(2,1.5),bty = "o",cex=0.9)
    title(sub = list("(d3)",cex=1.2,font=1),mgp=c(2,1,0))
    
    plot(zp,ARE,typ="b",lty=2,lwd=2,col=1,ann = F,xaxt="n",ylim = c(0,2))
    axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
    abline(h=1,col=1,lwd=2)
    title(xlab = expression(z),ylab = expression("ARE"),line = 2)
    title(sub = list(expression("(d4)"),cex=1.2,font=1),mgp=c(2,1,0))
    
    plot(tim,mu.ker,type="l",col=1,lwd=1.5,ann = F,xaxt="n",ylim = c(2,40))
    axis(1,0:10,0:10)
    title(xlab = expression(t),ylab = expression(mu(t)),line = 2)
    lines(tim,mu.ts,lty=2,lwd=2,col=1)
    lines(tim,8*sqrt(tim)+4,lwd=2,lty=4,col=1)
    legend("topleft",c("True curve","Proposed method","Two-satge method"),lwd = c(1.5,1.5,1.5),lty = c(4,1,2),col=c(1,1,1),bty = "o",cex = 0.9)
    title(sub = list(expression("(d5)"),cex=1.2,font=1),mgp=c(2,1,0))
    dev.off()
    
    postscript(file = paste(getwd(),"/fixkxq3.eps",sep = ""),onefile = FALSE,horizontal = FALSE)
    par(mfrow=c(3,2),mar= c(4, 4, 2, 2))
    plot(zk,beta1.ker,pch=".",cex=5,ann = F,xaxt="n",col=1,ylim = c(-0.2,1))
    axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
    title(xlab = "z", ylab = expression(paste(phi^(1),(z))),line = 2)
    lines(zk,2*zk/9,lwd=2,lty=1,col=1)
    legend("topleft",c("True curve","Average curve"),col = c(1,1),lty = c(1,3),lwd=c(2,2),bty = "o",cex = 0.9)
    title(sub = list("(f1)",cex=1.2,font=1),mgp=c(2,1,0))
    
    plot(zk,beta0.ker,pch=".",cex=5,ann = F,xaxt="n",col=1,ylim = c(0,1))
    axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
    title(xlab = "z", ylab = expression(paste(phi,(z))),line = 2)
    lines(zk,zk^2/9,lwd=2,lty=1,col=1)
    legend("topleft",c("True curve","Average curve"),col = c(1,1),lty = c(1,3),lwd=c(2,2),bty = "o",cex = 0.9)
    title(sub = list("(f2)",cex=1.2,font=1),mgp=c(2,1,0))
    
    plot(zk,msd.ker,pch=".",cex=5,col=1,ann = F,xaxt="n",ylim = c(0,1))
    axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
    title(xlab = "z",ylab = "S.E.",line = 2)
    lines(zk,esd.ker,lty=2,lwd=2.5,col=1)
    legend("top",c("ESE","MSE"),col = c(1,1),lty = c(2,3),lwd=c(1.5,2),bty = "o",cex=0.9)
    title(sub = list("(f3)",cex=1.2,font=1),mgp=c(2,1,0))
    
    plot(zk,Covpro,pch=".",cex=5,ann = F,xaxt="n",col=1,ylim = c(0.8,1))
    axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
    title(xlab = "z", ylab = expression("cov. prob."),line = 2)
    abline(h=0.95,col=1,lwd=2)
    title(sub = list("(f4)",cex=1.2,font=1),mgp=c(2,1,0))
    
    plot(tim,mu.ker,pch=".",col=1,cex=3,ann = F,xaxt="n",ylim = c(2,40))
    axis(1,seq(0,10,by=1),seq(0,10,by=1))
    title(xlab = expression(t),ylab = expression(mu(t)),line = 2)
    lines(tim,8*sqrt(tim)+4,lwd=2,lty=1,col=1)
    legend("topleft",c("True curve","Average curve"),col = c(1,1),lty = c(1,3),lwd=c(2,2),bty = "o",cex = 0.9)
    title(sub = list("(f5)",cex=1.2,font=1),mgp=c(2,1,0))
    dev.off()
    
