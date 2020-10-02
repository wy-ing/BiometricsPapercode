rm(list = ls())
setwd("/Users/yangwang/Documents/paper_code")
library(data.table)
library(intervals)
f.mu<-function(x){8*sqrt(x)+4}
f.phi<-function(x){x^2/9}
start<-Sys.time()
#########create a dataset
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

ker<-function(x,h){ re<-0.75*(1-(x/h)^2)*(abs(x/h)<1)/h }   ########epanechnikov kernel function

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
} #########newton-raphson iteration algorithm

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
}   ###########the estimates for shape function

cvscore.ts<-function(data,beta1,beta2,z,h,Ft){
  obstime<-data[,2]
  count<-data[,3]
  covar<-data[,4]
  sum11<-sum(ker(covar-z,h)^2*(count/Ft-exp(beta[1]+beta[2]*(covar-z)))^2)
  sum12<-sum(ker(covar-z,h)^2*(covar-z)*(count/Ft-exp(beta[1]+beta[2]*(covar-z)))^2)
  sum13<-sum(ker(covar-z,h)^2*(covar-z)^2*(count/Ft-exp(beta[1]+beta[2]*(covar-z)))^2)
  dfm<-matrix(c(sum11,sum12,sum12,sum13),nr=2)
  sum21<--sum(ker(covar-z,h)*exp(beta[1]+beta[2]*(covar-z)))
  sum22<--sum(ker(covar-z,h)*(covar-z)*exp(beta[1]+beta[2]*(covar-z)))
  sum23<--sum(ker(covar-z,h)*(covar-z)^2*exp(beta[1]+beta[2]*(covar-z)))
  jac<-matrix(c(sum21,sum22,sum22,sum23),nr=2)
  
  cvl1<-sum(diag(solve(jac)%*%dfm))
  cvl2<-sum(ker(covar-z,h)*((beta[1]+beta[2]*(covar-z))*(count/Ft)-exp(beta[1]+beta[2]*(covar-z))))
  return(cvl1+cvl2)
}   #########the cross-validation score function
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
}   #########the cross-validation approach for bandwidth selection
###################################simulation for 500 replications###################
sumlation.ts<-function(nr){
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
    
    All[[i]]<-list(hopt1,gamma1,Hbeta1,Muall1)
  }
  return(All)
}
nsub<-300
max.obs<-10
irate<-0.3
len<-10
beta<-c(0.01,0.01)
nr<-500
tsout<-sumlation.ts(nr)
#save(out,file=paste(getwd(),"/tsfixed.rda",sep="")) ####save the simulation outcomes
end<-Sys.time()
end-start
####################calculate the regression function estimators, the derivative function estimators, 
###############the estimated and empirical standard errors, the coverage probilities, and the baseline function estimators
  hopt1<-sapply(out,function(x)x[[1]])
  Beta1<-sapply(out,function(x)x[[2]])
  Beta0<-sapply(out,function(x)x[[3]])
  Mut<-sapply(out,function(x)x[[4]])
  beta1<-apply(Beta1,1,mean)
  beta0<-apply(Beta0,1,mean)
  esd<-apply(Beta1, 1, sd)
  mu<-apply(Mut,1,mean)
  tim<-seq(0.05,9.95,by=0.1)
  zk<-seq(0.04,2.98,by=0.06) 
#####plot the regression function estimators, the derivative function estimators,
#####the estimated and empirical standard errors, the coverage probility, and the baseline function estimators
  par(mfrow=c(3,2),mar= c(4, 4, 2, 2))
  plot(zk,beta1,pch=".",cex=5,ann = F,xaxt="n",col=1,ylim = c(-0.2,1))
  axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
  title(xlab = "z", ylab = expression(paste(phi^(1),(z))),line = 2)
  lines(zk,2*zk/9,lwd=2,lty=1,col=1)
  legend("topleft",c("True curve","Average curve"),col = c(1,1),lty = c(1,3),lwd=c(2,2),bty = "o",cex = 0.9)
  title(sub = list("(a1)",cex=1.2,font=1),mgp=c(2,1,0))
  
  plot(zk,beta0,pch=".",cex=5,ann = F,xaxt="n",col=1,ylim = c(0,1))
  axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
  title(xlab = "z", ylab = expression(paste(phi,(z))),line = 2)
  lines(zk,zk^2/9,lwd=2,lty=1,col=1)
  legend("topleft",c("True curve","Average curve"),col = c(1,1),lty = c(1,3),lwd=c(2,2),bty = "o",cex = 0.9)
  title(sub = list("(a2)",cex=1.2,font=1),mgp=c(2,1,0))
  
  plot(zk,esd,type = "l",lty=2,lwd=3,col=1,ann = F,xaxt="n",ylim = c(0,5))
  axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
  title(xlab = "z",ylab = "S.E.",line = 2)
  title(sub = list("(a3)",cex=1.2,font=1),mgp=c(2,1,0))
  
  plot(tim,mu,pch=".",col=1,cex=3,ann = F,xaxt="n",ylim = c(2,40))
  axis(1,0:10,0:10)
  title(xlab = expression(t),ylab = expression(mu(t)),line = 2)
  lines(tim,8*sqrt(tim)+4,lwd=2,lty=1,col=1)
  legend("topleft",c("True curve","Average curve"),col = c(1,1),lty = c(1,3),lwd=c(2,2),bty = "o",cex = 0.9)
  title(sub = list("(a4)",cex=1.2,font=1),mgp=c(2,1,0))
  
  