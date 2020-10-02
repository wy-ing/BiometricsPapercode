rm(list = ls())
setwd("/Users/yangwang/Documents/Biometrics-reviews/EMEQ code")
#setwd("/home/sam/Documents/WY/code")
library(data.table)
library(intervals)
#f.mu<-function(x){8*sqrt(x)+4}
f.mu<-function(x){x^2+4}
#f.phi<-function(x){x^2/9}
f.phi<-function(x){sin(pi*x/3)}

start<-Sys.time()

#create a dataset
data.gen<-function(
  nsub,                       #Sample size
  max.obs,                    #Maximal number of observation for each individual
  irate,                      #Rate for interarrival time
  len                         #length of study
)
{
  Z1<-round(runif(nsub,0,3),4)         #define covariate 1
  #C<-runif(nsub,0,len+2)      #generate censor times 
  #To creat a database starting from the first subject
  NumObs<-rep(nsub,0)
  Mobs<-rep(nsub,0)
  Mcoun<-rep(nsub,0)
  
  rexpf<-function(n,irate){
    u<-runif(n,0.001,0.999)     #1-u and u are same distribution; avoid generate the zero
    x<--log(u)/irate          #interarrival time x~exp(-irate*x) or x~1-exp(-irate*x),thus avoid zero data
    return(x)
  }
  
  x<-round(rexpf(max.obs,irate),4)   #generate time for the first subject
  while(x[1]>=len)            #if the first time larger than the length of study,generate time again
    x<-round(rexpf(max.obs,irate),4)
  for(j in 1:max.obs){
    if(sum(x[1:j])<len)
      k<-j
  }                           #choose the cumulative time less than the length of study
  indtimeObs<-cumsum(x[1:k])  #interarrival time
  NumObs[1]<-length(indtimeObs) #the number of observations for the first subject
  Mobs[1]<-indtimeObs[length(indtimeObs)]
  
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
  covar1<-rep(Z1[1],NumObs[1]) #covariate1 for subject1
  Mobstime<-rep(Mobs[1],NumObs[1])     #generate censor time for subject1
  Mcount<-rep(Mcoun[1],NumObs[1])
  
  #To expand the database by adding the subsequence observations
  for(i in 2:nsub){
    x<-round(rexpf(max.obs,irate),4)  #generate time for the i subject
    while(x[1]>=len)
      x<-round(rexpf(max.obs,irate),4)#if the first time larger than the length of study,generate time again
    for(j in 1:max.obs){
      if(sum(x[1:j])<len)
        k<-j
    }                                 #choose the cumulative time less than the length of study
    indtimeObs<-cumsum(x[1:k])        #interarrival time
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
    covar1<-c(covar1,rep(Z1[i],NumObs[i]))#combine all the covariate1
    Mobstime<-c(Mobstime,rep(Mobs[i],NumObs[i]))
    Mcount<-c(Mcount,rep(Mcoun[i],NumObs[i]))
  }
  cbind(subj,timeObs,count,covar1,Mobstime,Mcount,lambda.t)
  #cbind(subj,timeObs,count,covar1,censor,lambda.t)
}   #end of generate data


ker<-function(x,h){ re<-0.75*(1-(x/h)^2)*(abs(x/h)<1)/h }             #epanechnikov kernel function
beta.fun<-function(data,beta,z,h,Ft){
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
  }           #newton-raphson iteration algorithm
  return(Newton(fun,beta))
}

shape.fun<-function(data,tim){
  Insub<-data[,1]
  obstime<-data[,2]
  count<-data[,3]
  covar<-data[,4]
  mtime<-data[,5]
  mcount<-data[,6]
  unmtime<-unname(tapply(mtime,Insub,unique))
  unmcount<-unname(tapply(mcount,Insub,unique))
  ordtime<-unique(sort(obstime))  #### outcome more unbias
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
  }  #####itrate 6900 (1e-5); 3154 (1e-4)
  phat<-p
  Fhat<-sapply(tim,function(x){
    sum(phat[which(nordtime<=x)])
  })
  Flast<-c(Bmatrix%*%phat)
  return(list(Flast,Fhat))
}
cvscore1<-function(data,beta1,beta2,z,h,Ft){
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
}
choose.hb1<-function(data,z,Ft){
  hvalue<-seq(0.5,0.9,by=0.2)
  cvh<-function(z,h){
    beta0<-beta.fun(data,beta,z,h,Ft)
    cv<-cvscore1(data,beta0[1],beta0[2],z,h,Ft)
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
}

est.sumlation<-function(nr,h){
  All<-list()
  #tim<-seq(0.02,10,length=200)
  tim<-seq(0.05,9.95,by=0.1)
  zp<-seq(0.04,2.98,by=0.06)
  for(i in 1:nr){
    set.seed(100*i+i)
    print(100*i+i)
    data<-data.gen(nsub,max.obs,irate,len)
    Fvalue<-shape.fun(data,tim)
    Fk<-Fvalue[[1]]
    Fall<-Fvalue[[2]]
    obstime<-data[,2]
    Ft<-sapply(obstime,function(x){
      mean(Fall[order(abs(tim-x))[1:3]])
    })
    
    bhopt1<-sapply(zp,function(x){
      re<-choose.hb1(data,z=x,Ft)
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
    
    All[[i]]<-list(hopt1,gamma11,Hbeta1,Muall1)
  }
  return(All)
}

nsub<-300
max.obs<-6
irate<-0.3
len<-10
beta<-c(0.01,0.01)
nr<-500

out<-est.sumlation(nr,h)
save(out,file=paste(getwd(),"/twostagecvss3.rda",sep=""))
end<-Sys.time()
end-start

if(FALSE){
  hopt1<-sapply(out,function(x)x[[1]])
  Deta1<-sapply(out,function(x)x[[2]])
  Beta1<-sapply(out,function(x)x[[3]])
  Mut<-sapply(out,function(x)x[[4]])
  opt1<-apply(hopt1,1,mean)
  deta1<-apply(Deta1,1,mean)
  beta1<-apply(Beta1,1,mean)
  esd1<-apply(Deta1, 1, sd)
  mu1<-apply(Mut,1,mean)
tim<-seq(0.05,9.95,by=0.1)
#tim<-seq(0.02,10,length=100)
zp<-seq(0.04,2.98,by=0.06) 

plot(zp,opt1,pch=".",lty=3,cex=4,ann = F,xaxt="n",col=4,ylim = c(0,1))
axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
title(xlab = "z", ylab = expression(h),line = 2)
points(zp,hop,pch=".",lty=3,cex=4,col=5)
title(main = list("n=500,simu=500,random",cex=1.2,col="blue",font=3))

plot(zp,deta1,pch=".",lty=3,cex=4,ann = F,xaxt="n",col=4,ylim = c(0,1))
axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
title(xlab = "z", ylab = expression(paste(phi^(1),(z))),line = 2)
points(zp,ave1,pch=".",lty=3,cex=4,col=5)
#lines(zp,pi*cos(pi*zp/3)/3,lwd=2.5,lty=1,col=3)
lines(zp,zp*2/9,lwd=2.5,lty=1,col=3)
#legend("topleft",c("True curve","Average curve"),col = c(3,4),lty = c(1,3),lwd=c(2.5,3),bty = "n",cex=0.9)
legend("topleft",c("True curve","kernel","twostage"),col = c(3,5,4),lty = c(1,3,3),lwd=c(2.5,2.5,2.5),cex = 0.9,bty = "n")
#legend("bottomright",c("true","h=0.5","h=0.8"),col = c(3,4,2),lty = c(1,3,3),lwd=c(2,2,2),cex = 0.8)
title(main = list("n=500,simu=500,random",cex=1.2,col="blue",font=3))

plot(zp,beta1,pch=".",lty=3,cex=4,ann = F,xaxt="n",col=4,ylim = c(0,1))
axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
title(xlab = "z", ylab = expression(paste(phi,(z))),line = 2)
points(zp,ave0,pch=".",lty=3,cex=4,col=5)
#lines(zp,sin(pi*zp/3),lwd=2.5,lty=1,col=3)
lines(zp,zp^2/9,lwd=2.5,lty=1,col=3)
#legend("topleft",c("True curve","Average curve"),col = c(3,4),lty = c(1,3),lwd=c(2.5,3),bty = "n",cex=0.9)
legend("topleft",c("True curve","kernel","twostage"),col = c(3,5,4),lty = c(1,3,3),lwd=c(2.5,2.5,2.5),cex = 0.9,bty = "n")
#legend("bottomright",c("true","h=0.5","h=0.8"),col = c(3,4,2),lty = c(1,3,3),lwd=c(2,2,2),cex = 0.8)
title(main = list("n=500,simu=500,random",cex=1.2,col="blue",font=3))

plot(zp,esd1,pch=".",lty=3,cex=4,ann = F,xaxt="n",col=3,ylim = c(0,2))
axis(1,seq(0,10,by=1),seq(0,10,by=1))
title(xlab = "z", ylab = "S.E.",line = 2)
points(zp,esd,pch=".",lty=3,cex=4,col=4)
legend("top",c("kernel","twostage"),col = c(4,3),lty = c(3,3),lwd=c(2.5,2.5),bty = "n",cex=0.9)
title(main = list("n=500,simu=500,random",cex=1.2,col="blue",font=3))

plot(tim,mu1,pch=".",lty=3,cex=4,ann = F,xaxt="n",col=4,ylim = c(0,40))
axis(1,seq(0,10,by=1),seq(0,10,by=1))
title(xlab = "time", ylab = expression(paste(mu,(t))),line = 2)
points(tim1,amu1,pch=".",lty=3,cex=4,col=5)
#lines(tim,tim^2+4,lwd=2.5,lty=1,col=3)
lines(tim,8*sqrt(tim)+4,lwd=2.5,lty=1,col=3)
legend("topleft",c("True curve","kernel","twostage"),col = c(3,5,4),lty = c(1,3,3),lwd=c(2.5,2.5,2.5),cex = 0.9,bty = "n")
#legend("topleft",c("True curve","Average curve"),col = c(3,4),lty = c(1,3),lwd=c(2.5,3),bty = "n",cex=0.9)
title(main = list("n=500,simu=500,random",cex=1.2,col="blue",font=3))
}

