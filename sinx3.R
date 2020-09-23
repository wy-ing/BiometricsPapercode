
rm(list = ls())
setwd("/Users/yangwang/Documents/paper_code")
library(data.table)
library(interval)
f.mu<-function(x){8*sqrt(x)+4}
#f.mu<-function(x){x^2+4}
f.phi<-function(x){x^2/9}
#f.phi<-function(x){sin(pi*x/3)}

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
  
  rexpf<-function(n,irate){
    u<-runif(n,0.001,0.999)     #1-u and u are same distribution
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
  
  lambda.t<-lambda
  timeObs<-indtimeObs          #define observation time
  count<-indcount              #define observation count
  covar1<-rep(Z1[1],NumObs[1]) #covariate1 for subject1
  Mobstime<-rep(Mobs[1],NumObs[1])     #generate censor time for subject1
  
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
    lambda.t<-c(lambda.t,lambda)
    timeObs<-c(timeObs,indtimeObs)      #combine all the observation
    count<-c(count,indcount)            #combine all the count
    covar1<-c(covar1,rep(Z1[i],NumObs[i]))#combine all the covariate1
    Mobstime<-c(Mobstime,rep(Mobs[i],NumObs[i]))
  }
  cbind(subj,timeObs,count,covar1,Mobstime,lambda.t)
  #cbind(subj,timeObs,count,covar1,censor,lambda.t)
}   #end of generate data

###########simulation beta############

ker<-function(x,h){  re<-0.75*(1-(x/h)^2)*(abs(x/h)<1)/h }             #epanechnikov kernel function
#ker<-function(x,h){exp(-x^2/(2*h^2))/(h*sqrt(2*pi))}   ###### Gaussian kernel function

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
  }           #newton-raphson iteration algorithm
  return(Newton(fun,beta))
}

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
}

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
}

simulation.est<-function(nr){
  All<-list()
  for(i in 1:nr){
    set.seed(100*i+i)
    print(100*i+i)
    
    data<-data.gen(nsub,max.obs,irate,len)
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
    
    Muk<-rep()
    Var<-rep()
    
    hbeta0<-rep(0,50) 
    hbeta0[1]<-hbeta1[1]*zp[1]  
    for(j in 2:50){
      hbeta0[j]<-(hbeta1[j]+hbeta1[j-1])*(zp[j]-zp[j-1])/2
    }
    Hbeta<-cumsum(hbeta0)  
    
    beta_covar<-sapply(covar,function(x){
      ind<-order(abs(zp-x))[c(1,2)]
      betau<-mean(Hbeta[ind])
      return(betau)
    })
    
    Mut<-sapply(obstime,function(x){
      sum1<-sum(count*(obstime==x))
      sum2<-sum(exp(beta_covar)*(obstime==x))
      if(sum2==0){re<-0}else{re<-sum1/sum2}
    })
    
    Tim<-seq(0.1,9.9,by=0.2)
    
    mu.t<-sapply(Tim,function(x){
      index<-which(between(obstime,x-0.1,x+0.1))
      mut<-mean(Mut[index])
      return(mut)
    })
    
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
      a22<-sum(hopt[k]*ker(covar-zp[k],hopt[k])^2*Mut^2*exp(2*beta_covar)*(((covar-zp[k])/hopt[k])^2-V[2,])^2)/nsub  #### V 有点偏大
      
      b11<-sum(ker(covar-zp[k],hopt[k])*count*V[3,])/nsub
      b12<-sum(ker(covar-zp[k],hopt[k])*count*V[4,])/nsub
      b22<-sum(ker(covar-zp[k],hopt[k])*count*V[5,])/nsub
      
      A<-matrix(c(a11,a12,a12,a22),nr=2)
      B<-matrix(c(b11,b12,b12,b22),nr=2)    
      
      Sig1<-solve(B)%*%A%*%solve(B)
      Var[k]<-Sig1[1,1]/(nsub*hopt[k]^3)
  
    }
    All[[i]]<-list(hopt,Hbeta,hbeta1,Var,mu.t)
  }
  return(All)
}

nsub<-300
max.obs<-6
irate<-0.3
len<-10
beta<-c(0.001,0.001)
nr<-500

out<-simulation.est(nr)
save(out,file=paste(getwd(),"/sinx3.rda",sep=""))
end<-Sys.time()
end-start


if(FALSE){
  library(data.table)
  #tim<-seq(0.05,9.95,by=0.1)
  tim1<-seq(0.1,9.9,by=0.2)
  #zk<-seq(0.01,2.98,by=0.03)
  zk<-seq(0.04,2.98,by=0.06)
  #trub<-pi*cos(pi*zk/3)/3
  trub<-2*zk/9
  Hop<-sapply(out,function(x)x[[1]])  
  B0<-sapply(out,function(x)x[[2]])  
  B1<-sapply(out,function(x)x[[3]])  
  V1<-sapply(out,function(x)x[[4]])  
  Mu1<-sapply(out,function(x)x[[5]])
  hop<-apply(Hop,1,mean)
  ave0<-apply(B0,1,mean)             
  ave1<-apply(B1,1,mean)
  esd<-apply(B1,1,sd)
  msd<-apply(sqrt(V1),1,mean)
  amu1<-apply(Mu1,1,mean)
  
  CI<-lapply(1:50,function(k){          
    mapply(function(x,y){               
      lc<-mean(x)-1.96*sqrt(y)           
      rc<-mean(x)+1.96*sqrt(y)
      return(c(lc,rc))
    },B1[k,],V1[k,])
  })
  Count<-sapply(1:50,function(x){      
    data.table::between(trub[x],CI[[x]][1,],CI[[x]][2,],incbounds = FALSE)   
  })                                                             
  Covpro<-apply(Count,2,function(x)sum(x)/500)
  
  ind<-seq(3,48,by=3)
  zp<-zk[ind]
  MCbias<-round((ave1-2*zk/9)*100,3)[ind]
  MCbias1<-round((deta1-2*zk/9)*100,3)[ind]
  #MCbias<-round((ave1-pi*cos(pi*zk/3)/3)*100,3)[ind]
  #MCbias1<-round((deta1-pi*cos(pi*zk/3)/3)*100,3)[ind]
  MCsd<-round(esd[ind],3)
  MCsd1<-round(esd1[ind],3)
  #Thsd<-round(msd[ind],3)
  #Thsd1<-round(msd1[ind],3)
  ARE<-MCsd1/MCsd
  
  if(FALSE){
    par(mfrow=c(3,2),mar= c(4, 4, 2, 2))
    plot(zk,ave1,pch=".",cex=4.5,ann = F,xaxt="n",col=1,ylim = c(-0.2,1))
    axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
    title(xlab = "z", ylab = expression(paste(phi^(1),(z))),line = 2)
    lines(zk,deta1,lty=2,lwd=2,col=1)
    #lines(zk,pi*cos(pi*zk/3)/3,lwd=3,lty=1,col=1)
    lines(zk,2*zk/9,lwd=2,lty=1,col=1)
    legend("topleft",c("True curve","Proposed method","Two-satge method"),col = c(1,1,1),lty = c(1,3,2),lwd=c(1.5,1.8,1.5),bty = "o",cex = 0.9)
    #legend("topleft",c("True curve","Average curve"),col = c(1,1),lty = c(1,3),lwd=c(2,2),bty = "n",cex = 0.9)
    title(sub = list("(e1')",cex=1.2,font=1),mgp=c(2,1,0))
    
    plot(zk,ave0,pch=".",cex=4.5,ann = F,xaxt="n",col=1,ylim = c(0,1))
    axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
    title(xlab = "z", ylab = expression(paste(phi,(z))),line = 2)
    lines(zk,beta1,lty=2,lwd=2,col=1)
    #lines(zk,sin(pi*zk/3),lwd=3,lty=1,col=1)
    lines(zk,zk^2/9,lwd=2,lty=1,col=1)
    legend("topleft",c("True curve","Proposed method","Two-satge method"),col = c(1,1,1),lty = c(1,3,2),lwd=c(1.5,1.8,1.5),bty = "o",cex = 0.9)
    #legend("topleft",c("True curve","Average curve"),col = c(1,1),lty = c(1,3),lwd=c(2,2),bty = "n",cex = 0.9)
    title(sub = list("(e2')",cex=1.2,font=1),mgp=c(2,1,0))
    
    plot(zk,esd,type = "l",lty=3,lwd=3,col=1,ann = F,xaxt="n",ylim = c(0,3))
    axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
    title(xlab = "z",ylab = "S.E.",line = 2)
    #points(zk,esd,pch=".",cex=4.5,col=1)
    lines(zk,esd1,lty=2.5,lwd=2,col=1)
    legend("top",c("Proposed method","Two-satge method"),col = c(1,1),lty = c(3,2),lwd=c(2,1.5),bty = "o",cex=0.9)
    title(sub = list("(e3')",cex=1.2,font=1),mgp=c(2,1,0))
    
    plot(zp,ARE,typ="b",lty=2,lwd=2,col=1,ann = F,xaxt="n",ylim = c(0,2))
    axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
    abline(h=1,col=1,lwd=2)
    title(xlab = expression(z),ylab = expression("ARE"),line = 2)
    title(sub = list("(e4')",cex=1.2,font=1),mgp=c(2,1,0))
    
    plot(tim1,amu2,pch=".",cex=5,ann = F,xaxt="n",col=1,ylim = c(0,40))
    #plot(tim,amu2,type="l",col=1,lwd=2,ann = F,xaxt="n",ylim = c(2,40))
    axis(1,0:10,0:10)
    title(xlab = expression(t),ylab = expression(mu(t)),line = 2)
    lines(tim,mu1,lty=2,lwd=2,col=1)
    lines(tim,8*sqrt(tim)+4,lwd=2,lty=1,col=1)
    #lines(tim,tim^2+4,lwd=3,lty=1,col=1)
    legend("topleft",c("True curve","Proposed method","Two-satge method"),lwd = c(1.5,1.8,1.5),lty = c(1,3,2),col=c(1,1,1),bty = "o",cex = 0.9)
    #legend("topleft",c("True curve","Average curve"),col = c(1,1),lty = c(1,3),lwd=c(2,2),bty = "n",cex = 0.9)
    title(sub = list("(e5')",cex=1.2,font=1),mgp=c(2,1,0))

    
    #plot(zk,Covpro,lty=2,cex=0.3,ann = F,xaxt="n",col=3,ylim = c(0.5,1))
    plot(zk,Covpro,pch=".",cex=4.5,ann = F,xaxt="n",col=1,ylim = c(0.5,1))
    axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
    title(xlab = "z", ylab = expression("cov. prob."),line = 2)
    abline(h=0.95,col=1,lwd=2.5)
    title(sub = list("(b4)",cex=1.2,font=1),mgp=c(2,1,0))
  }
  par(mfrow=c(1,1))
  plot(zk,hop,pch=".",lty=3,cex=3,ann = F,xaxt="n",col=4,ylim = c(0,1))
  points(zk,hop1,pch=".",cex=3,lty=1,col=3)
  legend("bottomleft",c("n300","n500"),col = c(3,4),lty = c(3,3),lwd=c(3,3),bty = "n",cex=0.9)
  
  plot(zk,ave1,pch=".",lty=3,cex=4,ann = F,xaxt="n",col="blue",ylim = c(-1.2,1.2))
  #axis(1,0:3,0:3)
  axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
  title(xlab = "z", ylab = expression(paste(phi^(1),(z))),line = 2)
  points(zk,ave11,pch=".",col=5,cex=4)
  lines(zk,pi*cos(pi*zk/3)/3,lwd=3,lty=1,col=3)
  #lines(zk,2*zk/9,lwd=3,lty=1,col=3)
  #legend("topleft",c("True curve","Average curve"),col = c(3,4),lty = c(1,3),lwd=c(2.5,3),bty = "n",cex=0.9)
  legend("bottomleft",c("True curve","poisson","mixed-poisson"),col = c(3,4,5),lty = c(1,3,3),lwd=c(2.5,2.5,2.5),cex = 0.9,bty = "n")
  #legend("bottomright",c("true","h=0.5","h=0.8"),col = c(3,4,2),lty = c(1,3,3),lwd=c(2,2,2),cex = 0.8)
  title(main = list("n=300,simu=500,point=50",cex=1.2,col="blue",font=3))
  
  plot(zk,msd,pch = ".",lty=3,cex=4,col=4,ann = F,xaxt="n",ylim = c(0,2))
  axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
  title(xlab = "z",ylab = "S.E.",line = 2)
  points(zk,msd1,pch=".",cex=4,col=3)
  #points(zk,esd,pch=".",cex=4,col=2)
  #points(zk,esd1,pch=".",cex=4,col=5)
  #legend("top",c("MSE 300","ESE 300","MSE 500","ESE 500"),col = c(3,4,2,5),lty = c(3,3,3,3),lwd=c(2.5,2.5,2.5,2.5),bty = "n",cex = 0.9)
  legend("top",c("MSE-poisson","MSE-mixedpoisson"),col = c(4,3),lty = c(3,3),lwd=c(3,3),bty = "n",cex=0.9)
  #legend("top",c("ESE","MSE"),col = c(3,4),lty = c(3,3),lwd=c(3,3),bty = "n",cex=0.9)
  title(main = list("n=300,simu=500,point=50",cex=1.2,col="blue",font=1.5))
  
  plot(zk,ave0,pch=".",cex=4,ann = F,xaxt="n",col=4,ylim = c(0,1))
  axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
  title(xlab = "z", ylab = expression(paste(phi,(z))),line = 2)
  points(zk,ave01,pch=".",cex=4,col=5)
  lines(zk,sin(pi*zk/3),lwd=3,lty=1,col=3)
  #lines(zk,zk^2/9,lwd=3,lty=1,col=3)
  #legend("topleft",c("True curve","Average curve"),col = c(3,4),lty = c(1,3),lwd=c(2.5,3),bty = "n",cex = 0.9)
  legend("bottom",c("true curve","poisson","mixed_poisson"),col = c(3,4,5),lty = c(1,3,3),lwd=c(2.5,2.5,2.5),cex = 0.9,bty = "n")
  title(main = list("n=300,point=50,simu=500",cex=1.2,col="blue",font=3))
  
  #plot(zk,Covpro,lty=2,cex=0.3,ann = F,xaxt="n",col=3,ylim = c(0.5,1))
  plot(zk,Covpro,pch=".",cex=4,ann = F,xaxt="n",col=4,ylim = c(0.5,1))
  axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
  title(xlab = "z", ylab = expression("cov. prob."),line = 2)
  points(zk,Covpro1,col=3,pch=".",cex=4)
  abline(h=0.95,col=1,lwd=2)
  legend("bottom",c("poisson","mixedpoisson"),col = c(4,3),lty = c(3,3),lwd=c(3,3),bty = "n",cex=0.9)
  title(main = list("n=300,simu=500,point=50",cex=1.2,col="blue",font=3))
  
  plot(tim1,amu1,pch=".",col=4,cex=4,ann = F,xaxt="n",ylim = c(4,140))
  #axis(1,seq(0,10,length=10),seq(0,10,length=10))
  axis(1,0:10,0:10)
  title(xlab = expression(t),ylab = expression(mu(t)),line = 2)
  points(tim,amu11,pch=".",cex=4,col=5)
  #points(tim,amu2,pch=".",cex=4,col=2)
  #points(tim,amu21,pch=".",cex=4,col=6)
  lines(tim,tim^2+4,lwd=2.5,lty=1,col=3)
  lines(tim,8*sqrt(tim)+4,lwd=3,lty=1,col=3)
  #legend("topleft",c("true curve","poisson","mixed_poisson"),col = c(3,4,5),lty = c(1,3,3),lwd=c(2.5,2.5,2.5),cex = 0.9,bty = "n")
  #legend("topleft",c("True curve","Average curve"),col = c(3,4),lty = c(1,3),lwd=c(2.5,3),bty = "n",cex = 0.9)
  title(main = list("n=300,simu=500,point=50",cex=1.2,col="blue",font=3))

#zp=seq(0.28,2.92,by=0.24)
#index<-sapply(zp,function(x){which(zk==x)})
  ind<-seq(3,48,by=3)
  zp<-zk[ind]
  #MCbias<-round((ave1-2*zk/9)*100,3)[ind]
  #MCbias1<-round((ave11-2*zk/9)*100,3)[ind]
  MCbias<-round((ave1-pi*cos(pi*zk/3)/3)*100,3)[ind]
  MCbias1<-round((ave11-pi*cos(pi*zk/3)/3)*100,3)[ind]
  MCsd<-round(esd[ind],3)
  MCsd1<-round(esd1[ind],3)
  Thsd<-round(msd[ind],3)
  Thsd1<-round(msd1[ind],3)
  Mcp<-data.frame(zp,MCbias,MCsd,Thsd)
  names(Mcp)<-c("Point","Monte-Carlo Bias","Monte-Carlo s.d.","Theoretical s.d")
  Mcp1<-data.frame(zp,MCbias1,MCsd1,Thsd1)
  names(Mcp1)<-c("Point","Monte-Carlo Bias","Monte-Carlo s.d.","Theoretical s.d")
  #xtable(Mcp,digits=3,caption="Poisson process")
  xtable(Mcp,caption="Poisson process")  ###生成latex语言格式的表格
  xtable(Mcp1,caption="mixed-Poisson process")
  ARE<-Thsd/Thsd1
  
  plot(zp,ARE,typ="b",lty=2,lwd=2,col=4,ann = F,xaxt="n",ylim = c(0.8,1.2))
  axis(1,seq(0,3,by=0.5),seq(0,3,by=0.5))
  title(xlab = expression(z),ylab = expression(ARE),line = 2)
  title(main = list("n=300,simu=500,point=50",cex=1.2,col="blue",font=3))
}

