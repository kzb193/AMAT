## Title: Adaptive Microbiome Association Test (AMAT)
## Version: 0.1
## Author: Kalins Banerjee (kbanerjee@pennstatehealth.psu.edu)
## Date: May 9, 2021

###################################################

## Inputs

## y: a numeric vector of continuous or binary outcome.
## Z: a matrix or data frame representing OTU table. Rows are samples, and columns are OTUs.
## X: a matrix or data frame of additional covariates that are to be adjusted for. Rows are samples.
## B: number of permutations. Default is B=500.
## model: "c" for continuous outcome, and "b" for binary outcome.

## total: a numeric vector providing total reads in the entire community for each sample.
#         If the test involves the entire community of OTUs, keep this as NULL.
#         If the test involves a subset of the entire community (e.g. an upper level taxon), it is recommended to specify this.
#         Default is NULL.


## Output: a list with the followings

## p-value of AMAT.
## The testing subset. This is a subset of the column names of Z, and the elements are arranged in descending order based on the  corresponding sample distance correlations.
## A vector with sample distance correlations corresponding to the testing subset above.

########################################
library(Rfast)


#######################################
AMAT<-function(y,Z,X,B=500,model,total=NULL){
  
  
  ## Internal function 1
  ## Computes sum of powered score (SPU) statistics, and returns the corresponding absolute values
  
  pow.score<-function(z0,rule,resi0,phi){
    
    sco<- (1/phi)*( t(scale(z0[,rule]))%*%resi0 ) 
    
    M<-abs( c( sum(sco^2),sum(sco^3),sum(sco^4),sum(sco^8)
    ))
    
    return(M)}
  
  
  ## Internal function 2
  ## computes p-value for adaptive sum of powered score (aSPU) statistic
  
  pv<-function(T.SPU.obs,T.SPU.perm,B){
    
    T.SPU.all<-rbind(T.SPU.obs,T.SPU.perm)
    p.SPU.all<-(B- apply(T.SPU.all,2,rank) +1 )/B
    
    aSPU.obs<-min(p.SPU.all[1,])
    
    p.SPU.perm<-((B-1)- apply(T.SPU.perm,2,rank) +1 )/(B-1)
    
    aSPU.perm<-apply(p.SPU.perm,1,min)
    
    p.aSPU<-sum(aSPU.perm<=aSPU.obs)/B
    
    return(p.aSPU)
  }
  
  
  ############## 
  n<-length(y)
  T.SPU.obs<-matrix(NA,1,4)
  z0<-Z
  
  ###############
  
  if(model=="c"){
    
    if(is.null(X)){
      f.lm.1<-glm(y~ 1,family = "gaussian")
      
      my.res<-(y-  f.lm.1$fitted.values)
    }else{
      
     
      f.lm.1<-glm(y~ .,family = "gaussian",data=as.data.frame(X))
      
      
      my.res<-(y-  f.lm.1$fitted.values)
      
      
    }
    
    phi<-sigma(f.lm.1)^2
    
  }
  #########################
  
  if(model=="b"){
    
    if(is.null(X)){
      f.lm.1<-glm(y~ 1,family = "binomial")
      
      my.res<-(y-  f.lm.1$fitted.values)
      my.res.w<-residuals.glm(f.lm.1,type="working")
      
    }else{
      
      
      f.lm.1<-glm(y~ X,family = "binomial",data=as.data.frame(X))
      
      my.res<-(y-  f.lm.1$fitted.values)
      my.res.w<-residuals.glm(f.lm.1,type="working")
      
    }
    
    phi<-1
  }
  ##################################
  if(is.null(total)==T){
    z0.comp<-z0/apply(z0,1,sum)
  }else{
    z0.comp<-z0/total
  }
  
  
  
  z0.scaled<-scale(z0)
  
  ########################
  if(model=="c"){
    test<-t(sapply(1:B,function(b){
      permu<-sample(1:n,n)
      r0<-my.res[permu]
      
      c(apply(z0.scaled,2,function(x) Rfast::dcor(x, r0)$dcor),
        r0)
    }
    )
    )
    
    DC0<-test[,1:ncol(z0)]
    resi0<-test[,-(1:(ncol(z0)) )]
    
    
    DC<-apply(z0.scaled,2,function(x) Rfast::dcor(x, my.res)$dcor)
    
    cutoff<-apply(DC0,2,mean)
    
    rule<-which(DC> cutoff )
    
    if(length(rule)==0){
      rule<-which.max(DC)
    }
    
    #################################
    
    T.SPU.obs[1,]<-pow.score(z0.comp,rule,my.res,phi)  
    
    #################################
    perm.all<-t( sapply(1:B, function(b){
      
      rule<-which(DC0[b,]>  cutoff)
      if(length(rule)==0){
        rule<-which.max(DC0[b,])
      } 
      
      
      c(pow.score(z0.comp,rule,resi0[b,],phi)
        
      )
    }
    )
    
    )
    T.SPU.perm<-  perm.all[,1:4]
    
    ###################################
    vv1<-DC
    vv2<-cutoff
    
    ddx<-cbind(which(vv1>vv2),vv1[which(vv1>vv2)])
    
    
    if(dim(ddx)[1]!=1){
      v_out<-ddx[order(ddx[,2],decreasing = T),][,2]
    }
    
    if(dim(ddx)[1]==1){
      v_out<-ddx[,2]
    }
    
    
    if(dim(ddx)[1]==0){
      
      v_out<-DC[which.max(DC)]
    }
    
  }
  ################################
  ################################
  if(model=="b"){
    test<-t(sapply(1:B,function(b){
      permu<-sample(1:n,n)
      r0<-my.res[permu]
      
      r0.w<-my.res.w[permu]
      
      c(
        apply(z0.scaled,2,function(x) Rfast::dcor(x, r0.w)$dcor),
        r0,r0.w)
    }
    )
    )
    
    ul<-ncol(z0)
    
    DC0.w<-test[,1:ncol(z0)]
    
    resi0<-test[,(ul+1):(ul+n)]
    
    resi0.w<-test[,-(1:(ul+n) )]
    
    
    DC.w<-apply(z0.scaled,2,function(x) Rfast::dcor(x, my.res.w)$dcor)
    
    
    cutoff.w<-apply(DC0.w,2,mean)
    
    
    rule.w<-which(DC.w> cutoff.w )
    
    if(length(rule.w)==0){
      rule.w<-which.max(DC.w)
    }
    
    ####################################
    T.SPU.obs[1,]<-pow.score(z0.comp,rule.w,my.res,phi)  
    
    ####################################
    
    perm.all<-t( sapply(1:B, function(b){
      
      rule.w<-which(DC0.w[b,]>  cutoff.w)
      if(length(rule.w)==0){
        rule.w<-which.max(DC0.w[b,])
      } 
      
      
      c(
        
        pow.score(z0.comp,rule.w,resi0[b,],phi)
      )
      
      
    }
    )
    
    )
    T.SPU.perm<-  perm.all[,1:4]
    
    ###############################    
    vv1<-DC.w
    vv2<-cutoff.w
    
    ddx<-cbind(which(vv1>vv2),vv1[which(vv1>vv2)])
    
    if(dim(ddx)[1]!=1){
      v_out<-ddx[order(ddx[,2],decreasing = T),][,2]
    }
    
    if(dim(ddx)[1]==1){
      v_out<-ddx[,2]
    }
    
    
    if(dim(ddx)[1]==0){
      
      v_out<-DC.w[which.max(DC.w)]
    }
    
  }
  
  
  #####################################################################################
  my.list<-list(pv(T.SPU.obs,T.SPU.perm,B),names(v_out),v_out
                
  )  
  
  
  return(my.list)}
####################################################################################
