rm(list=ls()) 
#-------------------------------------------------------------------------------
#Packages
#-------------------------------------------------------------------------------
library(lhs)
library(EzGP)
library(laGP)
#serial job
f=function(l){
#-------------------------------------------------------------------------------
#initial data
#-------------------------------------------------------------------------------
  method_name="ARSD-LCB"
  print(method_name)
  load("initial.RData")
#-------------------------------------------------------------------------------
#n0 data
#-------------------------------------------------------------------------------
  filename= paste0("data/tradata","_",l,".csv",sep="")
  n0_dat=as.matrix(read.csv(file = filename))[,(2:(p+q+2))]
  model = EzGP_fit(n0_dat[,1:(p+q)],n0_dat[,(p+q+1)],p=p,q=q,m=m,tau=tau)
#-------------------------------------------------------------------------------
#measurements and time
#-------------------------------------------------------------------------------
  #M_c0
  M_c0_mean=c()
  #min|y-a|
  y_a=c()
  #time
  time=c()
#-------------------------------------------------------------------------------
#Criterion functions
#-------------------------------------------------------------------------------
  #ARSD
  ARSD_f <- function(ypred,mse,a){
    
    alpha=0.05
    mu=abs(ypred-a)
    beta_h=2*log(((pi^2)*((nrow(n0_dat))^2)*(prod(m)))/(6*(alpha)))
    LCB=(mu)-(sqrt(beta_h))*(sqrt(mse))
    UCB=(mu)+(sqrt(beta_h))*(sqrt(mse))
    est=list(LCB=LCB,UCB=UCB)
    return(est)
  }
  
  #LCB
  LCB_f=function(x,model,a){
    if(is.matrix(x)!=TRUE) x = matrix(x,nrow=1,ncol=p+q)
    pred= EzGP_predict(x[,1:(p+q)], model, MSE_on = 1)
    ypred = pred$Y_hat
    mse = pred$MSE
    for(ii in 1:length(mse)){ if(mse[ii]<0){
      mse[ii]=0
    }
    }
    rho=2
    mu=abs(ypred-a)
    LCB=(mu)-(rho)*(sqrt(mse))
    return(LCB)
  }
#-------------------------------------------------------------------------------
#Adaptive design
#-------------------------------------------------------------------------------
  for(i in 1:(nnew)){
    print(i)
#------------------------------------------------------
#Search space
#------------------------------------------------------
    if(q==1){
      xpred=matrix(NA,nrow=npred*prod(m),ncol=p)
      for(ij in 1:prod(m)){
        xpred[(((ij-1)*(npred))+1):(ij*npred),] = randomLHS(npred,p)
        
      }
      zpred=rep(1:m,each=npred)
    }else{   
      xpred=matrix(NA,nrow=npred*prod(m),ncol=p)
      zpred=matrix(NA,nrow=npred*prod(m),ncol=q)
      for(ij in 1:prod(m)){
        xpred[(((ij-1)*(npred))+1):(ij*npred),] = randomLHS(npred,p)
        zpred[(((ij-1)*(npred))+1):(ij*npred),] = apply(matrix(z_true[ij,],ncol=q,nrow=1),2,rep,npred)
      }
    }   
    XZ_pred= cbind(xpred, zpred)
    pred= EzGP_predict(XZ_pred[,1:(p+q)], model, MSE_on = 1)
    ypred = pred$Y_hat
    mse = pred$MSE
    for(ii in 1:length(mse)){ if(mse[ii]<0){
      mse[ii]=0
    }
    }
#------------------------------------------------------
#ARSD
#------------------------------------------------------
    ARSD = ARSD_f(ypred,mse,a)
    LCB_upperbound=min(ARSD$UCB)
    A_n=which(ARSD$LCB<=LCB_upperbound)
    XZ_pred_An=XZ_pred[A_n,]
    if(is.matrix(XZ_pred_An)!=TRUE) XZ_pred_An = matrix(XZ_pred_An,1,ncol=length(XZ_pred_An))
#------------------------------------------------------
#Criterion selection
#------------------------------------------------------ 
    start_time <- Sys.time()  
    LCB_An = LCB_f(XZ_pred_An[,1:(p+q)],model,a)
    XZ_new= matrix(XZ_pred_An[which.min(LCB_An),],nrow=1)
    end_time <- Sys.time()
    time[i]=as.numeric(end_time - start_time,unit = "secs")
#------------------------------------------------------
##add next point  
#------------------------------------------------------        
   Y_new=computer_simulator(XZ_new)
   new_input=cbind(XZ_new,Y_new)
#------------------------------------------------------
##update
#------------------------------------------------------
    n0_dat=rbind(n0_dat,new_input)
    model = EzGP_fit(n0_dat[,1:(p+q)],n0_dat[,(p+q+1)],p=p,q=q,m=m,tau = tau) 
#-------------------------------------------------------------------------------
#Measurements
#-------------------------------------------------------------------------------    
#min|y-a|    
    y_a[i]=n0_dat[nrow(n0_dat),(p+q+1)]
    
    # #predictions  
    pred= EzGP_predict(C_t_XZ, model, MSE_on = 0)
    Yhat = pred$Y_hat
    #M_c0
    M_c0_mean[i]=(1/nrow(M_c0_M))*(sum(abs(C_t_Y-Yhat)))
  
    
  }
  
  #-------------------------------------------------------------------------------
  #save results
  #-------------------------------------------------------------------------------   
  filename= paste0("M_c0_mean-n",n0,"_","N",N,"_",method_name,"_",round(a,digits = 3),"_",l,".rds",sep="")
  saveRDS(M_c0_mean,file=filename)
  filename= paste0("y_a-n",n0,"_","N",N,"_",method_name,"_",round(a,digits = 3),"_",l,".rds",sep="")
  saveRDS(y_a,file=filename)
  filename= paste0("time-n",n0,"_","N",N,"_",method_name,"_",round(a,digits = 3),"_",l,".rds",sep="")
  saveRDS(time,file=filename)
}

l=as.integer(commandArgs(trailingOnly = TRUE))
f(l)

