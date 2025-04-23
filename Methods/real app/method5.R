rm(list=ls()) 
#-------------------------------------------------------------------------------
#Packages
#-------------------------------------------------------------------------------
library(lhs)
library(EzGP)
library(laGP)
library(dplyr)
#serial job
f=function(l){
#-------------------------------------------------------------------------------
#initial data
#-------------------------------------------------------------------------------
  method_name="ARSD-LCB"
  print(method_name)
  load("initial.RData")
#-----------------------------------------------------------------------------
#real data
data_HPC_scaled=as.matrix(read.csv("data/data_HPC_scaled.csv")[,2:(p+q+3)])   
#-------------------------------------------------------------------------------
#measurements and time
#-------------------------------------------------------------------------------
  #M_c0
  M_c0_mean=c()
  #min|y-a|
  y_a=c()
  #time
  time=c()
  #count group
  count_RCC=c()
#-------------------------------------------------------------------------------
#n0_dat
#-------------------------------------------------------------------------------
  filename= paste0("data/n0_dat","_",l,".csv",sep="")
  n0_dat=as.matrix(read.csv(file = filename))[,(2:(p+q+3))]
  model = EzGP_fit(n0_dat[,1:(p+q)],n0_dat[,(p+q+1)],p=p,q=q,m=m,tau=tau)
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
    # ypred_min=ypred[which.min(LCB)]
    # mse_min=mse[which.min(LCB)]
    # est=list(LCB=LCB,ypred_min=ypred_min,mse_min=mse_min)
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
    XZ_pred=as.matrix(data_HPC_scaled[-n0_dat[,p+q+2],])
    id_HPC=c()
    W_dat=matrix(NA, ncol =p+q,nrow=npred*prod(m))
    x00=randomLHS(npred,p)
    W_dat[,1:p]=matrix(x00,prod(m)*npred, ncol(x00) ,byrow=T)
    W_dat[,(p+1):(p+q)]=(rep(1:m,each=npred))
    colnames(W_dat)=c("x1","x2","x3","x4","z")
    for(io in 1:prod(m)){
      W_LHD_z_i=as.matrix(data.frame(W_dat) %>% filter(z ==io ))
      data_z_i=as.matrix(data.frame(XZ_pred) %>% filter(Test ==io ))
      for(ip in 1:(table(W_dat[,(p+1):(p+q)])[io])){
        #find closest distance between W_LHD and data_z_i
        dist_init=distance(data_z_i[,1:p],matrix(W_LHD_z_i[ip,1:p],nrow=1))
        id_z=which.min(dist_init)
        id=data_z_i[id_z,p+q+2]
        data_z_i=data_z_i[-id_z,]
        id_HPC=c(id_HPC,id)
      }
    }
    
    XZ_pred=data_HPC_scaled[id_HPC,]
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
##update
#------------------------------------------------------
    n0_dat=rbind(n0_dat,XZ_new)
   #print(n0_dat)
     model = EzGP_fit(n0_dat[,1:(p+q)],n0_dat[,(p+q+1)],p=p,q=q,m=m,tau = tau) 
    
#-------------------------------------------------------------------------------
#Measurements
#-------------------------------------------------------------------------------   
#min|y-a|    
    y_a[i]=n0_dat[nrow(n0_dat),(p+q+1)]
     C_t_XZ=data_HPC_scaled[which(data_HPC_scaled[,1+p+q]>=a-epsilon & data_HPC_scaled[,1+p+q]<=a+epsilon),]
     if(is.matrix(C_t_XZ)!=TRUE) C_t_XZ = matrix(C_t_XZ,1,ncol=length(C_t_XZ))
     C_t_Y=as.numeric(C_t_XZ[,1+p+q])
     
    
    pred= EzGP_predict(C_t_XZ[,1:(p+q)], model, MSE_on = 1)
    Yhat = pred$Y_hat

    #M_c0
    M_c0_mean[i]=(1/nrow(C_t_XZ))*(sum(abs(C_t_Y-Yhat)))
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
