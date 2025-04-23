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
method_name="RCC-EI"
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
#RCC
RCC_f <- function(ypred,mse,a){
  alpha=0.05
  mu=abs(ypred-a)
  beta_h=2*log(((pi^2)*((nrow(n0_dat))^2)*(prod(m)))/(6*(alpha)))
  LCB=(mu)-(sqrt(beta_h))*(sqrt(mse))
  UCB=(mu)+(sqrt(beta_h))*(sqrt(mse))
  est=list(LCB=LCB,UCB=UCB)
  return(est)
}

SD_f=function(x,model,a){
  if(is.matrix(x)!=TRUE) x = matrix(x,nrow=1,ncol=p+q)
  pred= EzGP_predict(x[,1:(p+q)], model, MSE_on = 1)
  ypred = pred$Y_hat
  mse = pred$MSE
  for(ii in 1:length(mse)){ if(mse[ii]<0){
    mse[ii]=0
  }
  }
  SD=which.max(sqrt(mse))
  ypred_max=ypred[SD]
  mse_max=mse[SD]
  return(est=list(SD_id=SD,ypred_max=ypred_max,mse_max=mse_max))
}

  #EI
  EI_f=function(x,model,a){
    if(is.matrix(x)!=TRUE) x = matrix(x,nrow=1,ncol=p+q)
    pred= EzGP_predict(x[,1:(p+q)], model, MSE_on = 1)
    ypred = pred$Y_hat
    mse = pred$MSE
    for(ii in 1:length(mse)){ if(mse[ii]<0){
      mse[ii]=0
    }
    }
    shat = sqrt(mse)
    ep = 1.96*shat
    u1 = (a-ypred-ep)/shat
    u2 = (a-ypred+ep)/shat
    t1 = (ep^2-(ypred-a)^2-mse)*(pnorm(u2)-pnorm(u1))
    t2 = mse*(u2*dnorm(u2)-u1*dnorm(u1))
    t3 = 2*(ypred-a)*shat*(dnorm(u2)-dnorm(u1))
    EI = which.min(-(t1+t2+t3))
    ypred_max=ypred[EI]
    mse_max=mse[EI]
      return(est=list(EI_id=EI,ypred_max=ypred_max,mse_max=mse_max))
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
#RCC
#------------------------------------------------------
    RCC = RCC_f(ypred,mse,a)
    A_1=which( 0<RCC$LCB)
    A_2=which(0>=RCC$LCB )
    LCB_upperbound=min(RCC$UCB)
    
    XZ_pred_A1=XZ_pred[A_1,]
    if(is.matrix(XZ_pred_A1)!=TRUE) XZ_pred_A1 = matrix(XZ_pred_A1,1,ncol=length(XZ_pred_A1))
    XZ_pred_A2=XZ_pred[A_2,]
    if(is.matrix(XZ_pred_A2)!=TRUE) XZ_pred_A2 = matrix(XZ_pred_A2,1,ncol=length(XZ_pred_A2))
    
    
#A1min criterion
    ypred_A1 = ypred[A_1]
    mse_A1 = mse[A_1]
    A_min=which(RCC$LCB[A_1]<=LCB_upperbound)
    XZ_pred_Amin=XZ_pred_A1[A_min,]
    if(is.matrix(XZ_pred_Amin)!=TRUE) XZ_pred_Amin= matrix(XZ_pred_Amin,1,ncol=length(XZ_pred_Amin))

#select next point
    start_time <- Sys.time()
    if(length(A_min)==0|| length(A_2)==0){
      if(length(A_min)==0){
        EI_A2 = EI_f(XZ_pred_A2[,1:(p+q)],model,a)
        XZ_new = matrix(XZ_pred_A2[EI_A2$EI_id,],nrow=1)
        region_ID=2
      }
      if(length(A_2)==0){
        SD_Amin = SD_f(XZ_pred_Amin,model,a)
        XZ_new = matrix(XZ_pred_Amin[SD_Amin$SD_id,],nrow=1)
        region_ID=1
      }
    } else{
      #A1min
      SD_Amin = SD_f(XZ_pred_Amin[,1:(p+q)],model,a)
      x_new_Amin = matrix(c(XZ_pred_Amin[SD_Amin$SD_id,],SD_Amin$ypred_max,SD_Amin$mse_max),nrow=1)
      ratio_Amin=matrix(c((sqrt(x_new_Amin[,(p+q+4)])/max(abs(x_new_Amin[,(p+q+3)]-a),del)),1),nrow=1)
      
      #A2
      EI_A2 = EI_f(XZ_pred_A2[,1:(p+q)],model,a)
      x_new_A2 = matrix(c(XZ_pred_A2[EI_A2$EI_id,],EI_A2$ypred_max,EI_A2$mse_max),nrow=1)
      ratio_A2=matrix(c((sqrt(x_new_A2[,(p+q+4)])/max(abs(x_new_A2[,(p+q+3)]-a),del)),2),nrow=1)
      
      x_new_can=rbind(cbind(x_new_Amin,ratio_Amin),cbind(x_new_A2,ratio_A2))
      select_ID=x_new_can[which.max(x_new_can[,(p+q+5)]),1:(p+q+2)]
      region_ID=x_new_can[which.max(x_new_can[,(p+q+5)]),(p+q+6)]
      XZ_new = matrix(select_ID,nrow=1)
    }
    end_time <- Sys.time()
    time[i]=as.numeric(end_time - start_time,unit = "secs")
    count_RCC[i]=region_ID
#------------------------------------------------------
##update
#------------------------------------------------------
    n0_dat=rbind(n0_dat,XZ_new)
    model = EzGP_fit(n0_dat[,1:(p+q)],n0_dat[,(p+q+1)],p=p,q=q,m=m,tau = tau) 
    
#-------------------------------------------------------------------------------
#Measurements
#-------------------------------------------------------------------------------   
#min|y-a|    
  y_a[i]=n0_dat[nrow(n0_dat),(p+q+1)]
    
    #predictions  
    C_t_XZ=data_HPC_scaled[which(data_HPC_scaled[,1+p+q]>=a-epsilon & data_HPC_scaled[,1+p+q]<=a+epsilon),]
    if(is.matrix(C_t_XZ)!=TRUE) C_t_XZ = matrix(C_t_XZ,1,ncol=length(C_t_XZ))
    C_t_Y=as.numeric(C_t_XZ[,1+p+q])
    
    pred= EzGP_predict(C_t_XZ[,1:(p+q)], model, MSE_on = 0)
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
  filename= paste0("count_RCC-n",n0,"_","N",N,"_",method_name,"_",round(a,digits = 3),"_",l,".rds",sep="")
  saveRDS(count_RCC,file=filename)
}

l=as.integer(commandArgs(trailingOnly = TRUE))
f(l)
