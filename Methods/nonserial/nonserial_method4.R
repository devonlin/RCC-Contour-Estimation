rm(list=ls()) 
#-------------------------------------------------------------------------------
#Packages
#-------------------------------------------------------------------------------
library(lhs)
library(EzGP)
library(laGP)
#-------------------------------------------------------------------------------
#initial data
#-------------------------------------------------------------------------------
  method_name="EI"
  print(method_name)
  load("initial.RData")
#-------------------------------------------------------------------------------
#measurements and time
#-------------------------------------------------------------------------------
#M_c0 
  M_c0_mean=matrix(NA,nrow=nnew,ncol=nsim)
  #min|y-a|
  y_a=matrix(NA,nrow=nnew,ncol=nsim)
  #time
  time=matrix(NA,nrow=nnew,ncol=nsim)
#-------------------------------------------------------------------------------
#simulations
#-------------------------------------------------------------------------------
for(l in 1:nsim){
print(l)#print nsim
#-------------------------------------------------------------------------------
 #n0 data read same initial for each method
#-------------------------------------------------------------------------------
filename= paste0("data/tradata","_",l,".csv",sep="")
n0_dat=as.matrix(read.csv(file = filename))[,(2:(p+q+2))]
model = EzGP_fit(n0_dat[,1:(p+q)],n0_dat[,(p+q+1)],p=p,q=q,m=m,tau=tau)
#-----------------------------------------------------------------------------
#Criterion functions
#-----------------------------------------------------------------------------
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
      EI_optim = t1+t2+t3
      return(-EI_optim)
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
#------------------------------------------------------
#Criterion selection
#------------------------------------------------------    
start_time <- Sys.time()  
EI = EI_f(XZ_pred[,1:(p+q)],model,a)
XZ_new = matrix(XZ_pred[which.min(EI),],nrow=1)
end_time <- Sys.time()
time[i,l]=as.numeric(end_time - start_time,unit = "secs")
#------------------------------------------------------
#add next point  
#------------------------------------------------------        
Y_new=computer_simulator(XZ_new)
new_input=cbind(XZ_new,Y_new)
#------------------------------------------------------
#update
#------------------------------------------------------
n0_dat=rbind(n0_dat,new_input)
model = EzGP_fit(n0_dat[,1:(p+q)],n0_dat[,(p+q+1)],p=p,q=q,m=m,tau = tau) 
#-------------------------------------------------------------------------------
#Measurements
#-------------------------------------------------------------------------------    
#min|y-a|    
    y_a[i,l]=n0_dat[nrow(n0_dat),(p+q+1)]
    
    #predictions  
    pred= EzGP_predict(C_t_XZ, model, MSE_on = 0)
    Yhat = pred$Y_hat
    #M_c0
    M_c0_mean[i,l]=(1/nrow(M_c0_M))*(sum(abs(C_t_Y-Yhat)))
    
  }
  } 
#-------------------------------------------------------------------------------
#save results
#-------------------------------------------------------------------------------   
filename= paste0("M_c0_mean-n",n0,"_","N",N,"_",method_name,"_",round(a,digits = 3),".rds",sep="")
saveRDS(M_c0_mean,file=filename)
filename= paste0("y_a-n",n0,"_","N",N,"_",method_name,"_",round(a,digits = 3),".rds",sep="")
saveRDS(y_a,file=filename)
filename= paste0("time-n",n0,"_","N",N,"_",method_name,"_",round(a,digits = 3),".rds",sep="")
saveRDS(time,file=filename)