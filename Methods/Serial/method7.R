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
method_name="one-shot"
print(method_name)
load("initial.RData")
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
#data
#-------------------------------------------------------------------------------
x=randomLHS(N,p)
if(q==1){ 
  z_true=matrix(1:m,ncol=q)
  z=matrix(apply(z_true,2,rep,N)[1:N,],ncol=q)
}else{   
  z=apply(z_true,2,rep,N)[1:N,]
}   
id_z=sample(1:N,N)
z=z[id_z,]
XZ= cbind(x, z)
y= computer_simulator(XZ)
n0_dat=cbind(XZ,y)
#-------------------------------------------------------------------------------
#Measurements
#-------------------------------------------------------------------------------    
#min|y-a|     
y_a=n0_dat[,(p+q+1)]
  for(i in 1:nnew){
    print(i)
    model = EzGP_fit(n0_dat[1:(n0+i),1:(p+q)],n0_dat[1:(n0+i),(p+q+1)],p=p,q=q,m=m,tau = tau) 
    
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

