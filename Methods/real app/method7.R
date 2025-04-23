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
method_name="one-shot"
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
id_HPC=c()
W_dat=matrix(NA, ncol =p+q,nrow=N)
  x00=randomLHS(N/(prod(m)),p)
  W_dat[,1:p]=matrix(x00,N, ncol(x00) ,byrow=T)
  W_dat[,(p+1):(p+q)]=(rep(1:m,each=N/(prod(m))))
  colnames(W_dat)=c("x1","x2","x3","x4","z")
  for(io in 1:prod(m)){
    W_LHD_z_i=as.matrix(data.frame(W_dat) %>% filter(z ==io ))
    data_z_i=as.matrix(data.frame(data_HPC_scaled) %>% filter(Test ==io ))
    for(ip in 1:(table(W_dat[,(p+1):(p+q)])[io])){
      #find closest distance between W_dat and data_z_i
      dist_init=distance(data_z_i[,1:p],matrix(W_LHD_z_i[ip,1:p],nrow=1))
      id_z=which.min(dist_init)
      id=data_z_i[id_z,p+q+2]
      data_z_i=data_z_i[-id_z,]
      id_HPC=c(id_HPC,id)
    }
  }
  n0_dat=as.matrix(data_HPC_scaled[id_HPC,])
  Id_sample_N=sample(1:N)
  n0_dat=n0_dat[Id_sample_N,]

#-------------------------------------------------------------------------------
#Measurements
#-------------------------------------------------------------------------------    
#min|y-a|     
  y_a=n0_dat[,(p+q+1)]
  for(i in 1:nnew){
    print(i)
    model = EzGP_fit(n0_dat[1:(n0+i),1:(p+q)],n0_dat[1:(n0+i),(p+q+1)],p=p,q=q,m=m,tau = tau) 
    #-------------------------------------------------------------------------------
    #Measurements
    #-------------------------------------------------------------------------------   
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
}

l=as.integer(commandArgs(trailingOnly = TRUE))
f(l)
