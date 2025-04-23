

rm(list=ls())
nsim=30
n0=9
N_total=54
a=2.6
nmethod=7
P=N_total-n0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
method_name="RCC-EI"
method_1=matrix(NA,nrow=P,ncol=nsim)
for(i in 1:nsim){
  filename= paste0("M_c0_mean-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3),"_",i,".rds",sep="")
 #filename= paste0("time-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3),"_",i,".rds",sep="")
  output=readRDS(filename)
  method_1[,i]=output
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#2
method_name="RCC"
method_2=matrix(NA,nrow=P,ncol=nsim)
for(i in 1:nsim){
  filename= paste0("M_c0_mean-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3),"_",i,".rds",sep="")
#filename= paste0("time-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3),"_",i,".rds",sep="")
  output=readRDS(filename)
  method_2[,i]=output
}


#3
method_name="ECL"
method_3=matrix(NA,nrow=P,ncol=nsim)
for(i in 1:nsim){
  filename= paste0("M_c0_mean-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3),"_",i,".rds",sep="")
  #filename= paste0("time-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3),"_",i,".rds",sep="")
  output=readRDS(filename)
  method_3[,i]=output
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#4
method_name="EI"
method_4=matrix(NA,nrow=P,ncol=nsim)
for(i in 1:nsim){
  filename= paste0("M_c0_mean-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3),"_",i,".rds",sep="")
 # filename= paste0("time-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3),"_",i,".rds",sep="")
  output=readRDS(filename)
  method_4[,i]=output
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#5
method_name="ARSD-LCB"
method_5=matrix(NA,nrow=P,ncol=nsim)
for(i in 1:nsim){
  filename= paste0("M_c0_mean-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3),"_",i,".rds",sep="")
  #filename= paste0("time-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3),"_",i,".rds",sep="")
  output=readRDS(filename)
  method_5[,i]=output
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

method_name="LCB"
method_6=matrix(NA,nrow=P,ncol=nsim)
for(i in 1:nsim){
  filename= paste0("M_c0_mean-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3),"_",i,".rds",sep="")
 #filename= paste0("time-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3),"_",i,".rds",sep="")
  output=readRDS(filename)
  method_6[,i]=output
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

method_name="one-shot"
method_7=matrix(NA,nrow=P,ncol=nsim)
for(i in 1:nsim){
  filename= paste0("M_c0_mean-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3),"_",i,".rds",sep="")
  #filename= paste0("time-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3),"_",i,".rds",sep="")
  output=readRDS(filename)
  method_7[,i]=output
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#M_c0
#Ex2 Ex3
#N=c(18,27,36,45)
#real data
N=c(105,110,115,120)

methods=rep(c("RCC-EI","RCC","ECL","EI","ARSD-LCB","LCB","one-shot"), each=length(N))
kkk=1:50
nsim=50
P22=length(N)
P2=N
M_c0_mean_mat=array(NA,dim=c(P22,nsim,nmethod))
M_c0_mean_mat[,,1]=method_1[P2,kkk]#RCC-EI
M_c0_mean_mat[,,2]=method_2[P2,kkk]#RCC
M_c0_mean_mat[,,3]=method_3[P2,kkk]#ECL
M_c0_mean_mat[,,4]=method_4[P2,kkk]#EI
M_c0_mean_mat[,,5]=method_5[P2,kkk]#ARSD-LCB
M_c0_mean_mat[,,6]=method_6[P2,kkk]#LCB
M_c0_mean_mat[,,7]=method_7[P2,kkk]#oneshot

M_c0_mean_center=NULL
for(i in 1:nmethod){
  for(j in 1:length(P2)){
    M_c0_mean_center=  c(M_c0_mean_center, mean(M_c0_mean_mat[j,,i]))
  }
}

df_M_c0 = data.frame(samplesize = rep(N,nmethod), M_c0= M_c0_mean_center, Methods=methods)


filename= paste0("M_c0_final","_","n",n0,"_","N",N_total,"_",round(a,digits = 3),".Rdata",sep="")
save.image(filename)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#time
PPP=120
methods=c("RCC-EI","RCC","ECL","EI","ARSD-LCB","LCB")
time_secs=round(colMeans(cbind(colSums(method_1[1:ppp,1:nsim]),colSums(method_2[1:PPP,]),colSums(method_3[1:PPP,]),
                               colSums(method_4[1:PPP,]), colSums(method_5[1:PPP,]),colSums(method_6[1:PPP,]))),digits = 2)
names(time_secs)=methods
time_secs
filename= paste0("time","_","n",n0,"_","N",N_total,"_",round(a,digits = 3),".Rdata",sep="")
save.image(filename)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



