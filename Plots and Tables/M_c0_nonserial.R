
rm(list=ls())
nsim=50
n0=9
N=N_total=24
a=0.5
nmethod=7
P=N-n0

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
method_name="RCC-EI"
method_1=matrix(NA,nrow=P,ncol=nsim)
  filename= paste0("M_c0_mean-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3),".rds",sep="")
  #filename= paste0("time-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3) ,".rds",sep="")
  method_1=readRDS(filename)[,1:nsim]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#2
method_name="RCC"
method_2=matrix(NA,nrow=P,ncol=nsim)
  filename= paste0("M_c0_mean-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3),".rds",sep="")
  #filename= paste0("time-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3) ,".rds",sep="")
  method_2=readRDS(filename)[,1:nsim]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#3
method_name="ECL"
method_3=matrix(NA,nrow=P,ncol=nsim)
  filename= paste0("M_c0_mean-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3) ,".rds",sep="")
 #filename= paste0("time-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3) ,".rds",sep="")
  method_3=readRDS(filename)[,1:nsim]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#4
method_name="EI"
method_4=matrix(NA,nrow=P,ncol=nsim)
  filename= paste0("M_c0_mean-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3) ,".rds",sep="")
  #filename= paste0("time-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3) ,".rds",sep="")
  method_4=readRDS(filename)[,1:nsim]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#5
method_name="ARSD-LCB"
method_5=matrix(NA,nrow=P,ncol=nsim)
  filename= paste0("M_c0_mean-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3) ,".rds",sep="")
  #filename= paste0("time-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3) ,".rds",sep="")
  method_5=readRDS(filename)[,1:nsim]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

method_name="LCB"
method_6=matrix(NA,nrow=P,ncol=nsim)
  filename= paste0("M_c0_mean-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3) ,".rds",sep="")
  #filename= paste0("time-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3) ,".rds",sep="")
  method_6=readRDS(filename)[,1:nsim]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

method_name="one-shot"
method_7=matrix(NA,nrow=P,ncol=nsim)
filename= paste0("M_c0_mean-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3) ,".rds",sep="")
#filename= paste0("time-n",n0,"_","N",N_total,"_",method_name,"_",round(a,digits = 3) ,".rds",sep="")
method_7=readRDS(filename)[,1:nsim]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#M_c0
#Ex1
#N_seq=c(3,6,9,12)
#N_seq=c(2,4,6,8,10)

methods=rep(c("RCC-EI","RCC","ECL","EI","ARSD-LCB","LCB","one-shot"), each=length(N_seq))

kkk=1:50
nsim=50
P22=length(N_seq)
P2=N_seq
M_c0_mean_mat=array(NA,dim=c(P22,nsim,nmethod))
M_c0_mean_mat[,,1]=method_1[P2,kkk]#RCC-EI-dis
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

df_M_c0 = data.frame(samplesize = rep(N_seq,nmethod), M_c0= M_c0_mean_center, Methods=methods)


filename= paste0("M_c0_final","_","n",n0,"_","N",N_total,"_",round(a,digits = 3),".Rdata",sep="")
save.image(filename)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#time
ppp=12
methods=c("RCC-EI","RCC","ECL","EI","ARSD-LCB","LCB")
time_secs=round(colMeans(cbind(colSums(method_1[1:ppp,1:nsim]),colSums(method_2[1:ppp,1:nsim]),colSums(method_3[1:ppp,1:nsim]),
                               colSums(method_4[1:ppp,1:nsim]),colSums(method_5[1:ppp,1:nsim]),colSums(method_6[1:ppp,1:nsim]))),digits = 2)
names(time_secs)=methods
time_secs
filename= paste0("time","_","n",n0,"_","N",N_total,"_",round(a,digits = 3),".Rdata",sep="")
save.image(filename)
