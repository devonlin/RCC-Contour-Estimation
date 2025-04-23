

P=4
nmethod=7
#one-shot
one_shot=df_M_c0[(6*P+1):(7*P),2]
one_shot=round(one_shot,digit=4)
one_shot

#ARSD
ARSD=df_M_c0[(4*P+1):(5*P),2]
ARSD=round(ARSD,digit=4)
R5=round(one_shot/ARSD,digits = 2)

#RCC
RCC=df_M_c0[(1*P+1):(2*P),2]
RCC=round(RCC,digit=4)
R2=round(one_shot/RCC,digits=2)


#RCC_EI
RCC_EI=df_M_c0[(1):P,2]
RCC_EI=round(RCC_EI,digit=4)
R1=round(one_shot/RCC_EI,digits=2)

#ECL
ECL=df_M_c0[(2*P+1):(3*P),2]
ECL=round(ECL,digit=4)
R3=round(one_shot/ECL,digits=2)
R3

#EI
EI=df_M_c0[(3*P+1):(4*P),2]
EI=round(EI,digit=4)
R4=round(one_shot/EI,digits=2)

#LCB
LCB=df_M_c0[(5*P+1):(6*P),2]
LCB=round(LCB,digit=4)
R6=round(one_shot/LCB,digits=2)




M_c0_value=matrix(NA,ncol=7,nrow=P)
colnames(M_c0_value)=c("oneshot","ARSD","RCC","RCC_EI","ECL","EI","LCB")
M_c0_value[,1]=one_shot
M_c0_value[,2]=ARSD
M_c0_value[,3]=RCC
M_c0_value[,4]=RCC_EI
M_c0_value[,5]=ECL
M_c0_value[,6]=EI
M_c0_value[,7]=LCB
M_c0_value





relative_eff=matrix(NA,ncol=6,nrow=P)
colnames(relative_eff)=c("ARSD","RCC","RCC_EI","ECL","EI","LCB")
relative_eff[,1]=R5
relative_eff[,2]=R2
relative_eff[,3]=R1
relative_eff[,4]=R3
relative_eff[,5]=R4
relative_eff[,6]=R6

M_c0_table=cbind(M_c0_value[,1:2],relative_eff[,1],M_c0_value[,3],
                 relative_eff[,2],M_c0_value[,4],relative_eff[,3],M_c0_value[,5],relative_eff[,4],
                 M_c0_value[,6],relative_eff[,5],M_c0_value[,7],relative_eff[,6])
colnames(M_c0_table)=c("oneshot","ARSD","R5","RCC","R2","RCC_EI","R1","ECL","R3","EI","R4","LCB","R6")

library(xtable)
xtable(M_c0_table,digits = 4)
M_c0_table


