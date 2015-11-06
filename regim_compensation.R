
source("E:/USU_Research_work/TOPNET PROJECT/MODEL COMPARISON/From_Sulochan/Sixteen_variables/Functions for R/getData.r")
source("E:/USU_Research_work/TOPNET PROJECT/MODEL COMPARISON/From_Sulochan/Sixteen_variables/Functions for R/WY_conv.R")
source("E:/USU_Research_work/TOPNET PROJECT/MODEL COMPARISON/From_Sulochan/Sixteen_variables/Functions for R/PCM.R")
source("E:/USU_Research_work/TOPNET PROJECT/MODEL COMPARISON/From_Sulochan/Sixteen_variables/Functions for R/Q7MaxMin.R")
source("E:/USU_Research_work/TOPNET PROJECT/MODEL COMPARISON/From_Sulochan/Sixteen_variables/Functions for R/timings.R")
source("E:/USU_Research_work/TOPNET PROJECT/MODEL COMPARISON/From_Sulochan/Sixteen_variables/Functions for R/TPTH.R")
source("E:/USU_Research_work/TOPNET PROJECT/MODEL COMPARISON/From_Sulochan/Sixteen_variables/Functions for R/flowrev.R")
source("E:/USU_Research_work/TOPNET PROJECT/MODEL COMPARISON/From_Sulochan/Sixteen_variables/Functions for R/BFI_whole.R")
source("E:/USU_Research_work/TOPNET PROJECT/MODEL COMPARISON/From_Sulochan/Sixteen_variables/Functions for R/COVyr.R")
source("E:/USU_Research_work/TOPNET PROJECT/MODEL COMPARISON/From_Sulochan/Sixteen_variables/Functions for R/MEANyr.R")
source("E:/USU_Research_work/TOPNET PROJECT/MODEL COMPARISON/From_Sulochan/Sixteen_variables/Functions for R/zeroday.R")
source("E:/USU_Research_work/TOPNET PROJECT/MODEL COMPARISON/From_Sulochan/Sixteen_variables/Functions for R/q167.R")
source("E:/USU_Research_work/TOPNET PROJECT/MODEL COMPARISON/From_Sulochan/Sixteen_variables/Functions for R/q167_lp3.R")
source("E:/USU_Research_work/TOPNET PROJECT/MODEL COMPARISON/From_Sulochan/Sixteen_variables/Functions for R/COV.R")
source("E:/USU_Research_work/TOPNET PROJECT/MODEL COMPARISON/From_Sulochan/Sixteen_variables/Functions for R/qmean_whole.R")
source("E:/USU_Research_work/TOPNET PROJECT/MODEL COMPARISON/From_Sulochan/Sixteen_variables/Functions for R/qmean_whole.R")
source("E:/USU_Research_work/TOPNET PROJECT/MODEL COMPARISON/From_Sulochan/Sixteen_variables/Functions for R/SI.R")
source("E:/USU_Research_work/TOPNET PROJECT/MODEL COMPARISON/From_Sulochan/Sixteen_variables/Functions for R/flw_puls_evnt_5_95.R")

library(hydroTSM)
library(ggplot2)
library(hydroGOF)
library(reshape2)
library(plyr)
###RUN MODEL ###


watsed=c("A1","A2","B1","B21","C1","C22")
stream_source=c("Daymet_RUN","PRISM_RUN","Maurer_RUN")
# for ( i in 1:length(watsed)){
#   for ( j in 1:length(stream_source)) {
#     
#     setwd(paste("E:\\USU_Research_work\\RAINFALL_COMPARISON\\",watsed[i],"_watershed\\",stream_source[j],sep=""))
#     system(paste("topnet_modified"))
#   
#   }
# }


banum=c(29,113,25,135,131,43)
regime_function=function(flow){
  stnvar=matrix(NA,nrow=1,ncol=58)
  col_name=c("BFI","Slp_BFI","Pval_BFI",
             "CoV_whole",
             "CoV_yr","Slp_CoV","Pval_CoV",
             "Qmean_whole",
             "Qmean_yr","Slp_Qmean","Pval_Qmean",
             "Zero","Slp_Zero","Pval_Zero",
             "BFF_Ln",
             "FD_Ln","Slp_FD_Ln","Pval_FD_Ln",
             "BFF_Lp3",
             "FD_Lp3","Slp_FD_Lp3","Pval_FD_Lp3",
             "P","C","M",
             "T25","Slp_T25","Pval_T25",
             "T50","Slp_T50","Pval_T50",
             "T75","Slp_T75","Pval_T75",
             "Pk_time","H1",
             "Q7min","Slp_Q7min","Pval_Q7min",
             "Q7max","Slp_Q7max","Pval_Q7max",
             "FR","Slp_FR","Pval_FR",
             "SI","Slp_SI","Pval_SI",
             "HFE","Slp_HFE","Pval_HFE",
             "LFE","Slp_LFE","Pval_LFE",
             "ZFE","Slp_ZFE","Pval_ZFE",
             "Stn_No")
  colnames(stnvar)=col_name
  #flow=obsflow
  WYmat=WY_conv.R(flow) 
  
  timing=t25t50t75(WYmat)
  t25=as.numeric(timing[1])
  t50=as.numeric(timing[4])
  t75=as.numeric(timing[7])
  slope_t25=as.numeric(timing[2])
  slope_t50=as.numeric(timing[5])
  slope_t75=as.numeric(timing[8])
  pval_t25=as.numeric(timing[3])
  pval_t50=as.numeric(timing[6])
  pval_t75=as.numeric(timing[9])
  # if matrix of T25, T50 and T75 is needed:
  # matT25=unlist(timing[10])
  # matT50=unlist(timing[11])
  # matT75=unlist(timing[12])
  stnvar[26]=t25
  stnvar[27]=slope_t25
  stnvar[28]=pval_t25
  stnvar[29]=t50
  stnvar[30]=slope_t50
  stnvar[31]=pval_t50
  stnvar[32]=t75
  stnvar[33]=slope_t75
  stnvar[34]=pval_t75
  
  pkhar1=tpth(WYmat)
  pktime=as.numeric(pkhar1[1])
  har1=as.numeric(pkhar1[2])
  stnvar[35]=pktime
  stnvar[36]=har1
  
  max7min7=Q7MaxMin.R(WYmat)
  Q7min=as.numeric(max7min7[1])
  slope_Q7min=as.numeric(max7min7[2])
  pval_Q7min=as.numeric(max7min7[3])
  
  Q7max=as.numeric(max7min7[5])
  slope_Q7max=as.numeric(max7min7[6])
  pval_Q7max=as.numeric(max7min7[7])
  stnvar[37]=Q7min
  stnvar[38]=slope_Q7min
  stnvar[39]=pval_Q7min
  stnvar[40]=Q7max
  stnvar[41]=slope_Q7max
  stnvar[42]=pval_Q7max
  
  fr_list= flowrev(WYmat)
  fr=as.numeric(fr_list[1])
  slope_fr=as.numeric(fr_list[2])
  pval_fr=as.numeric(fr_list[3])
  # if matrix of yearly Flow Reversal is needed:
  # mat=unlist(fr_list[4])
  stnvar[43]=fr
  stnvar[44]=slope_fr
  stnvar[45]=pval_fr
  
  bfi_list=BFI_whole.R(WYmat)
  BFI=as.numeric((bfi_list)[1])
  slope_BFI=as.numeric((bfi_list)[2])
  pval_BFI=as.numeric((bfi_list)[3])
  stnvar[1]=BFI
  stnvar[2]=slope_BFI
  stnvar[3]=pval_BFI
  
  # if matrix of yearly BFI is needed:
  # mat=unlist(bfi_list[4])
  
  cov=COV.R(WYmat)
  stnvar[4]=cov  
  
  cov_list=COVyr.R(WYmat)
  cov_yr=as.numeric(cov_list[1])
  slope_cov=as.numeric(cov_list[2])
  pval_cov=as.numeric(cov_list[3])
  # if matrix of yearly CoV is needed:
  # mat=unlist(cov_list[4])
  stnvar[5]=cov_yr
  stnvar[6]=slope_cov
  stnvar[7]=pval_cov
  
  qmean1=qmean_whole.R(WYmat)
  stnvar[8]=qmean1
  
  qmean_list=MEANyr.R(WYmat)
  qmean_yr=as.numeric(qmean_list[1])
  slope_qmean=as.numeric(qmean_list[2])
  pval_qmean=as.numeric(qmean_list[3])
  # if matrix of yearly Qmean is needed:
  # mat=unlist(qmean_list[4])
  stnvar[9]=qmean_yr
  stnvar[10]=slope_qmean
  stnvar[11]=pval_qmean
  
  zero_list=zeroday.R(WYmat)
  zero_yr=as.numeric(zero_list[1])
  slope_zero=as.numeric(zero_list[2])
  pval_zero=as.numeric(zero_list[3])
  # if matrix of yearly Qmean is needed:
  # mat=unlist(zero_list[4])
  stnvar[12]=zero_yr
  stnvar[13]=slope_zero
  stnvar[14]=pval_zero
  
  bff_ln_list=q167.R(WYmat) ## bank Full Flow = Q1.67 using log-Normal
  bff_ln=as.numeric(bff_ln_list[1])
  flddur_ln=as.numeric(bff_ln_list[2]) ## Flood Duration = Q1.67 using Log-Normal
  slope_flddur=as.numeric(bff_ln_list[3])
  pval_flddur=as.numeric(bff_ln_list[4])
  # if matrix of yearly Flood Duration is needed:
  # mat=unlist(bff_ln_list[5])
  stnvar[15]=bff_ln
  stnvar[16]=flddur_ln
  stnvar[17]=slope_flddur
  stnvar[18]=pval_flddur
  
  bff_lp3_list=q167_lp3.R(WYmat) ## bank Full Flow = Q1.67 using log-Pearson III
  bff_lp3=as.numeric(bff_lp3_list[1])
  flddur_lp3=as.numeric(bff_lp3_list[2]) ## Flood Duration = Q1.67 using Log-Pearson III
  slope_flddur_lp3=as.numeric(bff_lp3_list[3])
  pval_flddur_lp3=as.numeric(bff_lp3_list[4])
  # if matrix of yearly Flood Duration is needed:
  # mat=unlist(bff_lp3_list[5])
  stnvar[19]=bff_lp3
  stnvar[20]=flddur_lp3
  stnvar[21]=slope_flddur_lp3
  stnvar[22]=pval_flddur_lp3
  
  pcm_list=PCM(flow)
  p=as.numeric(pcm_list[1])
  c=as.numeric(pcm_list[2])
  m=as.numeric(pcm_list[3])
  # if contingency table (bin matrix) is needed:
  # vec = unlist(pcm_list[4])
  # mat = matrix(vec,nrow=7,ncol=12)
  stnvar[23]=p
  stnvar[24]=c
  stnvar[25]=m
  
  SI_list=SI.R(WYmat)
  SI=as.numeric(SI_list[1])
  slope_SI=as.numeric(SI_list[2])
  pval_SI=as.numeric(SI_list[3])
  # if matrix of yearly Qmean is needed:
  # mat=unlist(SI_list[4])
  stnvar[46]=SI
  stnvar[47]=slope_SI
  stnvar[48]=pval_SI
  
  FE_list=flw_puls_evnt_5_95.R(WYmat)
  LFE=as.numeric(FE_list[1])
  HFE=as.numeric(FE_list[2])
  ZFE=as.numeric(FE_list[3])
  
  slope_LFE=as.numeric(FE_list[7])
  pval_LFE=as.numeric(FE_list[8])
  slope_HFE=as.numeric(FE_list[9])
  pval_HFE=as.numeric(FE_list[10])
  slope_ZFE=as.numeric(FE_list[11])
  pval_ZFE=as.numeric(FE_list[12])
  
  # if matrix of yearly HFE,LFE and ZFE is needed:
  # matLFE=unlist(FE_list[4])
  # matHFE=unlist(FE_list[5])
  # matZFE=unlist(FE_list[6])
  
  stnvar[49]=HFE
  stnvar[50]=slope_HFE
  stnvar[51]=pval_HFE
  
  stnvar[52]=LFE
  stnvar[53]=slope_HFE
  stnvar[54]=pval_HFE
  
  stnvar[55]=ZFE
  stnvar[56]=slope_ZFE
  stnvar[57]=pval_ZFE
  return(stnvar)
  
  
}

watsed=c("A1","A2","B1","B21","C1","C22")
stream_source=c("Daymet_RUN","PRISM_RUN","Maurer_RUN")
sr=c('DY','MR','PR')
SNVAR_ALL=matrix(NA,nrow=1,ncol=58)
CE_all=matrix(NA,nrow=16,ncol=1)

sixvar=matrix(NA,nrow=1,ncol=16)
CE=matrix(NA,nrow=16,ncol=1)
for ( j in 1:length(stream_source)){
    
    for ( k in 1:3) {
  
      SNVAR=matrix(NA,nrow=6,ncol=16)
      obsvar=matrix(NA,nrow=6,ncol=16)
      for ( i in 1:length(watsed)) {
       
    setwd(paste("F:\\USU_Research_work_update_March_30_2014\\RAINFALL_COMPARISON\\",watsed[i],"_watershed\\",stream_source[j],sep=""))
        
    targetdir <- paste("F:\\USU_Research_work_update_March_30_2014\\RAINFALL_COMPARISON\\",watsed[i],"_watershed\\",stream_source[j],sep="")
    
    origindir <- paste("F:\\USU_Research_work_update_March_30_2014\\RAINFALL_COMPARISON\\TOPINP\\",watsed[i],"Watershed\\",sr[k],sep="")
    file.copy(paste (origindir, "topinp.dat", sep = "/"), paste(targetdir, "topinp.dat" , sep = "/"), overwrite=TRUE)
    origindir1 <-"E:\\exePrograms"
    file.copy(paste (origindir1, "topnet_modified.exe", sep = "/"), paste(targetdir, "topnet_modified.exe" , sep = "/"), overwrite=TRUE,recursive=FALSE,copy.mode = TRUE)
    
    system(paste("topnet_modified"))
    
    time=seq(as.Date("2000/1/1"), as.Date("2010/12/31"), "day")
    time1=seq(as.Date("2001/10/1"), as.Date("2010/9/30"), "day")
    date_match=match(time1,time) # get overlap time interval
    
    
    sf=scan('streamflow_calibration.dat', what="")
    sf1=sf[seq(21,length(sf),1)] ## need to change theis things later
    obs_flow=(as.numeric(sf1[seq(1,length(sf1),3)])) ##need to change this later
    time_all_observ= seq(as.Date("1980/1/1"), as.Date("2012/12/31"), "day") #create time series
    date_overlap_obs=match(time1,time_all_observ) # get overlap time interval
    obsflow=data.frame(date=time1,q=obs_flow[date_overlap_obs])
    
    SN1=regime_function(obsflow)
    
    ff=scan('FlowAtStreamNodes_cms.txt', what="")
    l=banum[i]+2#for A2
    ff1=ff[seq(l,length(ff),1)] ## need to change theis things later
    simu_flow=matrix(as.numeric(ff1[seq(2,length(ff1),l-1)])) ##need to change this later
    simflow=data.frame(date=time1,q=simu_flow[date_match])
    simflow[simflow<0.001]=0.001
    
    SN2=regime_function(simflow)
    SN=rbind(SN1,SN2)
    SNVAR_ALL=rbind(SNVAR_ALL,SN)
    SN1=matrix(as.numeric(SN1[1,1:57]),nrow=1,ncol=57)
    SN2=matrix(as.numeric(SN2[1,1:57]),nrow=1,ncol=57)
    var_order=c(4,8,15,16,23,24,25,29,35,37,40,43,46,49,52,55)
    #obs_seq=seq(1)
    #mod_seq=seq(2,36,6)
    #sixteen_VAR=rv_proj1[,var_order];
    #observed=sixteen_VAR[1,]
    #model=sixteen_VAR[2,]
    obsvar[i,]=SN1[var_order]
    SNVAR[i,]=SN2[var_order]
   # sixvar=rbind(sixvar,SNVAR)
    }
    
   for ( m in 1:16){
     
     CE[m]=NSE(obsvar[,m],SNVAR[,m]) }
    
   CE_all=cbind(CE_all,CE)
   rm(obsvar)
   rm(SNVAR)
}
  
}


#write.csv(stnvar,file="A1_climate_change_regime.xls")

rv_proj1=matrix(as.numeric(SNVAR_ALL[2:37,1:57]),nrow=36,ncol=57)
var_order=c(4,8,15,16,23,24,25,29,35,37,40,43,46,49,52,55)
obs_seq=seq(1,36,6)
day_seq=seq(2,36,6)
prism_seq=seq(4,36,6)
maurer_seq=seq(6,36,6)

sixteen_VAR=rv_proj1[,var_order];
observed=sixteen_VAR[obs_seq,]
daymet=sixteen_VAR[day_seq,]
prism=sixteen_VAR[prism_seq,]
maurer=sixteen_VAR[maurer_seq,]

CE=matrix(NA,nrow=3,ncol=16)
for ( i in 1:16){
  
  CE[1,i]=NSE(daymet[,i],observed[,i])
  CE[2,i]=NSE(prism[,i],observed[,i])
  CE[3,i]=NSE(maurer[,i],observed[,i])
}

CE=t(CE)

porder=c(2,1,11,10,3,12,8,5,6,7,14,15,9,13,4,16)
obs1=observed[,porder];daymet1=daymet[,porder];prism1=prism[,porder];maurer1=maurer[,porder]


#########scatter plot to see how it looks like plot#########
col_name=c("BFI","Slp_BFI","Pval_BFI",
           "CoV_whole",
           "CoV_yr","Slp_CoV","Pval_CoV",
           "Qmean_whole",
           "Qmean_yr","Slp_Qmean","Pval_Qmean",
           "Zero","Slp_Zero","Pval_Zero",
           "BFF",
           "FLDD","Slp_FD_Ln","Pval_FD_Ln",
           "BFF_Lp3",
           "FD_Lp3","Slp_FD_Lp3","Pval_FD_Lp3",
           "P","C","M",
           "T25","Slp_T25","Pval_T25",
           "T50","Slp_T50","Pval_T50",
           "T75","Slp_T75","Pval_T75",
           "Pk_time","H1",
           "Q7min","Slp_Q7min","Pval_Q7min",
           "Q7max","Slp_Q7max","Pval_Q7max",
           "FR","Slp_FR","Pval_FR",
           "SI","Slp_SI","Pval_SI",
           "HFE","Slp_HFE","Pval_HFE",
           "LFE","Slp_LFE","Pval_LFE",
           "ZFE","Slp_ZFE","Pval_ZFE",
           "Stn_No")
col_names=matrix(col_name)
var_order=c(4,8,15,16,23,24,25,29,35,37,40,43,46,49,52,55)
names=col_names[var_order]


######serial based on performance:

sdnames=c("a1","a2","a3","a4","a5","a6","a7","a8","b1","b2","b3","b4","b5","b6","b7","b8")
ranames=names[porder]##(mean,cov,q7max,q7min,bankful,fr,T50,p,c,m,HFE,LFE,peak,SI,FD,zfe)
ranames[1]="Qmean";ranames[2]="CV";ranames[13]="Peaktime"
rg_names=data.frame(name=matrix(rep(sdnames,6)))
rg_names=data.frame(name=matrix(rep(sdnames,6)))

obs_data=data.frame(obs=matrix(t(obs1),96,1))
daymet_data=data.frame(daymet=matrix(t(daymet1),96,1))
prism_data=data.frame(prism=matrix(t(prism1),96,1))
maurer_data=data.frame(maurer=matrix(t(maurer1),96,1))
watsed=c("A1","A2","B1","B21","C1","C21")
watershed1=matrix((rep(watsed,each=16)))
a1=rep(0,16)
b1=c(100,8,450,3,550,150,250,1,1,0.3,15,6,300,0.5,15,1)
DATA=cbind(obs_data,daymet_data,prism_data,maurer_data,rg_names,watershed1)


###########comapre observed and Daymet ###################

lm_eqn = function(df){
  m = lm(obs~daymet, data=df)
  a=df$obs;b=df$daymet
  d=NSE(a,b)
  eq <- list(f = format(summary(m)$r.squared, digits = 2),t=format(d,digits=2))
  
  
  c(eq = paste("NSE=",as.character(as.expression(eq$t)),",R2=",as.character(as.expression(eq$f)),sep=""))               
}


DATA2 <- DATA
labeldata2 <- ddply(DATA2,.(name),lm_eqn )
##label position for daymet vs observed

x=c(20,2.5,150,3.5,200,55,170,0.5,0.35,0.12,5,1,100,-0.1,2,1)
y=c(70,6,375,14,570,135,255,0.6,0.6,0.3,14,4,310,0.67,15,0.9)
p=c(24,21,24,22,21, 22)
s=c(6,6,6,6,6,6)

labeldata=cbind(labeldata2,x,y)
levels(DATA2$name)[levels(DATA2$name)==sdnames] <- ranames
cbPalette <- c("blue", "blue", "gray60", "gray60", "gray60", "red")
ggplot(DATA2,aes(obs,daymet)) + geom_point(size=4,shape =p,colour=cbPalette,fill=cbPalette)+
  geom_smooth(method="lm",colour="black",linetype="dotted",size = 1.1) + 
  xlab("Observe") + ylab("Daymet") +
  #ggtitle("Stream flow Regime Variables") + 
  geom_text(data=labeldata, aes(x=x, y=y, label=eq,family="Times", fontface="bold"), parse=FALSE) +
  facet_wrap(~name,scales="free")+
  
  theme_bw()+theme(legend.position="right")+
  theme(strip.text.x = element_text(size=14, face="bold"),
        axis.text.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=16,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=16,angle=90,hjust=.5,vjust=.5,face="bold"))

ggsave("regime_estimate.tiff", width=12, height=13, dpi=500)
dev.off()




###########comapre observed and prism ###################

lm_eqn = function(df){
  m = lm(obs~prism, data=df)
  a=df$obs;b=df$prism
  d=NSE(a,b)
  eq <- list(f = format(summary(m)$r.squared, digits = 2),t=format(d,digits=2))
  
  
  c(eq = paste("NSE=",as.character(as.expression(eq$t)),",R2=",as.character(as.expression(eq$f)),sep=""))               
}

DATA2 <- DATA
labeldata2 <- ddply(DATA2,.(name),lm_eqn )
##label position for daymet vs observed

x=c(20,2.5,150,3.5,200,55,170,0.5,0.35,0.12,5,1,100,-0.1,2,1)
y=c(70,6,375,14,570,135,255,0.6,0.6,0.3,14,4,310,0.67,15,0.9)
p=c(24,21,24,22,21, 22)
s=c(6,6,6,6,6,6)

labeldata=cbind(labeldata2,x,y)
levels(DATA2$name)[levels(DATA2$name)==sdnames] <- ranames
cbPalette <- c("blue", "blue", "gray60", "gray60", "gray60", "red")
ggplot(DATA2,aes(obs,prism)) + geom_point(size=4,shape =p,colour=cbPalette,fill=cbPalette)+
  geom_smooth(method="lm",colour="black",linetype="dotted",size = 1.1) + 
  xlab("Observe") + ylab("Daymet") +
  #ggtitle("Stream flow Regime Variables") + 
  geom_text(data=labeldata, aes(x=x, y=y, label=eq,family="Times", fontface="bold"), parse=FALSE) +
  facet_wrap(~name,scales="free")+
  
  theme_bw()+theme(legend.position="right")+
  theme(strip.text.x = element_text(size=14, face="bold"),
        axis.text.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=16,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=16,angle=90,hjust=.5,vjust=.5,face="bold"))




###########comapre observed and maurer ###################


lm_eqn = function(df){
  m = lm(obs~maurer, data=df)
  a=df$obs;b=df$maurer
  d=NSE(a,b)
  eq <- list(f = format(summary(m)$r.squared, digits = 2),t=format(d,digits=2))
  
  
  c(eq = paste("NSE=",as.character(as.expression(eq$t)),",R2=",as.character(as.expression(eq$f)),sep=""))               
}


DATA2 <- DATA
labeldata2 <- ddply(DATA2,.(name),lm_eqn )
##label position for daymet vs observed

x=c(20,2.5,150,3.5,200,55,170,0.5,0.35,0.12,5,1,100,-0.1,2,1)
y=c(70,6,375,14,570,135,255,0.6,0.6,0.3,14,4,310,0.67,15,0.9)
p=c(24,21,24,22,21, 22)
s=c(6,6,6,6,6,6)
#check


labeldata=cbind(labeldata2,x,y)
levels(DATA2$name)[levels(DATA2$name)==sdnames] <- ranames
cbPalette <- c("blue", "blue", "gray60", "gray60", "gray60", "red")
ggplot(DATA2,aes(obs,maurer)) + geom_point(size=4,shape =p,colour=cbPalette,fill=cbPalette)+
  geom_smooth(method="lm",colour="black",linetype="dotted",size = 1.1) + 
  xlab("Observe") + ylab("Maurer") +
  #ggtitle("Stream flow Regime Variables") + 
  geom_text(data=labeldata, aes(x=x, y=y, label=eq,family="Times", fontface="bold"), parse=FALSE) +
  facet_wrap(~name,scales="free")+
  
  theme_bw()+theme(legend.position="right")+
  theme(strip.text.x = element_text(size=14, face="bold"),
        axis.text.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=16,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=16,angle=90,hjust=.5,vjust=.5,face="bold"))







