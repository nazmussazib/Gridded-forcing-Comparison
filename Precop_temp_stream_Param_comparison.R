

library(RNCEP)
library(zoo)
library(hydroTSM)
library(ggplot2)

###determine monthly rainfall  changes########

rainfall=function(file){
  rain=read.table(file,skip=3)
  dates_all=seq(as.Date("1980/1/1"), as.Date("2012/12/31"), "day")
  dates_bl=seq(as.Date("1980/1/1"), as.Date("2010/12/31"), "day")
  dates_match=match(dates_bl,dates_all)
  rainfall=rain[dates_match,(1:(ncol(rain)-2))]
  rain_mon=matrix(colMeans(data.frame(monthlyfunction(rainfall, FUN=sum, na.rm=TRUE,dates=dates_bl))/31))
  rain_ann=matrix(rowMeans(data.frame(annualfunction(rainfall, FUN=sum, na.rm=TRUE,dates=dates_bl))/31))
  rain_all=rbind(rain_ann,rain_mon)
}

rainfall_mon_ann=matrix(NA,nrow=13,ncol=1)
watsed=c("A1","A2","B1","B21","C1","C22")
rainfall_source=c("Daymet_RUN","PRISM_RUN","Maurer_RUN")
for ( i in 1:length(watsed)){
  for ( j in 1:length(rainfall_source)) {
   
  setwd(paste("E:\\USU_Research_work\\Other_Work\\RAINFALL_COMPARISON\\",watsed[i],"_watershed\\",rainfall_source[j],sep=""))
  rainfile=list.files(path =paste("E:\\USU_Research_work\\Other_Work\\RAINFALL_COMPARISON\\",watsed[i],"_watershed\\",rainfall_source[j],sep=""),pattern ="rain.dat")
  mr=sapply(rainfile,rainfall)
  rainfall_mon_ann=cbind(rainfall_mon_ann,mr)
 
}
print (i)
}


## differnece in annual means, stdv and annua coorealtuon


rainfall_stat_ann=function(file){
  rain=read.table(file,skip=3)
  dates_all=seq(as.Date("1980/1/1"), as.Date("2012/12/31"), "day")
  dates_bl=seq(as.Date("1980/1/1"), as.Date("2010/12/31"), "day")
  dates_match=match(dates_bl,dates_all)
  rainfall=matrix(rowMeans(rain[dates_match,(1:(ncol(rain)-2))]))
  rain_ann_mean=matrix((daily2annual(rainfall, FUN=sum, na.rm=TRUE,dates=dates_bl)))
  rain_ann_std=matrix((daily2annual(rainfall, FUN=sd, na.rm=TRUE,dates=dates_bl))) 
  rain_ann_max=matrix((daily2annual(rainfall, FUN=max, na.rm=TRUE,dates=dates_bl))) 
  rain_all=rbind(rain_ann_mean,rain_ann_std,rain_ann_max)
}




rainfall_ann_stat=matrix(NA,nrow=93,ncol=1)
watsed=c("A1","A2","B1","B21","C1","C22")
rainfall_source=c("Daymet_RUN","PRISM_RUN","Maurer_RUN")
for ( i in 1:length(watsed)){
  for ( j in 1:length(rainfall_source)) {
    
    setwd(paste("E:\\USU_Research_work\\Other_Work\\RAINFALL_COMPARISON\\",watsed[i],"_watershed\\",rainfall_source[j],sep=""))
    rainfile=list.files(path =paste("E:\\USU_Research_work\\Other_Work\\RAINFALL_COMPARISON\\",watsed[i],"_watershed\\",rainfall_source[j],sep=""),pattern ="rain.dat")
    mr=sapply(rainfile,rainfall_stat_ann)
    rainfall_ann_stat=cbind(rainfall_ann_stat,mr)
    
  }
  print (i)
}

rainfall_ann_stat=rainfall_ann_stat[,-1]
daymet_stat_prcp=rainfall_ann_stat[,seq(1,18,3)]
prism_stat_prcp=rainfall_ann_stat[,seq(2,18,3)]
Maurer_stat_prcp=rainfall_ann_stat[,seq(3,18,3)]
diff_D_M=(1-daymet_stat_prcp/Maurer_stat_prcp)*100
diff_P_M=(1-prism_stat_prcp/Maurer_stat_prcp)*100


diff_D_M1=matrix(diff_D_M,93*6,1)
diff_P_M1=matrix(diff_P_M,93*6,1)
diff_val=rbind(diff_D_M1,diff_P_M1)

diff_src=c("Daymet","PRISM")
dif_src=matrix((rep(diff_src,each=558)))

wats=c("A1","A2","B1","B21","C1","C21")
wats1=matrix((rep(watsed,each=93)))
wats2=matrix((rep(wats1,2)))

var=c("sum","std","max")
var1=matrix((rep(var,each=31)))
var2=matrix((rep(var1,12)))


x <- data.frame(
  values = diff_val,
  scenario=  dif_src,
  time = var2,
  place = wats2
)

# compare different sample populations across various temperatures
ggplot(x, aes(x = place, y = values, fill = scenario)) +
  geom_boxplot(outlier.shape = NA) +geom_hline(yintercept = 0.0,colour="black",size=1)+ #stat_summary(aes(group=scenario), fun.y=mean, geom="line")+
  facet_wrap(~ time,scales = "free")+theme_bw()+
  #theme(legend.position="none")+
  scale_fill_manual(name = "Emission Scenario", values = c("slateblue", "red"))+
  theme(strip.text.x = element_text(size=14, face="bold"),
        axis.text.x = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=16,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=16,angle=90,hjust=.5,vjust=.5,face="bold"))

##########monthly change #################


rainfall_stat_mon=function(file){
  rain=read.table(file,skip=3)
  dates_all=seq(as.Date("1980/1/1"), as.Date("2012/12/31"), "day")
  dates_bl=seq(as.Date("1995/1/1"), as.Date("2005/12/31"), "day")
  dates_match=match(dates_bl,dates_all)
  rainfall=rowMeans(rain[dates_match,(1:(ncol(rain)-2))])
  z <- vector2zoo(rainfall, as.Date(dates_bl))
  dsea=c("DJF","MAM","JJA","SON")
  
  rain_mon_mean=matrix(NA,nrow=11,ncol=4)
  
  for ( i in 1:4){
  rain_mon_mean[,i]=dm2seasonal(z,season=dsea[i],FUN=sum,na.rm=TRUE)
  }
  rain_mon_mean1=matrix(rain_mon_mean,11*4,1)
  
}

rainfall_mon_stat=matrix(NA,nrow=11*4,ncol=1)
watsed=c("A1","A2","B1","B21","C1","C22")
rainfall_source=c("Daymet_RUN","PRISM_RUN","Maurer_RUN")
for ( i in 1:length(watsed)){
  for ( j in 1:length(rainfall_source)) {
    
    setwd(paste("E:\\USU_Research_work\\Other_Work\\RAINFALL_COMPARISON\\",watsed[i],"_watershed\\",rainfall_source[j],sep=""))
    rainfile=list.files(path =paste("E:\\USU_Research_work\\Other_Work\\RAINFALL_COMPARISON\\",watsed[i],"_watershed\\",rainfall_source[j],sep=""),pattern ="rain.dat")
    mr=sapply(rainfile,rainfall_stat_mon)
    rainfall_mon_stat=cbind(rainfall_mon_stat,mr)
    
  }
  print (i)
}

rainfall_mon_stat=rainfall_mon_stat[,-1]


daymet_stat_prcp_mon=rainfall_mon_stat[,seq(1,18,3)]
prism_stat_prcp_mon=rainfall_mon_stat[,seq(2,18,3)]
Maurer_stat_prcp_mon=rainfall_mon_stat[,seq(3,18,3)]
diff_D_M_mon=(1-daymet_stat_prcp_mon/Maurer_stat_prcp_mon)*100
diff_P_M_mon=(1-prism_stat_prcp_mon/Maurer_stat_prcp_mon)*100


diff_D_M_mon1=matrix(diff_D_M_mon,11*4*6,1)
diff_P_M_mon1=matrix(diff_P_M_mon,11*4*6,1)
diff_val=rbind(diff_D_M_mon1,diff_P_M_mon1)

diff_src=c("Daymet","PRISM")
dif_src=matrix((rep(diff_src,each=11*4*6)))

wats=c("A1","A2","B1","B21","C1","C21")
wats1=matrix((rep(watsed,each=11*4*1)))
wats2=matrix((rep(wats1,2)))

##for temperature and precipitation
season1=matrix(rep(c('s1', 's2', 's3', 's4'),each=11))
season2=matrix((rep(season1,12)))
##for streamflow


diff_val[is.na(diff_val)]=0
diff_val[is.infinite(diff_val)]=0
diff_val[diff_val<(-100)]=-200




x <- data.frame(
  values = diff_val,
  scenario=  dif_src,
  time = season2,
  place = wats2
)

# compare different sample populations across various temperatures
ggplot(x, aes(x = time, y = values, fill = scenario)) +
  geom_boxplot(outlier.shape = NA) +geom_hline(yintercept = 0.0,colour="black",size=1)+ #stat_summary(aes(group=scenario), fun.y=mean, geom="line")+
  facet_wrap(~ place,scales = "free")+theme_bw()+
  #theme(legend.position="none")+
  scale_fill_manual(name = "Emission Scenario", values = c("slateblue", "red"))+
  theme(strip.text.x = element_text(size=14, face="bold"),
        axis.text.x = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=16,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=16,angle=90,hjust=.5,vjust=.5,face="bold"))



























###plotting the differences of rainfall (taking DAYMET AS REFERENCE)

mon_rainfall=matrix(rainfall_mon_ann[2:13,2:19],12*18,1)
rainfall_src=c("Daymet","PRISM","Maurer")
ran_src=matrix((rep(rainfall_src,each=12)))
rain_src=matrix((rep(ran_src,6)))
wats=c("A1","A2","B1","B21","C1","C21")
wats1=matrix((rep(watsed,each=36)))


##for precipitation
mon11=matrix(rep(c('J1','J2', 'J3', 'J4', 'J5','J6','J7', 'J8', 'J9', 'k1','k2','k3'),18))


data <- data.frame(
  values = mon_rainfall,
  scenario= rain_src,
  time = mon11,
  place =wats1
)

ggplot(data, aes(time,values,group=scenario,color=scenario)) + geom_line() + facet_wrap(~place,scales = "free") +
  geom_line(data = data, aes(time,values))+geom_point()+theme_bw()+
  theme(legend.position="none")+
  scale_fill_manual(name = "Scenario", values = c("white", "grey"))+
  theme(strip.text.x = element_text(size=14, face="bold"),
        axis.text.x = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="bold"))+
  


ggplot(x, aes(x = time, y = values, fill = scenario)) +
  geom_line() + #stat_summary(aes(group=scenario), fun.y=mean, geom="line")+
  facet_wrap(~ place,scales = "free")+
  xlab("Season") + ylab("Rainfall Change (%)")+ggtitle(" Year 2076-2095 compare to 1986-2005 ")

xlab("Season") + ylab(expression(paste(bold("Temperature Change(")^bold("0"),bold("C)"),sep="")))+ggtitle(" Year 2076-2095 compare to 1986-2005 ")

scale_x_discrete(breaks=c("Ann", "s1", "s2", "s3", "s4"), labels=c("Ann", 'Winter', 'Spring', 'Summer', 'Fall'))+



########determiing the tmax, tmin and tran differnces #######



temperature=function(file){
  temp=read.table(file,skip=3)
  dates_all=seq(as.Date("1980/1/1"), as.Date("2012/12/31"), "day")
  dates_bl=seq(as.Date("1980/1/1"), as.Date("2012/12/31"), "day")
  dates_match=match(dates_bl,dates_all)
  temperature=temp[dates_match,(1:(ncol(temp)-2))]
  tmax=(temperature[1:length(dates_bl),seq(1,dim(temperature)[2],3)])
  tmin=(temperature[1:length(dates_bl),seq(2,dim(temperature)[2],3)])
  tran=tmax-tmin
  #avgT=0.5*(tmax+tmin)
  temp_mon_max=matrix(colMeans(data.frame(monthlyfunction(tmax, FUN=mean, na.rm=TRUE,dates=dates_bl))))
  temp_ann_max=matrix(rowMeans(data.frame(annualfunction(tmax, FUN=mean, na.rm=TRUE,dates=dates_bl))))
  temp_all_max=rbind(temp_ann_max,temp_mon_max)
  
  temp_mon_min=matrix(colMeans(data.frame(monthlyfunction(tmin, FUN=mean, na.rm=TRUE,dates=dates_bl))))
  temp_ann_min=matrix(rowMeans(data.frame(annualfunction(tmin, FUN=mean, na.rm=TRUE,dates=dates_bl))))
  temp_all_min=rbind(temp_ann_min,temp_mon_min)
  
  temp_mon_ran=matrix(colMeans(data.frame(monthlyfunction(tran, FUN=mean, na.rm=TRUE,dates=dates_bl))))
  temp_ann_ran=matrix(rowMeans(data.frame(annualfunction(tran, FUN=mean, na.rm=TRUE,dates=dates_bl))))
  temp_all_ran=rbind(temp_ann_ran,temp_mon_ran)
  
  tempall=cbind(temp_all_max,temp_all_min,temp_all_ran)
  return(tempall)
  
  }



temperature_mon_ann_all=matrix(NA,nrow=39,ncol=1)

watsed=c("A1","A2","B1","B21","C1","C22")
temp_source=c("Daymet_RUN","PRISM_RUN","Maurer_RUN")


for ( i in 1:length(watsed)){
  for ( j in 1:length(temp_source)) {

    setwd(paste("E:\\USU_Research_work\\Other_Work\\RAINFALL_COMPARISON\\",watsed[i],"_watershed\\",rainfall_source[j],sep=""))
    tempfile=list.files(path =paste("E:\\USU_Research_work\\Other_Work\\RAINFALL_COMPARISON\\",watsed[i],"_watershed\\",rainfall_source[j],sep=""),pattern ="tmaxtmintdew.dat")
    temp_r=sapply(tempfile,temperature)
    temperature_mon_ann_all=cbind(temperature_mon_ann_all,temp_r)
    
  }
}



########plottign temperature comparuiosn #####








########determiing monthly streamflow comparison with the observed streamflow data #######












###parameter comparison for different watershed ########


watsed=c("A1","A2","B1","B21","C1","C22")
rainfall_source=c("Daymet_RUN","PRISM_RUN","Maurer_RUN")

setwd("F:\\USU_Research_work_update_March_30_2014\\RAINFALL_COMPARISON\\A1_watershed\\Daymet_RUN")
dd=read.table('topinp.dat',skip=7)
data=dd[,1]


param=matrix(NA,17)
for ( i in 1:length(watsed)){
  for ( j in 1:length(rainfall_source)) {
    
    setwd(paste("F:\\USU_Research_work_update_March_30_2014\\RAINFALL_COMPARISON\\",watsed[i],"_watershed\\",rainfall_source[j],sep=""))
    tempfile=list.files(path =paste("F:\\USU_Research_work_update_March_30_2014\\RAINFALL_COMPARISON\\",watsed[i],"_watershed\\",rainfall_source[j],sep=""),pattern ="topinp.dat")
    dd=read.table(tempfile,skip=7)
    data=dd[,1]
    param=cbind(param,data)
    
  }
}

param_need=matrix(param[c(1,2,3,4,5,9,10,16),2:19],8*18,1)
rainfall_src=c("Daymet","PRISM","Maurer")
ran_src=matrix((rep(rainfall_src,each=8)))
rain_src=matrix((rep(ran_src,6)))
wats=c("A1","A2","B1","B21","C1","C21")
wats1=matrix((rep(watsed,each=24)))


##for precipitation
parm1=matrix(rep(c('J1','J2', 'J3', 'J4', 'J5','J6','J7', 'J8'),18))


data <- data.frame(
  values = param_need,
  scenario= rain_src,
  time = parm1,
  place =wats1
)

ggplot(data, aes(place,values,group=scenario,color=scenario,shape=scenario)) + facet_wrap(~time,scales = "free") +
  geom_line(data = data, aes(place,values))+geom_point()+theme_bw()+
  theme(legend.position="none")+
  scale_fill_manual(name = "Scenario", values = c("white", "grey"))+
  theme(strip.text.x = element_text(size=14, face="bold"),
        axis.text.x = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="bold"))+
  
  
  
  ggplot(x, aes(x = time, y = values, fill = scenario)) +
  geom_line() + #stat_summary(aes(group=scenario), fun.y=mean, geom="line")+
  facet_wrap(~ place,scales = "free")+
  xlab("Season") + ylab("Rainfall Change (%)")+ggtitle(" Year 2076-2095 compare to 1986-2005 ")

xlab("Season") + ylab(expression(paste(bold("Temperature Change(")^bold("0"),bold("C)"),sep="")))+ggtitle(" Year 2076-2095 compare to 1986-2005 ")

scale_x_discrete(breaks=c("Ann", "s1", "s2", "s3", "s4"), labels=c("Ann", 'Winter', 'Spring', 'Summer', 'Fall'))+
  
####streamflwo monyhly  comparioson for different watersheds ###############

banum=c(29,113,25,135,131,43)
watsed=c("A1","A2","B1","B21","C1","C22")
stream_source=c("Daymet_RUN","PRISM_RUN","Maurer_RUN")
SNVAR_ALL=matrix(NA,nrow=12,ncol=1)
for ( i in 1:length(watsed)){
  for ( j in 1:length(stream_source)) {
   
    setwd(paste("F:\\USU_Research_work_update_March_30_2014\\RAINFALL_COMPARISON\\",watsed[i],"_watershed\\",stream_source[j],sep=""))
    time=seq(as.Date("2000/1/1"), as.Date("2010/12/31"), "day")
    time1=seq(as.Date("2000/1/1"), as.Date("2010/12/31"), "day")
    date_match=match(time1,time) # get overlap time interval

    sf=scan('streamflow_calibration.dat', what="")
    sf1=sf[seq(21,length(sf),1)] ## need to change theis things later
    obs_flow=(as.numeric(sf1[seq(1,length(sf1),3)])) ##need to change this later
    time_all_observ= seq(as.Date("1980/1/1"), as.Date("2012/12/31"), "day") #create time series
    date_overlap_obs=match(time1,time_all_observ) # get overlap time interval
    obsflow=data.frame(date=time1,q=obs_flow[date_overlap_obs])
    
    SN1=matrix(colMeans(data.frame(monthlyfunction(matrix(obsflow$q), FUN=mean, na.rm=TRUE,dates=time1))))
    
    ff=scan('FlowAtStreamNodes_cms.txt', what="")
    l=banum[i]+2#for A2
    ff1=ff[seq(l,length(ff),1)] ## need to change theis things later
    simu_flow=matrix(as.numeric(ff1[seq(2,length(ff1),l-1)])) ##need to change this later
    simflow=data.frame(date=time1,q=simu_flow[date_match])
    simflow[simflow<0.001]=0.001
    
    SN2=matrix(colMeans(data.frame(monthlyfunction(matrix(simflow$q), FUN=mean, na.rm=TRUE,dates=time1))))
    SN=cbind(SN1,SN2)
    SNVAR_ALL=cbind(SNVAR_ALL,SN)
    
    
  }
  
}

rm_col=c(1,4,6,10,12,16,18,20,24,28,30,34,36) ##remove duplicate observed data
SNVAR_ALL=SNVAR_ALL[ ,-rm_col]  
mon_srflow=matrix(SNVAR_ALL,12*24,1)
srflow_src=c("Observed","Daymet","PRISM","Maurer")
srf_src=matrix((rep(srflow_src,each=12)))
streamflow_src=matrix((rep(srf_src,6)))
wats=c("A1","A2","B1","B21","C1","C22")
wats1=matrix((rep(watsed,each=48)))


##for precipitation
mon11=matrix(rep(c('J1','J2', 'J3', 'J4', 'J5','J6','J7', 'J8', 'J9', 'k1','k2','k3'),24))


data <- data.frame(
  values = mon_srflow,
  scenario= streamflow_src,
  time = mon11,
  place =wats1
)


##calcualte NSE and write in the plot
eq1=matrix('str',1,1)
eq2=matrix('str',1,1)
eq3=matrix('str',1,1)
obs_ord=c(1,5,9,13,17,21)
for ( i in 1:6){
  obs=matrix(SNVAR_ALL[,obs_ord[i]])
  sim1=matrix(SNVAR_ALL[,obs_ord[i]+1])
  sim2=matrix(SNVAR_ALL[,obs_ord[i]+2])
  sim3=matrix(SNVAR_ALL[,obs_ord[i]+3])
  d1=NSE(obs,sim1)
  d2=NSE(obs,sim2)
  d3=NSE(obs,sim3)
  eq <- list(t1=format(d1,digits=2),t2=format(d2,digits=2),t3=format(d3,digits=2))
  txt1=c(paste("Daymet=",as.character(as.expression(eq$t1)),sep="")) 
  txt2=c(paste("PRISM=",as.character(as.expression(eq$t2)),sep="")) 
  txt3=c(paste("Maurer=",as.character(as.expression(eq$t3)),sep="")) 
  eq1=rbind(eq1,txt1)
  eq2=rbind(eq2,txt2)
  eq3=rbind(eq3,txt3)
  }
eq1=eq1[-1,]
eq1=data.frame(eq=eq1,place=wats)
eq2=eq2[-1,]
eq2=data.frame(eq=eq2,place=wats)
eq3=eq3[-1,]
eq3=data.frame(eq=eq3,place=wats)

x1=c(5,5,5,5,5,5)
y1=c(3,300,2,7,60,40)
y2=c(2,250,1.7,5,50,30)
y3=c(1,220,1.5,4,40,20)
labeldata1=data.frame(cbind(eq1,x,y1))
labeldata2=data.frame(cbind(eq2,x,y2))
labeldata3=data.frame(cbind(eq3,x,y3))

ggplot(data, aes(time,values)) + 
geom_line(data = data, aes(time,values,group=scenario,color=scenario))+geom_point()+
#geom_text(data=labeldata1, aes(x=x, y=y1, label=eq,family="Times", fontface="bold"), parse=FALSE)+ 
#geom_text(data=labeldata2, aes(x=x, y=y2, label=eq,family="Times", fontface="bold"), parse=FALSE)+  
#geom_text(data=labeldata3, aes(x=x, y=y3, label=eq,family="Times", fontface="bold"), parse=FALSE)+  
facet_wrap(~place,scales = "free") +
  
  
  


























