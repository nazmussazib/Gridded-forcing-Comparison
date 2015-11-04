

library(RNCEP)
library(zoo)
library(hydroTSM)
library(ggplot2)

###determine monthly rainfall  changes########

rainfall=function(file){
  rain=read.table(file,skip=3)
  dates_all=seq(as.Date("1980/1/1"), as.Date("2012/12/31"), "day")
  dates_bl=seq(as.Date("1980/1/1"), as.Date("2012/12/31"), "day")
  dates_match=match(dates_bl,dates_all)
  rainfall=rain[dates_match,(1:(ncol(rain)-2))]
  rain_mon=matrix(colMeans(data.frame(monthlyfunction(rainfall, FUN=sum, na.rm=TRUE,dates=dates_bl))/33))
  rain_ann=matrix(rowMeans(data.frame(annualfunction(rainfall, FUN=sum, na.rm=TRUE,dates=dates_bl))/33))
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

watsed=c("A1","A2","B1","B21","C1","C22")
stream_source=c("Daymet_RUN","PRISM_RUN","Maurer_RUN")
SNVAR_ALL=matrix(NA,nrow=12,ncol=1)
for ( i in 1:length(watsed)){
  for ( j in 1:length(stream_source)) {
   
    setwd(paste("F:\\USU_Research_work_update_March_30_2014\\RAINFALL_COMPARISON\\",watsed[i],"_watershed\\",stream_source[j],sep=""))
    time=seq(as.Date("2000/1/1"), as.Date("2010/12/31"), "day")
    time1=seq(as.Date("2001/1/1"), as.Date("2010/12/31"), "day")
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
wats=c("A1","A2","B1","B21","C1","C21")
wats1=matrix((rep(watsed,each=48)))


##for precipitation
mon11=matrix(rep(c('J1','J2', 'J3', 'J4', 'J5','J6','J7', 'J8', 'J9', 'k1','k2','k3'),24))


data <- data.frame(
  values = mon_srflow,
  scenario= streamflow_src,
  time = mon11,
  place =wats1
)

ggplot(data, aes(time,values,group=scenario,color=scenario)) + geom_line() + facet_wrap(~place,scales = "free") +
  geom_line(data = data, aes(time,values))+geom_point()+theme_bw()+
  #theme(legend.position="none")+
  scale_fill_manual(name = "Scenario", values = c("white", "grey"))+
  theme(strip.text.x = element_text(size=14, face="bold"),
        axis.text.x = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="plain"),
        axis.text.y = element_text(colour="grey20",size=12,angle=0,hjust=1,vjust=0,face="plain"),  
        axis.title.x = element_text(colour="grey20",size=12,angle=0,hjust=.5,vjust=0,face="bold"),
        axis.title.y = element_text(colour="grey20",size=12,angle=90,hjust=.5,vjust=.5,face="bold"))+
  
  
  


























