library(tidyverse)
library(ggthemes)
library(ggnewscale)

#Read in Michigan defoliation data
defo<-read.csv("MichiganDefoliationData.csv",header=TRUE)
defo<-defo[defo$Year<2023,]
defo$Hectares<-defo$Acres/2.471
g_defo<-ggplot(defo,aes(x=Year,y=Hectares/10^6))+theme_clean()+geom_point(size=3,color="darkmagenta")+geom_line(color="darkmagenta")+ylab("Defoliation (million hectares)")+
  scale_y_continuous(breaks=seq(0,1.5,by=0.25))+scale_x_continuous(limits=c(1980,2024))+
  theme(plot.title=element_text(size=28,face="bold"), axis.text=element_text(size=24),
        axis.title=element_text(size=28,face="bold"),plot.background = element_blank(),legend.key.size = unit(1.5, "cm"),
        legend.title=element_text(size=24,face="bold"),legend.text=element_text(size=24))

g_defo_sub<-ggplot(defo[defo$Year>=2005,],aes(x=Year,y=Hectares/(10^6)))+theme_clean()+geom_point(size=3,color="darkmagenta")+geom_line(color="darkmagenta")+
  scale_x_continuous(limits=c(2005,2026),breaks=seq(2006,2024,by=2))+ylab("Defoliation (million hectares)")+
  scale_y_continuous(breaks=seq(0,1.25,by=0.25))+
  theme(plot.title=element_text(size=28,face="bold"), axis.text=element_text(size=24),
        axis.title=element_text(size=28,face="bold"),plot.background = element_blank(),legend.key.size = unit(1.5, "cm"),
        legend.title=element_text(size=24,face="bold"),legend.text=element_text(size=24))

#Read in data collected from 2010-2012 (also provided in the folder Infection_Data, see file description there)
EXPdataO=read.csv("expdataO3.csv", header = TRUE, sep = ",")
EXPdataC=read.csv("expdataC3.csv", header = TRUE, sep = ",")
FERALdata=read.csv("feraldataALL2.csv", header = TRUE, sep = ",")

FERALdata<-FERALdata[!FERALdata$Pop %in% c(6,7),]
FERALdata$Fdata=FERALdata$fungus/(FERALdata$total+FERALdata$fungus+FERALdata$virus)

EXPdataC<-EXPdataC[!EXPdataC$Pop %in% c(6,7),]
EXPdataO<-EXPdataO[!EXPdataO$Pop %in% c(6,7),]

FERALdata$year<-2010
FERALdata$year[FERALdata$Pop==4]<-2011
FERALdata$year[FERALdata$Pop==5]<-2012
FERALdata$year[FERALdata$Pop==8]<-2011

EXPdataO$year<-2010
EXPdataO$year[EXPdataO$Pop==4]<-2011
EXPdataO$year[EXPdataO$Pop==5]<-2012
EXPdataO$year[EXPdataO$Pop==8]<-2011

EXPdataC$year<-2010
EXPdataC$year[EXPdataC$Pop==4]<-2011
EXPdataC$year[EXPdataC$Pop==5]<-2012
EXPdataC$year[EXPdataC$Pop==8]<-2011

#Calculate annual cumulative fraction infected (feral data, 2010-2012)
finf_season<-c()
vinf_season<-c()
year_array<-c()
for (i in 1:length(unique(FERALdata$Pop))){
  frame_sub<-FERALdata[FERALdata$Pop==unique(FERALdata$Pop)[i],]
  year_array[i]<-unique(FERALdata$year[FERALdata$Pop==unique(FERALdata$Pop)[i]])
  #if (max(frame_sub$Fdata)>0.1){
  finf_season[i]<-1-prod(1-frame_sub$Fdata)
  #}else{
  #  finf_season[i]<-NA
  #}
  #if (max(frame_sub$Vdata)>0.1){
  vinf_season[i]<-1-prod(1-frame_sub$Vdata)
  #}else{
  #  vinf_season[i]<-NA
  #}
}

feral_table<-data.frame(year_array,1:length(unique(FERALdata$Pop)),finf_season,vinf_season)
colnames(feral_table)[1:2]<-c("year","site")

#Read in data collected in 2021 (also provided in the folder Infection_Data, see file description there) and caulcated annual fraction infected
data21<-read.csv("MI_Observational_Sites_2021.csv",header=TRUE)
data21<-data21[data21$Site!="HPM",]

finf_season<-c()
vinf_season<-c()
year_array<-c()
for (i in 1:length(unique(data21$Site))){
  frame_sub<-data21[data21$Site==unique(data21$Site)[i],]
  #if (max(frame_sub$Frac.Fungus)>0.1){
  finf_season[i]<-1-prod(1-frame_sub$Frac.Fungus)
  #}else{
  #  finf_season[i]<-NA
  #}
  #if (max(frame_sub$Frac.Virus)>0.1){
  vinf_season[i]<-1-prod(1-frame_sub$Frac.Virus)
  #}else{
  #  vinf_season[i]<-NA
  #}
}
data21_table<-data.frame(2021,unique(data21$Site),finf_season,vinf_season)
colnames(data21_table)[1:2]<-c("year","site")

#Read in data collected in 2023 (also provided in the folder Infection_Data, see file description there) and caulcated annual fraction infected
data23<-read.csv("MI_Observational_Sites_2023.csv",header=TRUE)
finf_season<-c()
vinf_season<-c()
year_array<-c()
for (i in 1:length(unique(data23$Location))){
  frame_sub<-data23[data23$Location==unique(data23$Location)[i],]
  #if (max(frame_sub$Frac_F_Inf)>0.1){
  finf_season[i]<-1-prod(1-frame_sub$Frac_F_Inf)
  #}else{
  #  finf_season[i]<-NA
  #}
  #if (max(frame_sub$Frac_V_Inf)>0.1){
  vinf_season[i]<-1-prod(1-frame_sub$Frac_V_Inf)
  #}else{
  #  vinf_season[i]<-NA
  #}
}
data23_table<-data.frame(2023,unique(data23$Location),finf_season,vinf_season)
colnames(data23_table)[1:2]<-c("year","site")

#Calculate annual cumulative fraction infected (experimental data, 2010-2012)
cage_C<-as.data.frame(EXPdataC %>% group_by(Pop,year,Week) %>% summarise(week_fracf=mean(Idata)))
cage_O<-as.data.frame(EXPdataO %>% group_by(Pop,year,Week) %>% summarise(week_fracf=mean(Idata)))


finf_season<-c()
year_array<-c()
for (i in 1:length(unique(cage_C$Pop))){
  frame_sub<-cage_C[cage_C$Pop==unique(cage_C$Pop)[i],]
  year_array[i]<-unique(cage_C$year[cage_C$Pop==unique(cage_C$Pop)[i]])
  #if (max(frame_sub$week_fracf)>0.1){
  finf_season[i]<-1-prod(1-frame_sub$week_fracf)
  #}else{
  #  finf_season[i]<-NA
  #}
}

covercage_table<-data.frame(year_array,1:length(unique(cage_C$Pop)),finf_season)
colnames(covercage_table)<-c("year","site","finf_cover")

finf_season<-c()
year_array<-c()
for (i in 1:length(unique(cage_O$Pop))){
  frame_sub<-cage_O[cage_O$Pop==unique(cage_C$Pop)[i],]
  year_array[i]<-unique(cage_O$year[cage_O$Pop==unique(cage_C$Pop)[i]])
  #if (max(frame_sub$week_fracf)>0.1){
  finf_season[i]<-1-prod(1-frame_sub$week_fracf)
  #}else{
  #  finf_season[i]<-NA
  #}
}

opencage_table<-data.frame(year_array,1:length(unique(cage_O$Pop)),finf_season)
colnames(opencage_table)<-c("year","site","finf_open")

#Make a summary table of cumulative fractions infected from 2010-2012 data
Colin_table<-left_join(feral_table,covercage_table,by=c("year","site"))
Colin_table<-left_join(Colin_table,opencage_table,by=c("year","site"))
Colin_table$finf_season_mean<-(Colin_table$finf_season+Colin_table$finf_cover+Colin_table$finf_open)/3

Colin_table_sub<-Colin_table[,c("year","site","finf_season_mean","vinf_season")]
colnames(Colin_table_sub)[3]<-"finf_season"
Colin_table_sub$group<-"2010-2012"

#Combine 2010-2012 and 2021/2023 tables
data21_table$group<-"2021-2023"
data23_table$group<-"2021-2023"
alldata<-rbind(Colin_table_sub,data21_table,data23_table)


allsumm<-as.data.frame(alldata %>% group_by(group) %>% summarise(ave_f=mean(finf_season),ave_v=mean(vinf_season),se_f=sd(finf_season)/sqrt(n()),se_v=sd(vinf_season)/sqrt(n())))

#Rearrange the table for plotting
library(reshape2)
melt_mean<-melt(allsumm[,c("group","ave_f","ave_v")], id.vars = "group", variable.name = "quantity")
melt_se<-melt(allsumm[,c("group","se_f","se_v")], id.vars = "group", variable.name = "quantity")

colnames(melt_mean)[3]<-"mean"
colnames(melt_se)[3]<-"se"
melt_mean$quantity<-factor(melt_mean$quantity,levels=c("ave_f","ave_v"),labels=c("F","V"))
melt_se$quantity<-factor(melt_se$quantity,levels=c("se_f","se_v"),labels=c("F","V"))

long_summ<-left_join(melt_mean,melt_se,by=c("group","quantity"))


g_frac<-ggplot(long_summ)+theme_clean()+
  geom_point(aes(x=group,y=mean,color=quantity),size=5)+
  geom_line(aes(x=group,y=mean,group=quantity,color=quantity))+
  geom_errorbar(aes(x=group,y=mean,xmin=group,xmax=group,ymin=mean-se,ymax=mean+se,color=quantity,group=quantity),width=0.2)+
  ylab("Seasonal fraction infected")+xlab("Year")+
  scale_color_manual(values=c("F"="black","V=grey"))+
  scale_y_continuous(limits=c(0,0.75),breaks=seq(0,0.75,by=0.25))+
  theme(plot.title=element_text(size=28,face="bold"), axis.text=element_text(size=24),
        axis.title=element_text(size=28,face="bold"),plot.background = element_blank(),legend.key.size = unit(1.5, "cm"),
        legend.title=element_text(size=24,face="bold"),legend.text=element_text(size=24),legend.position="none")+
  geom_text(x=0.8,y=0.61,label="E. maimaiga",color="black",size=8)+
  geom_text(x=0.8,y=0.22,label="Virus",color="grey",size=8)+
  geom_point(data=alldata[alldata$site==1,],aes(x=group,y=finf_season),color="red",size=5)+
  geom_point(data=alldata[alldata$site==1,],aes(x=group,y=vinf_season),color="pink",size=5)+
  geom_text(x=0.8,y=0.745,label="South Only",color="red",size=8)+
  geom_text(x=0.8,y=0.03,label="South Only",color="pink",size=8)




library(gridExtra)

jpeg("EX_Figure9_data_compare_Michigan.jpg",height=6600,width=4500,res=300)
grid.arrange(g_defo,g_defo_sub,g_frac,nrow=3)
dev.off()


