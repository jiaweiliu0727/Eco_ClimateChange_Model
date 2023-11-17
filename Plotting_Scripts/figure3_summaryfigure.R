library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(ggthemes)

lit<-read.csv("cumfracinf_literature_final.csv") #Literature review of insect density and the corresponding fraction infected
fracinf_thres<-0.566667691 #Fraction infected threshold value to indicate whether there is an epizootic
lit<-lit[lit$cumulative_fracinf>fracinf_thres,]
lit$eggs_per_hectare[lit$eggs_per_hectare==0]<-7 #Some literature data has 0 egg masses/ha; turn them into very small values to take log
lit<-do.call("rbind", replicate(3, lit, simplify = FALSE)) #Replicate 3 times to adjust the height of the two y axis to similar ranges

#Read in a summary table of fraction infected and defoliation under different initial insect densities
table<-read.csv("fullHiENDMfRel_hdefof_0.05_summarytable.csv",header=TRUE)
table$group<-as.character(table$group)
table$group[table$group=="Fraction_Infected"]<-"Fraction Infected"
table$group<-factor(table$group,levels=c("Fraction Infected","Defoliation"))
table$insect_density<-table$insect_density*10000/250 #Larvae/m^2 into egg masses/ha

color<-c("Fraction infected mean"="#a2cd5a","Fraction infected 5th/95th %iles"="#00563f","Defoliation mean"="#cc397b","Defoliation 5th/95th %iles"="#8b357f","Density frequencies"="grey50")
linetype_list<-c("Fraction infected mean"=1,"Fraction infected 5th/95th %iles"=1,"Defoliation mean"=2,"Defoliation 5th/95th %iles"=2)


integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

table_norm<-table

# The values y>0 indicate the extent of the harmful effects of climate change to forests. 
# The decrease in fraction infected is on the same side as the increase in defoliation.
pd<-ggplot(table) + theme_bw() + theme_clean()+
  geom_histogram(data=lit,aes(x=eggs_per_hectare,fill="Density frequencies"),bins=28,alpha=0.2,color="black")+
  geom_point(data=table_norm[table_norm$group=="Fraction Infected",],aes(x=insect_density,y=-(mean),color="Fraction infected mean"),size=3)+
  geom_line(data=table_norm[table_norm$group=="Fraction Infected",],aes(x=insect_density,y=-(mean),color="Fraction infected mean",linetype="Fraction infected mean"),linewidth=1)+
  geom_point(data=table_norm[table_norm$group=="Fraction Infected",],aes(x=insect_density,y=-(X5th_pct),color="Fraction infected 5th/95th %iles"),size=3)+
  geom_line(data=table_norm[table_norm$group=="Fraction Infected",],aes(x=insect_density,y=-(X5th_pct),color="Fraction infected 5th/95th %iles",linetype="Fraction infected 5th/95th %iles"),linewidth=1)+
  geom_point(data=table_norm[table_norm$group=="Fraction Infected",],aes(x=insect_density,y=-(X95th_pct),color="Fraction infected 5th/95th %iles"),size=3)+
  geom_line(data=table_norm[table_norm$group=="Fraction Infected",],aes(x=insect_density,y=-(X95th_pct),color="Fraction infected 5th/95th %iles",linetype="Fraction infected 5th/95th %iles"),linewidth=1)+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)),limits=c(5,100000))+
  scale_y_continuous(limits=c(-2,46),breaks=integer_breaks(),sec.axis = sec_axis(trans=~./3,name="Number of cases",breaks=c(0,3,6,9,12,15)))+
  geom_point(data=table_norm[table_norm$group=="Defoliation",],aes(x=insect_density,y=(mean),color="Defoliation mean"),size=3)+
  geom_line(data=table_norm[table_norm$group=="Defoliation",],aes(x=insect_density,y=(mean),color="Defoliation mean",linetype="Defoliation mean"),linewidth=1)+
  geom_point(data=table_norm[table_norm$group=="Defoliation",],aes(x=insect_density,y=(X5th_pct),color="Defoliation 5th/95th %iles"),size=3)+
  geom_line(data=table_norm[table_norm$group=="Defoliation",],aes(x=insect_density,y=(X5th_pct),color="Defoliation 5th/95th %iles",linetype="Defoliation 5th/95th %iles"),linewidth=1)+
  geom_point(data=table_norm[table_norm$group=="Defoliation",],aes(x=insect_density,y=(X95th_pct),color="Defoliation 5th/95th %iles"),size=3)+
  geom_line(data=table_norm[table_norm$group=="Defoliation",],aes(x=insect_density,y=(X95th_pct),color="Defoliation 5th/95th %iles",linetype="Defoliation 5th/95th %iles"),linewidth=1)+
  xlab(expression(bold(paste("Initial insect density (egg masses/ha)")) ))+ylab(expression(bold("Change in infection/defoliation (%)")))+
  theme(axis.line.y.right = element_line(color = "grey30"), 
        axis.ticks.y.right = element_line(color = "grey30"),
        axis.text.y.right = element_text(color = "grey30"), 
        axis.title.y.right = element_text(color = "grey30"))+
  theme(plot.title=element_text(size=28,face="bold"), axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"),plot.background = element_blank(),
        legend.position = c(0.2,0.7),legend.title=element_blank(),legend.text=element_text(size=24),legend.background=element_blank(),
        strip.background =element_rect(fill="grey50"),strip.text = element_text(colour = 'white'), legend.key.width = unit(3,"cm"))+
  scale_color_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles","Density frequencies"),values=color,name="value type")+
  scale_fill_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles","Density frequencies"),values=color,name="value type")+ 
  scale_linetype_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles"),values=linetype_list,name="value type")


png(paste("figure3_summarytable.png",sep=""),height=600,width=1600)
grid.arrange(pd,nrow=1)

dev.off()   
