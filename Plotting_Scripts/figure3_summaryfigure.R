library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(ggthemes)


#Read in a summary table of fraction infected and defoliation under different initial insect densities (RCP 8.5, End of century)
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

table_base<-table
table_base$title<-"RCP 8.5, end of century"

# The values y>0 indicate the extent of the harmful effects of climate change to forests. 
# The decrease in fraction infected is on the same side as the increase in defoliation.
color<-c("Fraction infected mean"="#a2cd5a","Fraction infected 5th/95th %iles"="#00563f","Defoliation mean"="#cc397b","Defoliation 5th/95th %iles"="#8b357f","Density frequencies"="grey50")
linetype_list<-c("Fraction infected mean"=1,"Fraction infected 5th/95th %iles"=1,"Defoliation mean"=2,"Defoliation 5th/95th %iles"=2)




END<-ggplot(table_base) + theme_bw() + theme_clean()+
  geom_point(data=table_base[table_base$group=="Fraction Infected",],aes(x=insect_density,y=-(mean),color="Fraction infected mean"),size=3)+
  geom_line(data=table_base[table_base$group=="Fraction Infected",],aes(x=insect_density,y=-(mean),color="Fraction infected mean",linetype="Fraction infected mean"),linewidth=1)+
  geom_point(data=table_base[table_base$group=="Fraction Infected",],aes(x=insect_density,y=-(X5th_pct),color="Fraction infected 5th/95th %iles"),size=3)+
  geom_line(data=table_base[table_base$group=="Fraction Infected",],aes(x=insect_density,y=-(X5th_pct),color="Fraction infected 5th/95th %iles",linetype="Fraction infected 5th/95th %iles"),linewidth=1)+
  geom_point(data=table_base[table_base$group=="Fraction Infected",],aes(x=insect_density,y=-(X95th_pct),color="Fraction infected 5th/95th %iles"),size=3)+
  geom_line(data=table_base[table_base$group=="Fraction Infected",],aes(x=insect_density,y=-(X95th_pct),color="Fraction infected 5th/95th %iles",linetype="Fraction infected 5th/95th %iles"),linewidth=1)+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)),limits=c(6,100000))+
  scale_y_continuous(limits=c(-12,43),breaks=integer_breaks(),sec.axis = sec_axis(trans=~./1,name="Increase in defoliation (%)",breaks=c(-10,0,10,20,30,40)))+
  geom_point(data=table_base[table_base$group=="Defoliation",],aes(x=insect_density,y=(mean),color="Defoliation mean"),size=3)+
  geom_line(data=table_base[table_base$group=="Defoliation",],aes(x=insect_density,y=(mean),color="Defoliation mean",linetype="Defoliation mean"),linewidth=1)+
  geom_point(data=table_base[table_base$group=="Defoliation",],aes(x=insect_density,y=(X5th_pct),color="Defoliation 5th/95th %iles"),size=3)+
  geom_line(data=table_base[table_base$group=="Defoliation",],aes(x=insect_density,y=(X5th_pct),color="Defoliation 5th/95th %iles",linetype="Defoliation 5th/95th %iles"),linewidth=1)+
  geom_point(data=table_base[table_base$group=="Defoliation",],aes(x=insect_density,y=(X95th_pct),color="Defoliation 5th/95th %iles"),size=3)+
  geom_line(data=table_base[table_base$group=="Defoliation",],aes(x=insect_density,y=(X95th_pct),color="Defoliation 5th/95th %iles",linetype="Defoliation 5th/95th %iles"),linewidth=1)+
  xlab(expression(bold(paste("Initial insect density (egg masses/ha)")) ))+ylab(expression(bold("Relative decrease in infection (%)")))+
  theme(axis.line.y.right = element_line(color = "#cc397b"), 
        axis.ticks.y.right = element_line(color = "#cc397b"),
        axis.text.y.right = element_text(color = "#cc397b"), 
        axis.title.y.right = element_text(color = "#cc397b"))+
  theme(plot.title=element_text(size=28,face="bold"), axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"),plot.background = element_blank(),
        legend.position = c(0.2,0.795),legend.title=element_blank(),legend.text=element_text(size=24),legend.background=element_blank(),legend.key=element_blank(),
        strip.background =element_rect(fill="grey50"),strip.text = element_text(colour = 'white'), legend.key.width = unit(3,"cm"))+
  theme(axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 10)))+
  scale_color_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles","Density frequencies"),values=color,name="value type")+
  scale_fill_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles","Density frequencies"),values=color,name="value type")+ 
  scale_linetype_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles"),values=linetype_list,name="value type")+
  facet_grid(. ~ title)+theme(strip.text.x = element_text(size=28,face="bold"))+theme(panel.spacing.x = unit(2, "lines"))




#Mid-century
table<-read.csv("fullHiMIDMfRel_hdefof_0.05_summarytable.csv",header=TRUE)
table$group<-as.character(table$group)
table$group[table$group=="Fraction_Infected"]<-"Fraction Infected"
table$group<-factor(table$group,levels=c("Fraction Infected","Defoliation"))
table$insect_density<-table$insect_density*10000/250

color<-c("Fraction infected mean"="#a2cd5a","Fraction infected 5th/95th %iles"="#00563f","Defoliation mean"="#cc397b","Defoliation 5th/95th %iles"="#8b357f","Density frequencies"="grey50")
linetype_list<-c("Fraction infected mean"=1,"Fraction infected 5th/95th %iles"=1,"Defoliation mean"=2,"Defoliation 5th/95th %iles"=2)

table_alt<-table
table_alt$group<-factor(table_alt$group,levels=c("Fraction Infected","Defoliation"))
table_alt$title<-"RCP 8.5, middle of Century"



MID<-ggplot(table_alt) + theme_bw() + theme_clean()+
  geom_point(data=table_alt[table_alt$group=="Fraction Infected",],aes(x=insect_density,y=-(mean),color="Fraction infected mean"),size=3)+
  geom_line(data=table_alt[table_alt$group=="Fraction Infected",],aes(x=insect_density,y=-(mean),color="Fraction infected mean",linetype="Fraction infected mean"),linewidth=1)+
  geom_point(data=table_alt[table_alt$group=="Fraction Infected",],aes(x=insect_density,y=-(X5th_pct),color="Fraction infected 5th/95th %iles"),size=3)+
  geom_line(data=table_alt[table_alt$group=="Fraction Infected",],aes(x=insect_density,y=-(X5th_pct),color="Fraction infected 5th/95th %iles",linetype="Fraction infected 5th/95th %iles"),linewidth=1)+
  geom_point(data=table_alt[table_alt$group=="Fraction Infected",],aes(x=insect_density,y=-(X95th_pct),color="Fraction infected 5th/95th %iles"),size=3)+
  geom_line(data=table_alt[table_alt$group=="Fraction Infected",],aes(x=insect_density,y=-(X95th_pct),color="Fraction infected 5th/95th %iles",linetype="Fraction infected 5th/95th %iles"),linewidth=1)+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)),limits=c(6,100000))+
  scale_y_continuous(limits=c(-12,43),breaks=integer_breaks(),sec.axis = sec_axis(trans=~./1,name="Increase in defoliation (%)",breaks=c(-10,0,10,20,30,40)))+
  geom_point(data=table_alt[table_alt$group=="Defoliation",],aes(x=insect_density,y=(mean),color="Defoliation mean"),size=3)+
  geom_line(data=table_alt[table_alt$group=="Defoliation",],aes(x=insect_density,y=(mean),color="Defoliation mean",linetype="Defoliation mean"),linewidth=1)+
  geom_point(data=table_alt[table_alt$group=="Defoliation",],aes(x=insect_density,y=(X5th_pct),color="Defoliation 5th/95th %iles"),size=3)+
  geom_line(data=table_alt[table_alt$group=="Defoliation",],aes(x=insect_density,y=(X5th_pct),color="Defoliation 5th/95th %iles",linetype="Defoliation 5th/95th %iles"),linewidth=1)+
  geom_point(data=table_alt[table_alt$group=="Defoliation",],aes(x=insect_density,y=(X95th_pct),color="Defoliation 5th/95th %iles"),size=3)+
  geom_line(data=table_alt[table_alt$group=="Defoliation",],aes(x=insect_density,y=(X95th_pct),color="Defoliation 5th/95th %iles",linetype="Defoliation 5th/95th %iles"),linewidth=1)+
  xlab(expression(bold(paste("Initial insect density (egg masses/ha)")) ))+ylab(expression(bold("Relative decrease in infection (%)")))+
  theme(axis.line.y.right = element_line(color = "#cc397b"), 
        axis.ticks.y.right = element_line(color = "#cc397b"),
        axis.text.y.right = element_text(color = "#cc397b"), 
        axis.title.y.right = element_text(color = "#cc397b"))+
  theme(plot.title=element_text(size=28,face="bold"), axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"),plot.background = element_blank(),
        legend.position = "none",legend.title=element_blank(),legend.text=element_text(size=24),legend.background=element_blank(),legend.key=element_blank(),
        strip.background =element_rect(fill="grey50"),strip.text = element_text(colour = 'white'), legend.key.width = unit(3,"cm"))+
  theme(axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 10)))+
  scale_color_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles","Density frequencies"),values=color,name="value type")+
  scale_fill_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles","Density frequencies"),values=color,name="value type")+ 
  scale_linetype_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles"),values=linetype_list,name="value type")+
  facet_grid(. ~ title)+theme(strip.text.x = element_text(size=28,face="bold"))+theme(panel.spacing.x = unit(2, "lines"))
#scale_x_continuous(trans=log_trans(),breaks=round(exp(-5:-1),3),limits=c(0.004,0.38))+
#facet_wrap(~group)+theme(strip.text.x = element_text(size=24,face="bold"),strip.text.y = element_text(size=24,face="bold"))+theme(panel.spacing.x = unit(2, "lines"))



#RCP 4.5
table<-read.csv("fullLowENDMfRel_hdefof_0.05_summarytable.csv",header=TRUE)
table$group<-as.character(table$group)
table$group[table$group=="Fraction_Infected"]<-"Fraction Infected"
table$group<-factor(table$group,levels=c("Fraction Infected","Defoliation"))
table$insect_density<-table$insect_density*10000/250

color<-c("Fraction infected mean"="#a2cd5a","Fraction infected 5th/95th %iles"="#00563f","Defoliation mean"="#cc397b","Defoliation 5th/95th %iles"="#8b357f","Density frequencies"="grey50")
linetype_list<-c("Fraction infected mean"=1,"Fraction infected 5th/95th %iles"=1,"Defoliation mean"=2,"Defoliation 5th/95th %iles"=2)

table_alt<-table
table_alt$group<-factor(table_alt$group,levels=c("Fraction Infected","Defoliation"))
table_alt$title<-"RCP 4.5, end of century"



RCP4<-ggplot(table_alt) + theme_bw() + theme_clean()+
  geom_point(data=table_alt[table_alt$group=="Fraction Infected",],aes(x=insect_density,y=-(mean),color="Fraction infected mean"),size=3)+
  geom_line(data=table_alt[table_alt$group=="Fraction Infected",],aes(x=insect_density,y=-(mean),color="Fraction infected mean",linetype="Fraction infected mean"),linewidth=1)+
  geom_point(data=table_alt[table_alt$group=="Fraction Infected",],aes(x=insect_density,y=-(X5th_pct),color="Fraction infected 5th/95th %iles"),size=3)+
  geom_line(data=table_alt[table_alt$group=="Fraction Infected",],aes(x=insect_density,y=-(X5th_pct),color="Fraction infected 5th/95th %iles",linetype="Fraction infected 5th/95th %iles"),linewidth=1)+
  geom_point(data=table_alt[table_alt$group=="Fraction Infected",],aes(x=insect_density,y=-(X95th_pct),color="Fraction infected 5th/95th %iles"),size=3)+
  geom_line(data=table_alt[table_alt$group=="Fraction Infected",],aes(x=insect_density,y=-(X95th_pct),color="Fraction infected 5th/95th %iles",linetype="Fraction infected 5th/95th %iles"),linewidth=1)+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)),limits=c(6,100000))+
  scale_y_continuous(limits=c(-12,43),breaks=integer_breaks(),sec.axis = sec_axis(trans=~./1,name="Increase in defoliation (%)",breaks=c(-10,0,10,20,30,40)))+
  geom_point(data=table_alt[table_alt$group=="Defoliation",],aes(x=insect_density,y=(mean),color="Defoliation mean"),size=3)+
  geom_line(data=table_alt[table_alt$group=="Defoliation",],aes(x=insect_density,y=(mean),color="Defoliation mean",linetype="Defoliation mean"),linewidth=1)+
  geom_point(data=table_alt[table_alt$group=="Defoliation",],aes(x=insect_density,y=(X5th_pct),color="Defoliation 5th/95th %iles"),size=3)+
  geom_line(data=table_alt[table_alt$group=="Defoliation",],aes(x=insect_density,y=(X5th_pct),color="Defoliation 5th/95th %iles",linetype="Defoliation 5th/95th %iles"),linewidth=1)+
  geom_point(data=table_alt[table_alt$group=="Defoliation",],aes(x=insect_density,y=(X95th_pct),color="Defoliation 5th/95th %iles"),size=3)+
  geom_line(data=table_alt[table_alt$group=="Defoliation",],aes(x=insect_density,y=(X95th_pct),color="Defoliation 5th/95th %iles",linetype="Defoliation 5th/95th %iles"),linewidth=1)+
  xlab(expression(bold(paste("Initial insect density (egg masses/ha)")) ))+ylab(expression(bold("Relative decrease in infection (%)")))+
  theme(axis.line.y.right = element_line(color = "#cc397b"), 
        axis.ticks.y.right = element_line(color = "#cc397b"),
        axis.text.y.right = element_text(color = "#cc397b"), 
        axis.title.y.right = element_text(color = "#cc397b"))+
  theme(plot.title=element_text(size=28,face="bold"), axis.text=element_text(size=24),
        axis.title=element_text(size=24,face="bold"),plot.background = element_blank(),
        legend.position = "none",legend.title=element_blank(),legend.text=element_text(size=24),legend.background=element_blank(),legend.key=element_blank(),
        strip.background =element_rect(fill="grey50"),strip.text = element_text(colour = 'white'), legend.key.width = unit(3,"cm"))+
  theme(axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 10)))+
  scale_color_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles","Density frequencies"),values=color,name="value type")+
  scale_fill_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles","Density frequencies"),values=color,name="value type")+ 
  scale_linetype_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles"),values=linetype_list,name="value type")+
  facet_grid(. ~ title)+theme(strip.text.x = element_text(size=28,face="bold"))+theme(panel.spacing.x = unit(2, "lines"))
#scale_x_continuous(trans=log_trans(),breaks=round(exp(-5:-1),3),limits=c(0.004,0.38))+
#facet_wrap(~group)+theme(strip.text.x = element_text(size=24,face="bold"),strip.text.y = element_text(size=24,face="bold"))+theme(panel.spacing.x = unit(2, "lines"))



png(paste("fullRel_hdefof_summarytable_to_plot_3in1.png",sep=""),height=1800,width=1600)
grid.arrange(END,MID,RCP4,nrow=3)
#grid.arrange(g_IT)

dev.off()  
