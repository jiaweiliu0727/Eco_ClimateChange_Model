library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)
library(ggthemes)

integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

#Read in the summary statistics of the basic model and the model considering heat-driven larval mortality
table<-read.csv("fullHiENDMfRel_hdefof_0.05_summarytable_ndm.csv",header=TRUE)
table$group<-as.character(table$group)
table$group[table$group=="Fraction_Infected"]<-"Fraction Infected"
table$model<-factor(table$model,levels=c("normal","ndm","ndmh"))
table$group<-factor(table$group,levels=c("Fraction Infected","Defoliation"))
table$insect_density<-table$insect_density*10000/250

color<-c("Fraction infected mean"="#a2cd5a","Fraction infected 5th/95th %iles"="#00563f","Defoliation mean"="#cc397b","Defoliation 5th/95th %iles"="#8b357f","Density frequencies"="grey50")
linetype_list<-c("Fraction infected mean"=1,"Fraction infected 5th/95th %iles"=1,"Defoliation mean"=2,"Defoliation 5th/95th %iles"=2)

table_base<-table[table$model=="normal",]
table_base$group<-factor(table_base$group,levels=c("Fraction Infected","Defoliation"))
table_base$title<-"No Temperature-Driven Mortality"



base<-ggplot(table_base) + theme_bw() + theme_clean()+
  geom_point(data=table_base[table_base$group=="Fraction Infected",],aes(x=insect_density,y=-(mean),color="Fraction infected mean"),size=5)+
  geom_line(data=table_base[table_base$group=="Fraction Infected",],aes(x=insect_density,y=-(mean),color="Fraction infected mean",linetype="Fraction infected mean"),linewidth=2)+
  geom_point(data=table_base[table_base$group=="Fraction Infected",],aes(x=insect_density,y=-(X5th_pct),color="Fraction infected 5th/95th %iles"),size=5)+
  geom_line(data=table_base[table_base$group=="Fraction Infected",],aes(x=insect_density,y=-(X5th_pct),color="Fraction infected 5th/95th %iles",linetype="Fraction infected 5th/95th %iles"),linewidth=2)+
  geom_point(data=table_base[table_base$group=="Fraction Infected",],aes(x=insect_density,y=-(X95th_pct),color="Fraction infected 5th/95th %iles"),size=5)+
  geom_line(data=table_base[table_base$group=="Fraction Infected",],aes(x=insect_density,y=-(X95th_pct),color="Fraction infected 5th/95th %iles",linetype="Fraction infected 5th/95th %iles"),linewidth=2)+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)),limits=c(6,100000))+
  scale_y_continuous(limits=c(-2,43),breaks=integer_breaks(),sec.axis = sec_axis(trans=~./1,name="Increase in defoliation (%)",breaks=c(-10,0,10,20,30,40)))+
  geom_point(data=table_base[table_base$group=="Defoliation",],aes(x=insect_density,y=(mean),color="Defoliation mean"),size=5)+
  geom_line(data=table_base[table_base$group=="Defoliation",],aes(x=insect_density,y=(mean),color="Defoliation mean",linetype="Defoliation mean"),linewidth=2)+
  geom_point(data=table_base[table_base$group=="Defoliation",],aes(x=insect_density,y=(X5th_pct),color="Defoliation 5th/95th %iles"),size=5)+
  geom_line(data=table_base[table_base$group=="Defoliation",],aes(x=insect_density,y=(X5th_pct),color="Defoliation 5th/95th %iles",linetype="Defoliation 5th/95th %iles"),linewidth=2)+
  geom_point(data=table_base[table_base$group=="Defoliation",],aes(x=insect_density,y=(X95th_pct),color="Defoliation 5th/95th %iles"),size=5)+
  geom_line(data=table_base[table_base$group=="Defoliation",],aes(x=insect_density,y=(X95th_pct),color="Defoliation 5th/95th %iles",linetype="Defoliation 5th/95th %iles"),linewidth=2)+
  xlab(expression(bold(paste("Initial insect density (egg masses/ha)")) ))+ylab(expression(bold("Relative decrease\n  in infection (%)")))+
  theme(axis.line.y.right = element_line(color = "#cc397b"), 
        axis.ticks.y.right = element_line(color = "#cc397b"),
        axis.text.y.right = element_text(color = "#cc397b"), 
        axis.title.y.right = element_text(color = "#cc397b"))+
  theme(plot.title=element_text(size=32,face="bold"), axis.text=element_text(size=30),
        axis.title=element_text(size=30,face="bold"),plot.background = element_blank(),
        legend.position = c(0.23,0.85),legend.title=element_blank(),legend.text=element_text(size=26),legend.background=element_blank(),legend.key=element_blank(),
        strip.background =element_rect(fill="grey50"),strip.text = element_text(colour = 'white'), legend.key.width = unit(3,"cm"))+
  theme(axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 10)))+
  scale_color_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles","Density frequencies"),values=color,name="value type")+
  scale_fill_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles","Density frequencies"),values=color,name="value type")+ 
  scale_linetype_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles"),values=linetype_list,name="value type")+
  facet_grid(. ~ title)+theme(strip.text.x = element_text(size=36,face="bold"))+theme(panel.spacing.x = unit(2, "lines"))+theme(plot.margin = unit(c(1,1,1,3), "cm"))
#scale_x_continuous(trans=log_trans(),breaks=round(exp(-5:-1),3),limits=c(0.004,0.38))+
#facet_wrap(~group)+theme(strip.text.x = element_text(size=30,face="bold"),strip.text.y = element_text(size=30,face="bold"))+theme(panel.spacing.x = unit(2, "lines"))



table_ndm<-table[table$model=="ndm",]
table_ndm$group<-factor(table_ndm$group,levels=c("Fraction Infected","Defoliation"))
table_ndm$title<-"Intermediate Temperature-Driven Mortality"

ndm<-ggplot(table_ndm) + theme_bw() + theme_clean()+
  geom_point(data=table_ndm[table_ndm$group=="Fraction Infected",],aes(x=insect_density,y=-(mean),color="Fraction infected mean"),size=5)+
  geom_line(data=table_ndm[table_ndm$group=="Fraction Infected",],aes(x=insect_density,y=-(mean),color="Fraction infected mean",linetype="Fraction infected mean"),linewidth=2)+
  geom_point(data=table_ndm[table_ndm$group=="Fraction Infected",],aes(x=insect_density,y=-(X5th_pct),color="Fraction infected 5th/95th %iles"),size=5)+
  geom_line(data=table_ndm[table_ndm$group=="Fraction Infected",],aes(x=insect_density,y=-(X5th_pct),color="Fraction infected 5th/95th %iles",linetype="Fraction infected 5th/95th %iles"),linewidth=2)+
  geom_point(data=table_ndm[table_ndm$group=="Fraction Infected",],aes(x=insect_density,y=-(X95th_pct),color="Fraction infected 5th/95th %iles"),size=5)+
  geom_line(data=table_ndm[table_ndm$group=="Fraction Infected",],aes(x=insect_density,y=-(X95th_pct),color="Fraction infected 5th/95th %iles",linetype="Fraction infected 5th/95th %iles"),linewidth=2)+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)),limits=c(6,100000))+
  scale_y_continuous(limits=c(-2,43),breaks=integer_breaks(),sec.axis = sec_axis(trans=~./1,name="Increase in defoliation (%)",breaks=c(-10,0,10,20,30,40)))+
  geom_point(data=table_ndm[table_ndm$group=="Defoliation",],aes(x=insect_density,y=(mean),color="Defoliation mean"),size=5)+
  geom_line(data=table_ndm[table_ndm$group=="Defoliation",],aes(x=insect_density,y=(mean),color="Defoliation mean",linetype="Defoliation mean"),linewidth=2)+
  geom_point(data=table_ndm[table_ndm$group=="Defoliation",],aes(x=insect_density,y=(X5th_pct),color="Defoliation 5th/95th %iles"),size=5)+
  geom_line(data=table_ndm[table_ndm$group=="Defoliation",],aes(x=insect_density,y=(X5th_pct),color="Defoliation 5th/95th %iles",linetype="Defoliation 5th/95th %iles"),linewidth=2)+
  geom_point(data=table_ndm[table_ndm$group=="Defoliation",],aes(x=insect_density,y=(X95th_pct),color="Defoliation 5th/95th %iles"),size=5)+
  geom_line(data=table_ndm[table_ndm$group=="Defoliation",],aes(x=insect_density,y=(X95th_pct),color="Defoliation 5th/95th %iles",linetype="Defoliation 5th/95th %iles"),linewidth=2)+
  xlab(expression(bold(paste("Initial insect density (egg masses/ha)")) ))+ylab(expression(bold("Relative decrease\n  in infection (%)")))+
  theme(axis.line.y.right = element_line(color = "#cc397b"), 
        axis.ticks.y.right = element_line(color = "#cc397b"),
        axis.text.y.right = element_text(color = "#cc397b"), 
        axis.title.y.right = element_text(color = "#cc397b"))+
  theme(plot.title=element_text(size=32,face="bold"), axis.text=element_text(size=30),
        axis.title=element_text(size=30,face="bold"),plot.background = element_blank(),
        legend.position = "none",legend.title=element_blank(),legend.text=element_text(size=30),legend.background=element_blank(),legend.key=element_blank(),
        strip.background =element_rect(fill="grey50"),strip.text = element_text(colour = 'white'), legend.key.width = unit(3,"cm"))+
  theme(axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 10)))+
  scale_color_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles","Density frequencies"),values=color,name="value type")+
  scale_fill_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles","Density frequencies"),values=color,name="value type")+ 
  scale_linetype_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles"),values=linetype_list,name="value type")+
  facet_grid(. ~ title)+theme(strip.text.x = element_text(size=36,face="bold"))+theme(panel.spacing.x = unit(2, "lines"))+theme(plot.margin = unit(c(1,1,1,3), "cm"))
#scale_x_continuous(trans=log_trans(),breaks=round(exp(-5:-1),3),limits=c(0.004,0.38))+
#facet_wrap(~group)+theme(strip.text.x = element_text(size=30,face="bold"),strip.text.y = element_text(size=30,face="bold"))+theme(panel.spacing.x = unit(2, "lines"))

table_ndmh<-table[table$model=="ndmh",]
table_ndmh$group<-factor(table_ndm$group,levels=c("Fraction Infected","Defoliation"))
table_ndmh$title<-"High Temperature-Driven Mortality"

ndmh<-ggplot(table_ndmh) + theme_bw() + theme_clean()+
  geom_point(data=table_ndmh[table_ndmh$group=="Fraction Infected",],aes(x=insect_density,y=-(mean),color="Fraction infected mean"),size=5)+
  geom_line(data=table_ndmh[table_ndmh$group=="Fraction Infected",],aes(x=insect_density,y=-(mean),color="Fraction infected mean",linetype="Fraction infected mean"),linewidth=2)+
  geom_point(data=table_ndmh[table_ndmh$group=="Fraction Infected",],aes(x=insect_density,y=-(X5th_pct),color="Fraction infected 5th/95th %iles"),size=5)+
  geom_line(data=table_ndmh[table_ndmh$group=="Fraction Infected",],aes(x=insect_density,y=-(X5th_pct),color="Fraction infected 5th/95th %iles",linetype="Fraction infected 5th/95th %iles"),linewidth=2)+
  geom_point(data=table_ndmh[table_ndmh$group=="Fraction Infected",],aes(x=insect_density,y=-(X95th_pct),color="Fraction infected 5th/95th %iles"),size=5)+
  geom_line(data=table_ndmh[table_ndmh$group=="Fraction Infected",],aes(x=insect_density,y=-(X95th_pct),color="Fraction infected 5th/95th %iles",linetype="Fraction infected 5th/95th %iles"),linewidth=2)+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x)),limits=c(6,100000))+
  scale_y_continuous(limits=c(-2,43),breaks=integer_breaks(),sec.axis = sec_axis(trans=~./1,name="Increase in defoliation (%)",breaks=c(-10,0,10,20,30,40)))+
  geom_point(data=table_ndmh[table_ndmh$group=="Defoliation",],aes(x=insect_density,y=(mean),color="Defoliation mean"),size=5)+
  geom_line(data=table_ndmh[table_ndmh$group=="Defoliation",],aes(x=insect_density,y=(mean),color="Defoliation mean",linetype="Defoliation mean"),linewidth=2)+
  geom_point(data=table_ndmh[table_ndmh$group=="Defoliation",],aes(x=insect_density,y=(X5th_pct),color="Defoliation 5th/95th %iles"),size=5)+
  geom_line(data=table_ndmh[table_ndmh$group=="Defoliation",],aes(x=insect_density,y=(X5th_pct),color="Defoliation 5th/95th %iles",linetype="Defoliation 5th/95th %iles"),linewidth=2)+
  geom_point(data=table_ndmh[table_ndmh$group=="Defoliation",],aes(x=insect_density,y=(X95th_pct),color="Defoliation 5th/95th %iles"),size=5)+
  geom_line(data=table_ndmh[table_ndmh$group=="Defoliation",],aes(x=insect_density,y=(X95th_pct),color="Defoliation 5th/95th %iles",linetype="Defoliation 5th/95th %iles"),linewidth=2)+
  xlab(expression(bold(paste("Initial insect density (egg masses/ha)")) ))+ylab(expression(bold("Relative decrease\n  in infection (%)")))+
  theme(axis.line.y.right = element_line(color = "#cc397b"), 
        axis.ticks.y.right = element_line(color = "#cc397b"),
        axis.text.y.right = element_text(color = "#cc397b"), 
        axis.title.y.right = element_text(color = "#cc397b"))+
  theme(plot.title=element_text(size=32,face="bold"), axis.text=element_text(size=30),
        axis.title=element_text(size=30,face="bold"),plot.background = element_blank(),
        legend.position = "none",legend.title=element_blank(),legend.text=element_text(size=30),legend.background=element_blank(),legend.key=element_blank(),
        strip.background =element_rect(fill="grey50"),strip.text = element_text(colour = 'white'), legend.key.width = unit(3,"cm"))+
  theme(axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 10)))+
  scale_color_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles","Density frequencies"),values=color,name="value type")+
  scale_fill_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles","Density frequencies"),values=color,name="value type")+ 
  scale_linetype_manual(breaks=c("Fraction infected mean","Fraction infected 5th/95th %iles","Defoliation mean","Defoliation 5th/95th %iles"),values=linetype_list,name="value type")+
  facet_grid(. ~ title)+theme(strip.text.x = element_text(size=36,face="bold"))+theme(panel.spacing.x = unit(2, "lines"))+theme(plot.margin = unit(c(1,1,1,3), "cm"))
#scale_x_continuous(trans=log_trans(),breaks=round(exp(-5:-1),3),limits=c(0.004,0.38))+
#facet_wrap(~group)+theme(strip.text.x = element_text(size=30,face="bold"),strip.text.y = element_text(size=30,face="bold"))+theme(panel.spacing.x = unit(2, "lines"))


jpeg(paste("EX_Figure8_summary_ndm.jpg",sep=""),height=6600,width=5800,res=300)
grid.arrange(base,ndm,ndmh,nrow=3)
#grid.arrange(g_IT)

dev.off()   
