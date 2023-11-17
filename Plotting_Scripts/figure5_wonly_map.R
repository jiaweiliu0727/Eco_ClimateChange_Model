#Function to extract and align the legend 
align_legend <- function(p, hjust = 0.5)
{
  # extract legend
  g <- cowplot::plot_to_gtable(p)
  grobs <- g$grobs
  legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
  legend <- grobs[[legend_index]]
  
  # extract guides table
  guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
  
  # there can be multiple guides within one legend box  
  for (gi in guides_index) {
    guides <- legend$grobs[[gi]]
    
    # add extra column for spacing
    # guides$width[5] is the extra spacing from the end of the legend text
    # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
    # both sides, we get an aligned legend
    spacing <- guides$width[5]
    guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
    guides$widths[6] <- (1-hjust)*spacing
    title_index <- guides$layout$name == "title"
    guides$layout$l[title_index] <- 2
    
    # reconstruct guides and write back
    legend$grobs[[gi]] <- guides
  }
  
  # reconstruct legend and write back
  g$grobs[[legend_index]] <- legend
  g
}

#Function that finds bounds for the diverging color scale.
#The middle point of a color palette is 0.5. It should match the difference of 0%. 0 and 1 are the darkest colors at the two sides.
#h is the color for the highest value, and b is the color for the lowest value
ret_lims <- function(data_vec, all_data, inv = T) {
  vals <- range(all_data)
  higher_mag <- which(abs(vals) == max(abs(vals)))
  if (length(higher_mag) == 2) { #If the upper and lower bound have the same absolute value , e.g. the range is [-5,5], then use the whole color palette.
    if (!inv) {
      b = 0
      h = 1
    } else {
      b = 1
      h = 0
    }
  } else {
    perc <- abs(vals[-higher_mag]) / abs(vals[higher_mag]) #The fraction of the "smaller" side comapred to the "larger" side
    # The "larger" side with a larger absolute value uses the deepest color. The "smaller" side use the "perc"-fold of the deepest color.
    if (!inv) {
      if (vals[1]<0) { b <- ifelse(higher_mag == 1, 0, 0.5 - 0.5 * perc) }
      else{ b <- ifelse(higher_mag == 1, 0, 0.5 + 0.5 * perc) }
      if (vals[2]>0) { h <- ifelse(higher_mag == 1, 0.5 + 0.5 * perc, 1) }
      else{ h <- ifelse(higher_mag == 1, 0.5 - 0.5 * perc, 1) }
    } else {
      if (vals[1]<0) { b <- ifelse(higher_mag == 1, 1, 0.5 + 0.5 * perc) }
      else{ b <- ifelse(higher_mag == 1, 0, 0.5 - 0.5 * perc) }
      if (vals[2]>0) { h <- ifelse(higher_mag == 1, 0.5 - 0.5 * perc, 0) }
      else { h <- ifelse(higher_mag == 1, 0.5 + 0.5 * perc, 0) }
    }
  }
  if (all (data_vec == all_data)) {
    b_curr <- b
    h_curr <- h
  } else { #In this case, multiple data sets share one color bar, so we need to adjust the color range of each data set
    r_vec <- range(data_vec)
    prop <- diff(r_vec) / diff(vals)
    if (prop > 1) {
      print(
        "Are the data subset and all data in the right order? Range of subset greater than that of all data"
      )
      b_curr = 0
      h_curr = 0
    }
    
    if (!inv) {
      h_curr=h-(vals[2]-r_vec[2])/diff(vals)*(h-b)
      b_curr=b+(r_vec[1]-vals[1])/diff(vals) *(h-b)
    } else {
      h_curr=h+(vals[2]-r_vec[2])/diff(vals)*(b-h)
      b_curr=b-(r_vec[1]-vals[1])/diff(vals) *(b-h)
    }
  }
  
  return(c(b_curr, h_curr))
}



library(ggplot2)
library(gridExtra)
library(scales)
library(RColorBrewer)
library(scico)
library(cowplot)


#Read in the file indicating the USA+Canada(Ontario+Quebec) range. The csv file could be obtained by unzipping the rar file.
#USA data obtained by usa=map_data("state",ylim=c(30,70), xlim=c(-92,-50))
#Canada data obtained from https://gadm.org/download_country.html. Here we use the level1 data in version 3.6.
usacan<-read.csv("usacan_gmregion.csv", header = TRUE, sep = ",")

#Get the map of the great lakes
lakes=map_data("lakes")
lakes<-lakes[lakes$region=="Great Lakes",]


## Read in the file including the summary of weather-only model outputs at each location
frame<-read.csv("wonly_highRCP_END_BC_s25_withnorth_hdefof_0.05.csv", header = TRUE, sep = ",")


frame$finfwd=frame$finfwd/frame$finfwh*100 #Calculate the relative difference

#Summary matrix (mean, 5th percentile, 95th percentile, pct.>0)
A=matrix(0,nrow=3,ncol=4)
a=ceiling(0.05*length(frame$finfwh))
b=length(frame$finfwh)-a

A[1,1]=mean(frame$finfwh)
A[1,2]=sort(frame$finfwh,decreasing=FALSE)[a]
A[1,3]=sort(frame$finfwh,decreasing=FALSE)[b]
A[1,4]=var(frame$finfwh)

A[2,1]=mean(frame$finfwf)
A[2,2]=sort(frame$finfwf,decreasing=FALSE)[a]
A[2,3]=sort(frame$finfwf,decreasing=FALSE)[b]
A[2,4]=var(frame$finfwf)

A[3,1]=mean(frame$finfwd)
A[3,2]=sort(frame$finfwd,decreasing=FALSE)[a]
A[3,3]=sort(frame$finfwd,decreasing=FALSE)[b]
A[3,4]=length(frame$finfwd[frame$finfwd>0])/length(frame$finfwd)*100


bottom=min(frame$finfwd)
top=max(frame$finfwd)
for (i in 1:length(frame$finfwd)){
  if (frame$finfwd[i]<bottom){
    frame$finfwd[i]=bottom
  }
  if (frame$finfwd[i]>top){
    frame$finfwd[i]=top
  }
}


memory.limit(size = 22000)
lims_fracinf <- ret_lims(frame$finfwd, frame$finfwd, inv = T)

finfw3<-ggplot() + geom_polygon(data = usacan, aes(x=long, y = lat, group = group),fill=NA,color="black") +coord_fixed(1.3)+xlab("Longitude") + ylab("Latitude")

finfw3<-finfw3+coord_map(xlim=c(-90.5,-65),ylim = c(35, 47))
finfw3<-finfw3+geom_tile(data=frame,aes(x=lon,y=lat,fill=finfwd),height=0.18,width=0.18)+scale_fill_scico(palette = "bam", begin = lims_fracinf[1], end = lims_fracinf[2])+ theme_bw()+ guides(fill = guide_colorbar(title="Relative \nDifference (%)"))+theme(legend.title.align=0.5)
finfw3<-finfw3+geom_polygon(data = usacan, aes(x=long, y = lat, group = group),fill=NA,color="black") +coord_fixed(1.3) #Draw the boarder lines again in case some are covered by the color dots
finfw3<-finfw3+coord_map(xlim=c(-90.5,-65),ylim = c(35, 47))
finfw3<-finfw3+geom_polygon(data = lakes, aes(x=long, y = lat, group = group),fill="white",color="black") + coord_fixed(1.3) #Make sure that the lake area is in white color
finfw3<-finfw3+coord_map(xlim=c(-90.5,-65),ylim = c(35, 47))

finfw3<-finfw3+scale_x_continuous(breaks=seq(-100,-70,10))
finfw3<-finfw3+theme(plot.title=element_text(size=24,face="bold"), axis.text=element_text(size=20),
                     axis.title=element_text(size=24,face="bold"),
                     legend.title=element_text(size=16,face="bold"),
                     legend.text=element_text(size=16),legend.key.size=unit(2, "cm"))

#Drawing the histogram of relative difference and showing the mean and 5th/95th percentiles
frame<-read.csv("wonly_highRCP_END_BC_s25_withnorth_hdefof_0.05.csv", header = TRUE, sep = ",")
frame$finfwd=frame$finfwd/frame$finfwh*100
g<-ggplot(frame, aes(x=finfwd)) + geom_histogram(aes(y=..count../sum(..count..)*100),color="black", fill="white",breaks=seq(-30,20,1))+scale_x_continuous(limits=c(-30,20))+theme_classic()
g<-g+xlab("Relative Diff (%)")+ylab("Pct. sites (%)")
maxcount=max(ggplot_build(g)$data[[1]]$count)/length(frame$finfwd)*100   #showing the maximum count in histogram
g<-g+annotate("segment", x = signif(A[3,1],3), xend = signif(A[3,1],3), y = 0, yend = maxcount+2)
g<-g+annotate("segment", x = signif(A[3,2],3), xend = signif(A[3,2],3), y = 0, yend = maxcount+2, col="red")
g<-g+annotate("segment", x = signif(A[3,3],3), xend = signif(A[3,3],3), y = 0, yend = maxcount+2, col="red")


g<-g+annotate("text",x=A[3,1],y=maxcount+3.5,label = paste(signif(A[3,1],3)),size=5)
g<-g+annotate("text",col="red",x=A[3,2]-5,y=maxcount+3.5,label = paste(signif(A[3,2],3)),size=5)
g<-g+annotate("text",col="red",x=A[3,3]+4,y=maxcount+3.5,label = "0.37",size=5)

g<-g+annotate("text",x=A[3,1]+1,y=maxcount+6,label = paste("Pct. sites >0: ",signif(A[3,4],3),"%"),size=5)

g<-g+theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))

w25 <-ggdraw(align_legend(finfw3)) +draw_plot(g, x = 0.58, y = .24, width = .25, height = .25)



pdf("figure5.pdf", width = 12, height = 9)
grid.arrange(w25)

dev.off()
gc()






