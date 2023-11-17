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


# Read in the file including the summary of weather variables across epizootic times at each location.
# The locations are in a subset of all the 19800 grids, inside the usacan region, with defoliation observed in the past.
frame<-read.csv("epiweather_highRCP_END_BC_withnorth_whole_hdefof_0.05.csv", header = TRUE, sep = ",")

#Rainfall
frame$epiraindiff=frame$epiraindiff/frame$epirainh*100 #Calculate the relative difference
#Summary matrix (mean, 5th percentile, 95th percentile, pct.>0)
A=matrix(0,nrow=3,ncol=4)
a=ceiling(0.05*length(frame$epiraindiff))
b=length(frame$epiraindiff)-a

A[1,1]=mean(frame$epiraindiff)
A[1,2]=sort(frame$epiraindiff,decreasing=FALSE)[a]
A[1,3]=sort(frame$epiraindiff,decreasing=FALSE)[b]
A[1,4]=length(frame$epiraindiff[frame$epiraindiff>0])/length(frame$epiraindiff)*100

A[2,1]=mean(frame$epirhdiff)
A[2,2]=sort(frame$epirhdiff,decreasing=FALSE)[a]
A[2,3]=sort(frame$epirhdiff,decreasing=FALSE)[b]
A[2,4]=length(frame$epirhdiff[frame$epirhdiff>0])/length(frame$epirhdiff)*100

A[3,1]=mean(frame$epiavetdiff)
A[3,2]=sort(frame$epiavetdiff,decreasing=FALSE)[a]
A[3,3]=sort(frame$epiavetdiff,decreasing=FALSE)[b]
A[3,4]=length(frame$epiavetdiff[frame$epiavetdiff>0])/length(frame$epiavetdiff)*100


bottom=min(frame$epiraindiff)
top=max(frame$epiraindiff)

for (i in 1:length(frame$epiraindiff)){
  if (frame$epiraindiff[i]>top){
    frame$epiraindiff[i]=top
  }
  if (frame$epiraindiff[i]<bottom){
    frame$epiraindiff[i]=bottom
  }
}
#Use the ret_lims function for the fill column
#Change inv (T/F) in order to reverse the color scale but keep the correct color balancing
#The color scale I've used in these plots (bam) goes from purple (0) to green (1), so we want to set inv=T for color to get purple with increasing values
lims_rainfall <- ret_lims(frame$epiraindiff, frame$epiraindiff, inv = T)

f<-ggplot() + geom_polygon(data = usacan, aes(x=long, y = lat, group = group),fill=NA,color="black")+xlab("Longitude") + ylab("Latitude")

f<-f+coord_map(xlim=c(-90.5,-65),ylim = c(35, 47))

#line that uses the lims_rainfall to actually apply the color scale adjusted for thew relatiuve magnitude of min and max using scale_fill_scico. begin and end should always be lims[1] and lims[2], no need to change anything but the lims vector name here.
f<-f+geom_tile(data=frame,aes(x=lon,y=lat,fill=epiraindiff),height=0.18,width=0.18)+
  scale_fill_scico(palette = "bam", begin = lims_rainfall[1], end = lims_rainfall[2])  + 
  theme_bw() + 
  guides(fill = guide_colorbar(title="Relative \nDifference (%)"))+
  theme(legend.title.align=0.5)

f<-f+geom_polygon(data = usacan, aes(x=long, y = lat, group = group),fill=NA,color="black") +coord_fixed(1.3) #Draw the boarder lines again in case some are covered by the color dots
f<-f+coord_map(xlim=c(-90.5,-65),ylim = c(35, 47))
f<-f+geom_polygon(data = lakes, aes(x=long, y = lat, group = group),fill="white",color="black") + coord_fixed(1.3) #Make sure that the lake area is in white color
f<-f+coord_map(xlim=c(-90.5,-65),ylim = c(35, 47))
f<-f+ggtitle("Rainfall")+theme(plot.title = element_text(hjust = 0.5))
f<-f+theme(plot.title=element_text(size=24,face="bold"), axis.text=element_text(size=20),
           axis.title=element_text(size=24,face="bold"),
           legend.title=element_text(size=16,face="bold"),
           legend.text=element_text(size=16),legend.key.size=unit(2, "cm"))

#Drawing the histogram of relative difference and showing the mean and 5th/95th percentiles
frame<-read.csv("epiweather_highRCP_END_BC_withnorth_whole_hdefof_0.05.csv", header = TRUE, sep = ",")
frame$epiraindiff=frame$epiraindiff/frame$epirainh*100
g<-ggplot(frame, aes(x=epiraindiff)) + geom_histogram(aes(y=..count../sum(..count..)*100),color="black", fill="white",bins=50)+theme_classic()
g<-g+xlab("Relative.Diff (%)")+ylab("Pct. sites (%)")+scale_y_continuous(limits=c(0,14))
maxcount=max(ggplot_build(g)$data[[1]]$count)/length(frame$epiraindiff)*100   #showing the maximum count in histogram
g<-g+annotate("segment", x = signif(A[1,1],3), xend = signif(A[1,1],3), y = 0, yend = maxcount+2)
g<-g+annotate("segment", x = signif(A[1,2],3), xend = signif(A[1,2],3), y = 0, yend = maxcount+2, col="red")
g<-g+annotate("segment", x = signif(A[1,3],3), xend = signif(A[1,3],3), y = 0, yend = maxcount+2, col="red")


g<-g+annotate("text",x=A[1,1]+4.5,y=maxcount+3.5,label = paste(signif(A[1,1],3)),size=5)
g<-g+annotate("text",col="red",x=A[1,2]-5,y=maxcount+3.5,label = paste(signif(A[1,2],3)),size=5)
g<-g+annotate("text",col="red",x=A[1,3]+5,y=maxcount+3.5,label = paste(signif(A[1,3],3)),size=5)

g<-g+annotate("text",x=A[3,1]+10,y=maxcount+6.5,label = paste("Pct. sites >0: ",signif(A[1,4],3),"%"),size=5)
g<-g+theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))
r <-ggdraw(align_legend(f)) +draw_plot(g, x = 0.565, y = .20, width = .25, height = .25)

#Relative humidity
frame<-read.csv("epiweather_highRCP_END_BC_withnorth_whole_hdefof_0.05.csv", header = TRUE, sep = ",")
frame$epirhdiff=frame$epirhdiff/frame$epirhh*100
#Summary matrix
A=matrix(0,nrow=3,ncol=4)
a=ceiling(0.05*length(frame$epiraindiff))
b=length(frame$epiraindiff)-a

A[1,1]=mean(frame$epiraindiff)
A[1,2]=sort(frame$epiraindiff,decreasing=FALSE)[a]
A[1,3]=sort(frame$epiraindiff,decreasing=FALSE)[b]
A[1,4]=length(frame$epiraindiff[frame$epiraindiff>0])/length(frame$epiraindiff)*100

A[2,1]=mean(frame$epirhdiff)
A[2,2]=sort(frame$epirhdiff,decreasing=FALSE)[a]
A[2,3]=sort(frame$epirhdiff,decreasing=FALSE)[b]
A[2,4]=length(frame$epirhdiff[frame$epirhdiff>0])/length(frame$epirhdiff)*100

A[3,1]=mean(frame$epiavetdiff)
A[3,2]=sort(frame$epiavetdiff,decreasing=FALSE)[a]
A[3,3]=sort(frame$epiavetdiff,decreasing=FALSE)[b]
A[3,4]=length(frame$epiavetdiff[frame$epiavetdiff>0])/length(frame$epiavetdiff)*100


bottom=min(frame$epirhdiff)
top=max(frame$epirhdiff)


for (i in 1:length(frame$epirhdiff)){
  if (frame$epirhdiff[i]>top){
    frame$epirhdiff[i]=top
  }
  if (frame$epirhdiff[i]<bottom){
    frame$epirhdiff[i]=bottom
  }
}

lims_rh <- ret_lims(frame$epirhdiff, frame$epirhdiff, inv = T)
f<-ggplot() + geom_polygon(data = usacan, aes(x=long, y = lat, group = group),fill=NA,color="black")+xlab("Longitude") + ylab("Latitude")

f<-f+coord_map(xlim=c(-90.5,-65),ylim = c(35, 47))
f<-f+geom_tile(data=frame,aes(x=lon,y=lat,fill=epirhdiff),height=0.18,width=0.18) + scale_fill_scico(palette = "bam", begin = lims_rh[1], end = lims_rh[2]) + theme_bw() + guides(fill = guide_colorbar(title="Relative \nDifference (%)"))+theme(legend.title.align=0.5)
f<-f+geom_polygon(data = usacan, aes(x=long, y = lat, group = group),fill=NA,color="black") +coord_fixed(1.3)
f<-f+coord_map(xlim=c(-90.5,-65),ylim = c(35, 47))
f<-f+geom_polygon(data = lakes, aes(x=long, y = lat, group = group),fill="white",color="black") + coord_fixed(1.3)
f<-f+coord_map(xlim=c(-90.5,-65),ylim = c(35, 47))
f<-f+ggtitle("Relative Humidity")+theme(plot.title = element_text(hjust = 0.5))


f<-f+theme(plot.title=element_text(size=24,face="bold"), axis.text=element_text(size=20),
           axis.title=element_text(size=24,face="bold"),
           legend.title=element_text(size=16,face="bold"),
           legend.text=element_text(size=16),legend.key.size=unit(2, "cm"))

frame<-read.csv("epiweather_highRCP_END_BC_withnorth_whole_hdefof_0.05.csv", header = TRUE, sep = ",")
frame$epirhdiff=frame$epirhdiff/frame$epirhh*100
g<-ggplot(frame, aes(x=epirhdiff)) + geom_histogram(aes(y=..count../sum(..count..)*100),color="black", fill="white",bins=50)+theme_classic()
g<-g+xlab("Relative.Diff (%)")+ylab("Pct. sites (%)")+scale_y_continuous(limits=c(0,14.5))
maxcount=max(ggplot_build(g)$data[[1]]$count)/length(frame$epiraindiff)*100   #showing the maximum count in histogram
g<-g+annotate("segment", x = signif(A[2,1],3), xend = signif(A[2,1],3), y = 0, yend = maxcount+1.5)
g<-g+annotate("segment", x = signif(A[2,2],3), xend = signif(A[2,2],3), y = 0, yend = maxcount+1.5, col="red")
g<-g+annotate("segment", x = signif(A[2,3],3), xend = signif(A[2,3],3), y = 0, yend = maxcount+1.5, col="red")




g<-g+annotate("text",x=A[2,1]+1.2,y=maxcount+3.5,label = paste(signif(A[2,1],3)),size=5)
g<-g+annotate("text",col="red",x=A[2,2],y=maxcount+3.5,label = "-6.10",size=5)
g<-g+annotate("text",col="red",x=A[2,3]+0.5,y=maxcount+3.5,label = "1.20",size=5)

g<-g+annotate("text",x=A[2,1]+1.8,y=maxcount+6,label = paste("Pct. sites >0: ",signif(A[2,4],3),"%"),size=5)
g<-g+theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))
rh <-ggdraw(align_legend(f)) +draw_plot(g, x = 0.565, y = .20, width = .25, height = .25)


#Average temperature
frame<-read.csv("epiweather_highRCP_END_BC_withnorth_whole_hdefof_0.05.csv", header = TRUE, sep = ",")
frame$epiavetdiff=frame$epiavetdiff/frame$epiaveth*100
#Summary matrix
A=matrix(0,nrow=3,ncol=4)
a=ceiling(0.05*length(frame$epiraindiff))
b=length(frame$epiraindiff)-a

A[1,1]=mean(frame$epiraindiff)
A[1,2]=sort(frame$epiraindiff,decreasing=FALSE)[a]
A[1,3]=sort(frame$epiraindiff,decreasing=FALSE)[b]
A[1,4]=length(frame$epiraindiff[frame$epiraindiff>0])/length(frame$epiraindiff)*100

A[2,1]=mean(frame$epirhdiff)
A[2,2]=sort(frame$epirhdiff,decreasing=FALSE)[a]
A[2,3]=sort(frame$epirhdiff,decreasing=FALSE)[b]
A[2,4]=length(frame$epirhdiff[frame$epirhdiff>0])/length(frame$epirhdiff)*100

A[3,1]=mean(frame$epiavetdiff)
A[3,2]=sort(frame$epiavetdiff,decreasing=FALSE)[a]
A[3,3]=sort(frame$epiavetdiff,decreasing=FALSE)[b]
A[3,4]=length(frame$epiavetdiff[frame$epiavetdiff>0])/length(frame$epiavetdiff)*100

a=ceiling(0.005*length(frame$epiavetdiff))
b=length(frame$epiavetdiff)-a
bottom<-sort(frame$epiavetdiff,decreasing=FALSE)[a]
top<-sort(frame$epiavetdiff,decreasing=FALSE)[b]
bottom=min(frame$epiavetdiff)
top=max(frame$epiavetdiff)

for (i in 1:length(frame$epiavetdiff)){
  if (frame$epiavetdiff[i]>top){
    frame$epiavetdiff[i]=top
  }
  if (frame$epiavetdiff[i]<bottom){
    frame$epiavetdiff[i]=bottom
  }
}


lims_epi <- ret_lims(frame$epiavetdiff, frame$epiavetdiff, inv = T)
f<-ggplot() + geom_polygon(data = usacan, aes(x=long, y = lat, group = group),fill=NA,color="black")+xlab("Longitude") + ylab("Latitude")

f<-f+coord_map(xlim=c(-90.5,-65),ylim = c(35, 47))
f<-f+geom_tile(data=frame,aes(x=lon,y=lat,fill=epiavetdiff),height=0.18,width=0.18) + scale_fill_scico(palette = "bam", begin = lims_epi[1], end = lims_epi[2]) + theme_bw() + guides(fill = guide_colorbar(title="Relative \nDifference (%)"))+theme(legend.title.align=0.5)
f<-f+geom_polygon(data = usacan, aes(x=long, y = lat, group = group),fill=NA,color="black") +coord_fixed(1.3)
f<-f+coord_map(xlim=c(-90.5,-65),ylim = c(35, 47))
f<-f+geom_polygon(data = lakes, aes(x=long, y = lat, group = group),fill="white",color="black") + coord_fixed(1.3)
f<-f+coord_map(xlim=c(-90.5,-65),ylim = c(35, 47))
f<-f+ggtitle("Temperature")+theme(plot.title = element_text(hjust = 0.5))



f<-f+theme(plot.title=element_text(size=24,face="bold"), axis.text=element_text(size=20),
           axis.title=element_text(size=24,face="bold"),
           legend.title=element_text(size=16,face="bold"),
           legend.text=element_text(size=16),legend.key.size=unit(2, "cm"))

frame<-read.csv("epiweather_highRCP_END_BC_withnorth_whole_hdefof_0.05.csv", header = TRUE, sep = ",")
frame$epiavetdiff=frame$epiavetdiff/frame$epiaveth*100


g<-ggplot(frame, aes(x=epiavetdiff)) + geom_histogram(aes(y=..count../sum(..count..)*100),color="black", fill="white",bins=50)+theme_classic()
g<-g+xlab("Relative.Diff (%)")+ylab("Pct. sites (%)")+scale_y_continuous(limits=c(0,13))
maxcount=max(ggplot_build(g)$data[[1]]$count)/length(frame$epiraindiff)*100   #showing the maximum count in histogram
g<-g+annotate("segment", x = signif(A[3,1],3), xend = signif(A[3,1],3), y = 0, yend = maxcount+1.5)
g<-g+annotate("segment", x = signif(A[3,2],3), xend = signif(A[3,2],3), y = 0, yend = maxcount+1.5, col="red")
g<-g+annotate("segment", x = signif(A[3,3],3), xend = signif(A[3,3],3), y = 0, yend = maxcount+1.5, col="red")



g<-g+annotate("text",x=A[3,1],y=maxcount+3,label = paste(signif(A[3,1],3)),size=5)
g<-g+annotate("text",col="red",x=A[3,2]-1,y=maxcount+3,label = paste(signif(A[3,2],3)),size=5)
g<-g+annotate("text",col="red",x=A[3,3]+1,y=maxcount+3,label = paste(signif(A[3,3],3)),size=5)

g<-g+annotate("text",x=A[3,1],y=maxcount+6,label = paste("Pct. sites >0: ",signif(A[3,4],3),"%"),size=5)
g<-g+theme(axis.text=element_text(size=16),axis.title=element_text(size=16,face="bold"))
t <-ggdraw(align_legend(f))+draw_plot(g, x = 0.565, y = .20, width = .25, height = .25)

b<-ggplot() + theme_void() #Blank plot

pdf("figure1.pdf", width = 22, height = 15)
grid.arrange(r, rh, t,b,layout_matrix = rbind(c(1,1,2,2), c(4,3,3,4)))


dev.off()
gc()

