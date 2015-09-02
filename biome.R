###########
par(mfrow=c(21,37),mai=c(0,0,0,0), oma=c(0,0,0,0), lwd=0.4)
#par(mai=c(0,0,0,0), oma=c(0,0,0,0), lwd=0.4)
bins<-c(seq(0,21000, by=1000), seq(22000, 50000, by=2000))
bins_image <-c(seq(-500, 21500, by=1000), seq(23000, 49000, by=2000))
bins_length<- c(rep(1000, ))
library(raster)
library(rgdal)
# Upload the climate velocity files
setwd("/Users/afr/Desktop/Koppen_full_cor/")
#velocity_map_year <- stack("climate_change_velocity_perYear_tmp_prec_rescaled_everyKyrs_truncated_0_00005.grd")
kop_list <- dir()
#kop_list_50 <- kop_list[1:37]
koppen <- stack(kop_list)
# Define the directory where the files
Species_color <- read.delim(file.choose(), header=T, stringsAsFactors=F)
Lat_color <- read.delim(file.choose(), header=T, stringsAsFactors=F)
setwd("/Users/afr/Desktop/Evolution_ppt/Corr_result/")
image_climate_biome <- matrix(NA,ncol=42, nrow=21)
colnames(image_climate_biome) <- c("Q0", "Q25", "Q50", "Q75", "Q100" , seq(0, 21000, by=1000), seq(22000, 50000, by=2000))
rownames(image_climate_biome) <- Single_sp
#image_climate_rank <- image_climate
image_climate_n_biome <- image_climate
Single_sp <- unique(Species_color$Species)
#Single_sp <- as.vector(Full_DB_LL$Species)
Single_sp <- c("Mammuthus_primigenius")
for (s in seq_along(Single_sp)){
  temp_DB_climate <- Full_DB_LL[which(Full_DB_LL$Species==Single_sp[s]),]
  temp_points <- as.data.frame(matrix(nrow=nrow(temp_DB_climate), ncol=8))
  colnames(temp_points) <- c("Longitude", "Latitude", "Time_sample", "Time_bin", "Layer", "cell","Biome", "Rec_type")
  matrix_time <- matrix("numeric", nrow=38, ncol=2 )
  matrix_time[,1] <- c(seq(50000, 22000, by=-2000), seq(21000, -1000, by=-1000))
  matrix_time[,2] <- seq(25, 62, by=1)
  for(i in seq_along(temp_DB_climate$Median_Age)){
    for (j in 1:(dim(matrix_time)[1]-1)){
      if(temp_DB_climate$Median_Age[i]<= as.numeric(matrix_time[j]) & temp_DB_climate$Median_Age[i]> as.numeric(matrix_time[j+1])){
        temp_points[i,c(1,2,3,4,5,8)] <- c(temp_DB_climate$Longitude[i], temp_DB_climate$Latitude[i],temp_DB_climate$Median_Age[i], matrix_time[j,1], matrix_time[j,2], ifelse(nchar(temp_DB_climate$Sequence[i]) > 1, "Seq", "Fossil") )
      }
    }
  }
  
  for(k in seq_along(temp_points$Longitude)){
    temp_points[k,c(6,7)]<- extract(koppen, layer=as.numeric(temp_points[k,5]), nl=1, y=matrix(as.numeric(temp_points[k,c(1,2)]), nrow=1, ncol=2), cellnumbers=T)
  }
  if (pie == TRUE){
  for(h in seq_along(bins)){
    temp_hist <- hist((temp_points$Biome[temp_points$Time_bin == bins[h]]), breaks=seq(0,33, by=1), main=bins[h], ylim=c(0,33), plot=F)
######### for pies ###########
      if (sum(temp_hist$counts != 0) > 0){
        temp_pie <- pie(temp_hist$counts, col=colores_biomes, labels=NA, radius=1.05, border=ifelse(((temp_hist$counts != 0) > 1), yes="white", no=NA))
      #text(1,1,sum(temp_hist$counts), cex=1.2)
      }
      else {
        plot(1, type="n", axes=F, xlab="", ylab="")
      }
    }
  }
  points_plot <- temp_points[!is.na(temp_points$Biome),]
  bp_all_points <- boxplot(points_plot$Biome, plot=F)
######### per species box plots ###########
    if (box_plot == TRUE){
      plot.new()
      par(mfrow=c(1,1),mai=c(1,1,1,1), oma=c(0,0,0,0), lwd=0.4)
      bp <- boxplot(points_plot$Biome ~ as.numeric(points_plot$Time_bin), border="#EE760095", boxwex=.7, col="#EE760095", at=tick/1000, las=2, yaxt="n", frame=F, xaxt="n", xaxs="i", yaxs="i", cex=0.5, xlim=c(0,50), ylim=c(0,34))
      tick_labels <- c(0, rep(NA, 4), 5000, rep(NA, 4), 10000, rep(NA, 4), 15000, rep(NA, 4), 20000, rep(NA, 4), 25000, rep(NA, 4), 30000, rep(NA, 4), 35000, rep(NA, 4), 40000, rep(NA, 4), 45000, rep(NA, 4), 50000)
      axis(side=1, at=seq(0, 50, by=1), labels=tick_labels, las=2, cex=0.5)
      axis(side=2, at=seq(0, 33, by=1), cex=0.5)
      mtext("Biome", side=2, line=2, cex.lab=2)
    }
######### all species box plot ###########
  if (box_plot_all == TRUE){
    if (s == 1){
    plot.new()
    par(mfrow=c(1,1),mai=c(1,1,1,1), oma=c(0,0,0,0), lwd=0.4)
    bp <- boxplot(points_plot$Biome ~ as.numeric(points_plot$Time_bin), border=Species_color$Color[which(Species_color$Species == Single_sp[s])], boxwex=.7, col=Species_color$Color[which(Species_color$Species == Single_sp[s])], at=unique(as.numeric(points_plot$Time_bin))/1000, las=2, yaxt="n", frame=F, xaxt="n", xaxs="i", yaxs="i", cex=0.5, xlim=c(0,50), ylim=c(0,34))
    tick_labels <- c(0, rep(NA, 4), 5000, rep(NA, 4), 10000, rep(NA, 4), 15000, rep(NA, 4), 20000, rep(NA, 4), 25000, rep(NA, 4), 30000, rep(NA, 4), 35000, rep(NA, 4), 40000, rep(NA, 4), 45000, rep(NA, 4), 50000)
    axis(side=1, at=seq(0, 50, by=1), labels=tick_labels, las=2, cex=0.5)
    axis(side=2, at=seq(0, 33, by=1), cex=0.5)
    mtext("Biome", side=2, line=2, cex.lab=2)
    }
    if (s > 1){
      bp <- boxplot(add=T,points_plot$Biome ~ as.numeric(points_plot$Time_bin), border=Species_color$Color[which(Species_color$Species == Single_sp[s])], boxwex=.7, col=Species_color$Color[which(Species_color$Species == Single_sp[s])], yaxt="n",at=unique(as.numeric(points_plot$Time_bin))/1000, frame=F, xaxt="n", xaxs="i", yaxs="i", cex=0.5, xlim=c(0,50), ylim=c(0,34))
    }
  } 
}

###########
plot(1, type="n", axes=F, xlab="", ylab="", xlim=c(0,50000), ylim=c(0,33))
for (points in seq_along(points_plot$Biome)){
  for (lat in seq_along(Lat_color$latitude)[-length(Lat_color$latitude)]){
    if (as.numeric(points_plot$Latitude[points]) > Lat_color$latitude[lat] & as.numeric(points_plot$Latitude[points]) < Lat_color$latitude[lat+1]){
      #points(as.numeric(points_plot$Time_bin[points]), points_plot$Biome[points], col=paste(Lat_color$Color[lat], "99", sep=""), pch=19, lwd=7)
      points(as.numeric(points_plot$Time_bin[points]), points_plot$Biome[points], col=Lat_color$Color[lat], pch=19, cex=3)
      text(as.character(Lat_color$latitude[lat]), x=as.numeric(points_plot$Time_bin[points])+0.01, y=points_plot$Biome[points]+0.01, cex=0.7)
    }
  }
}
#######
plot(1, type="n", axes=F, xlab="", ylab="", xlim=c(0,34), ylim=c(20,90))
for (points in seq_along(points_plot$Biome)){
  for (lat in seq_along(Lat_color$latitude)[-length(Lat_color$latitude)]){
    #if (as.numeric(points_plot$Latitude[points]) > Lat_color$latitude[lat] & as.numeric(points_plot$Latitude[points]) < Lat_color$latitude[lat+1]){
      #points(as.numeric(points_plot$Time_bin[points]), points_plot$Biome[points], col=paste(Lat_color$Color[lat], "99", sep=""), pch=19, lwd=7)
      points(points_plot$Time_bin,as.numeric(points_plot$Latitude), col=terrain.colors(33), pch=19, cex=3)
      #text(as.character(Lat_color$latitude[lat]), x=as.numeric(points_plot$Time_bin[points])+0.01, y=points_plot$Biome[points]+0.01, cex=0.7)
    }
  }


plot(points_plot$Time_bin,as.numeric(points_plot$Latitude), col=ifelse(points_plot$Rec_type== "Seq", "#00688B", "#4F4F4F"), pch=19, cex=ifelse(points_plot$Rec_type== "Seq", 1.5, 0.5))



plot(, col=points_plot$Time_bin, pch=19, cex=3)
points(as.numeric(points_plot$Time_bin[points]), points_plot$Biome[points], col=Lat_color$Color[lat], pch=19, cex=3)



