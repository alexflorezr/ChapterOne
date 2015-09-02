
#########################################################################################################################
# Climate velocity
#########################################################################################################################

bins<-c(seq(0, 21000, by=1000), seq(22000, 50000, by=2000))
bins_image <-c(seq(-500, 21500, by=1000), seq(23000, 49000, by=2000))
bins_length<- c(rep(1000, ))
library(raster)
# Upload the climate velocity files
setwd("/Users/afr/Desktop/CHECK/Ditte_files/Velocity/paleo_cc_velocity/")
velocity_map_year <- stack("climate_change_velocity_perYear_tmp_prec_rescaled_everyKyrs_truncated_0_00005.grd")
# Define the directory where the files
setwd("/Users/afr/Desktop/Evolution_ppt/Corr_result/")
image_climate <- matrix(NA,ncol=42, nrow=21)
colnames(image_climate) <- c("Q0", "Q25", "Q50", "Q75", "Q100" , seq(0, 21000, by=1000), seq(22000, 50000, by=2000))
rownames(image_climate) <- Single_sp
image_climate_rank <- image_climate
image_climate_n <- image_climate
#Single_sp <- unique(Full_DB_LL$Species)
Single_sp <- c("Mammuthus_primigenius","Coelodonta_antiquitatis","Ovibos_moschatus", "Dicrostonyx_torquatus","Panthera_leo")
for (s in seq_along(Single_sp)){
  temp_DB_climate <- Full_DB_LL[which(Full_DB_LL$Species==Single_sp[s]),]
  temp_points <- as.data.frame(matrix(nrow=nrow(temp_DB_climate), ncol=7))
  colnames(temp_points) <- c("Longitude", "Latitude", "Time_sample", "Time_bin", "Layer", "cell","Velocity")
  matrix_time <- matrix("numeric", nrow=38, ncol=2 )
  matrix_time[,1] <- c(seq(50000, 22000, by=-2000), seq(21000, -1000, by=-1000))
  matrix_time[,2] <- seq(25, 62, by=1)
  for(i in seq_along(temp_DB_climate$Median_Age)){
    for (j in 1:(dim(matrix_time)[1]-1)){
      if(temp_DB_climate$Median_Age[i]<= as.numeric(matrix_time[j]) & temp_DB_climate$Median_Age[i]> as.numeric(matrix_time[j+1])){
        temp_points[i,c(1,2,3,4,5)] <- c(temp_DB_climate$Longitude[i], temp_DB_climate$Latitude[i],temp_DB_climate$Median_Age[i], matrix_time[j,1], matrix_time[j,2] )
      }
    }
  }
  for(k in seq_along(temp_points$Longitude)){
    temp_points[k,c(6,7)]<- extract(velocity_map_year, layer=as.numeric(temp_points[k,5]), nl=1, y=matrix(as.numeric(temp_points[k,c(1,2)]), nrow=1, ncol=2), cellnumbers=T)
  }
  points_plot <- temp_points[!is.na(temp_points$Velocity),]
  bp_all_points <- boxplot(points_plot$Velocity, plot=F)
  #boxplot(points_plot$Velocity)
  bp <- boxplot(points_plot$Velocity ~ as.numeric(points_plot$Time_bin), plot=F)
  tick <- as.vector(as.numeric(bp$names))
  # Check point
  sort(as.numeric(points_plot$Time_bin))
  # Save the velocity values for every species
### setwd("/Users/afr/Desktop/kk_temp/")
### write.table(x=temp_points, file=paste(tolower(Single_sp[s]), "_vel_bin.txt", sep=""), sep="\t",row.names=F)
  species_sp <- paste(strsplit(strsplit(tolower(Single_sp[s]), split = "_")[[1]][1], split = "")[[1]][1], strsplit(strsplit(tolower(Single_sp[s]), split = "_")[[1]][2], split = "")[[1]][1], sep="")
  bp <- boxplot(points_plot$Velocity ~ as.numeric(points_plot$Time_bin), main=Single_sp[s], border="#EE760095", boxwex=.7, col="#EE760095", at=tick/1000, las=2, yaxt="n", frame=F, xaxt="n",  cex=0.5, xlim=c(0,50))
  tick_labels <- c(0, rep(NA, 4), 5000, rep(NA, 4), 10000, rep(NA, 4), 15000, rep(NA, 4), 20000, rep(NA, 4), 25000, rep(NA, 4), 30000, rep(NA, 4), 35000, rep(NA, 4), 40000, rep(NA, 4), 45000, rep(NA, 4), 50000)
 axis(side=1, at=seq(0, 50, by=1), labels=tick_labels, las=2, cex=0.5)
 axis(side=4, at=seq(0, 1, by=0.2), cex=0.5)
#lines(as.numeric(bp$names)/1000, bp$stats[1,], col="red", lwd=2)
#lines(as.numeric(bp$names)/1000, bp$stats[2,], col="pink", lwd=2)
#lines(as.numeric(bp$names)/1000, bp$stats[3,], col="grey", lwd=2)
#lines(as.numeric(bp$names)/1000, bp$stats[4,], col="lightblue", lwd=2)
#lines(as.numeric(bp$names)/1000, bp$stats[5,], col="green", lwd=2)
  image_climate[s,as.vector(na.omit(match(bp$names, colnames(image_climate))))] <- bp$stats[3,]
  image_climate[s, c(1,2,3,4,5)] <-as.vector(bp_all_points$stats)
  image_climate_rank[s,as.vector(na.omit(match(bp$names, colnames(image_climate))))] <- rank(bp$stats[3,])
  image_climate_n[s,as.vector(na.omit(match(bp$names, colnames(image_climate))))] <- bp$n                                                                                    

# Plot velocity profile per time bin
setwd("/Users/afr/Desktop/Evolution_ppt")
par(mar=c(6,6,6,6))
bsp_raw <- na.omit(read.delim(paste(species_sp, "_BSP_50_data.txt", sep=""), header=T, stringsAsFactors=F, sep="\t"))
#bsp_raw <- na.omit(read.delim(paste(species_vector[s], "_BSP_50_data.txt", sep=""), header=T, stringsAsFactors=F, sep="\t"))
plot(bsp_raw$Time, log(bsp_raw$Median), main=Single_sp[s], ylim=c(4,16), lwd=3, type="l",frame.plot=F, lab=c(10,10,5), ylab="Population size", xlab="", col="#008B45", cex=0.5, xaxt="n")
lines(bsp_raw$Time, log(bsp_raw$Mean), lty=2)
par(new=TRUE)
temp_time_bsp <- as.data.frame(matrix("numeric", nrow=37, ncol=3), stringsAsFactors=F)
colnames(temp_time_bsp) <- c("Time_bin", "Average_popSize", "Slope")
temp_time_bsp$Time_bin <- c(seq(50000, 22000, by=-2000), seq(21000, 0, by=-1000))
# Estimate the average pop size
for (t in 1:(length(temp_time_bsp$Time_bin)-1)){
  temp_time_bsp$Average_popSize[t] <- mean.default(bsp_raw[which(bsp_raw$Time <= temp_time_bsp$Time_bin[t] & bsp_raw$Time >= temp_time_bsp$Time_bin[t+1]),5])
}
# Estimate the slope of the average pop size
for (t in 1:(length(temp_time_bsp$Time_bin)-2)){ 
  temp_time_bsp$Slope[t] <- ((log(as.numeric(temp_time_bsp$Average_popSize[t+1])) - log(as.numeric(temp_time_bsp$Average_popSize[t])))/ (as.numeric(temp_time_bsp$Time_bin[t+1]) - as.numeric(temp_time_bsp$Time_bin[t])))
}
bp <- boxplot(points_plot$Velocity ~ as.numeric(points_plot$Time_bin), border="#EE760095", boxwex=.7, col="#EE760095", at=tick/1000, las=2, yaxt="n", frame=F, xaxt="n",  cex=0.5, xlim=c(0,50))
tick_labels <- c(0, rep(NA, 4), 5000, rep(NA, 4), 10000, rep(NA, 4), 15000, rep(NA, 4), 20000, rep(NA, 4), 25000, rep(NA, 4), 30000, rep(NA, 4), 35000, rep(NA, 4), 40000, rep(NA, 4), 45000, rep(NA, 4), 50000)
axis(side=1, at=seq(0, 50, by=1), labels=tick_labels, las=2, cex=0.5)
axis(side=4, at=seq(0, 1, by=0.2), cex=0.5)
mtext("Velocity", side=4, line=3, cex.lab=2, col="#EE760095")
#text(tick/1000, rep(-0.001, length(tick)) ,labels=as.character(bp$n), cex=0.4)
lines(as.numeric(bp$names)/1000, bp$stats[3,], col="#EE7600", lwd=3)
par(new=TRUE)
plot(temp_time_bsp$Time_bin, temp_time_bsp$Slope, type="l", lwd=3, col="#595959",axes=F, ylab="", xlab="")
#axis(side=2, cex=0.5)
#mtext("Slope population size", side=4, line=3, cex.lab=1, col="red")
}


#########################################################################################################################
# Koppen biomes
#########################################################################################################################

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
setwd("/Users/afr/Desktop/Evolution_ppt/Corr_result/")
image_climate_biome <- matrix(NA,ncol=42, nrow=21)
colnames(image_climate_biome) <- c("Q0", "Q25", "Q50", "Q75", "Q100" , seq(0, 21000, by=1000), seq(22000, 50000, by=2000))
rownames(image_climate_biome) <- Single_sp
#image_climate_rank <- image_climate
image_climate_n_biome <- image_climate
Single_sp <- unique(Full_DB_LL$Species)
#Single_sp <- c("Mammuthus_primigenius")
for (s in seq_along(Single_sp)){
  temp_DB_climate <- Full_DB_LL[which(Full_DB_LL$Species==Single_sp[s]),]
  temp_points <- as.data.frame(matrix(nrow=nrow(temp_DB_climate), ncol=7))
  colnames(temp_points) <- c("Longitude", "Latitude", "Time_sample", "Time_bin", "Layer", "cell","Biome")
  matrix_time <- matrix("numeric", nrow=38, ncol=2 )
  matrix_time[,1] <- c(seq(50000, 22000, by=-2000), seq(21000, -1000, by=-1000))
  matrix_time[,2] <- seq(25, 62, by=1)
  for(i in seq_along(temp_DB_climate$Median_Age)){
    for (j in 1:(dim(matrix_time)[1]-1)){
      if(temp_DB_climate$Median_Age[i]<= as.numeric(matrix_time[j]) & temp_DB_climate$Median_Age[i]> as.numeric(matrix_time[j+1])){
        temp_points[i,c(1,2,3,4,5)] <- c(temp_DB_climate$Longitude[i], temp_DB_climate$Latitude[i],temp_DB_climate$Median_Age[i], matrix_time[j,1], matrix_time[j,2] )
      }
    }
  }
  for(k in seq_along(temp_points$Longitude)){
    temp_points[k,c(6,7)]<- extract(koppen, layer=as.numeric(temp_points[k,5]), nl=1, y=matrix(as.numeric(temp_points[k,c(1,2)]), nrow=1, ncol=2), cellnumbers=T)
  }
  points_plot <- temp_points[!is.na(temp_points$Biome),]
  bp_all_points <- boxplot(points_plot$Biome, plot=F)
  #boxplot(points_plot$Velocity)
  bp <- boxplot(points_plot$Biome ~ as.numeric(points_plot$Time_bin), plot=F)
  tick <- as.vector(as.numeric(bp$names))
  # Check point
  sort(as.numeric(points_plot$Time_bin))
  # Save the velocity values for every species
  ### setwd("/Users/afr/Desktop/kk_temp/")
  ### write.table(x=temp_points, file=paste(tolower(Single_sp[s]), "_vel_bin.txt", sep=""), sep="\t",row.names=F)
  species_sp <- paste(strsplit(strsplit(tolower(Single_sp[s]), split = "_")[[1]][1], split = "")[[1]][1], strsplit(strsplit(tolower(Single_sp[s]), split = "_")[[1]][2], split = "")[[1]][1], sep="")
  bp <- boxplot(points_plot$Biome ~ as.numeric(points_plot$Time_bin), main=Single_sp[s], border="#EE760095", boxwex=.7, col="#EE760095", at=tick/1000, las=2, yaxt="n", frame=F, xaxt="n",  cex=0.5, xlim=c(0,50))
  tick_labels <- c(0, rep(NA, 4), 5000, rep(NA, 4), 10000, rep(NA, 4), 15000, rep(NA, 4), 20000, rep(NA, 4), 25000, rep(NA, 4), 30000, rep(NA, 4), 35000, rep(NA, 4), 40000, rep(NA, 4), 45000, rep(NA, 4), 50000)
  axis(side=1, at=seq(0, 50, by=1), labels=tick_labels, las=2, cex=0.5)
  axis(side=2, at=seq(0, 32, by=1), cex=0.5)
  #lines(as.numeric(bp$names)/1000, bp$stats[1,], col="red", lwd=2)
  #lines(as.numeric(bp$names)/1000, bp$stats[2,], col="pink", lwd=2)
  #lines(as.numeric(bp$names)/1000, bp$stats[3,], col="grey", lwd=2)
  #lines(as.numeric(bp$names)/1000, bp$stats[4,], col="lightblue", lwd=2)
  #lines(as.numeric(bp$names)/1000, bp$stats[5,], col="green", lwd=2)
  image_climate_biome[s,as.vector(na.omit(match(bp$names, colnames(image_climate))))] <- bp$stats[3,]
  image_climate_biome[s, c(1,2,3,4,5)] <-as.vector(bp_all_points$stats)
  #image_climate_rank[s,as.vector(na.omit(match(bp$names, colnames(image_climate))))] <- rank(bp$stats[3,])
  image_climate_n_biome[s,as.vector(na.omit(match(bp$names, colnames(image_climate))))] <- bp$n                                                                                    
  
setwd("/Users/afr/Desktop/Evolution_ppt")
par(mar=c(6,6,6,6))
#bsp_raw <- na.omit(read.delim(paste(species_sp, "_BSP_50_data.txt", sep=""), header=T, stringsAsFactors=F, sep="\t"))
#bsp_raw <- na.omit(read.delim(paste(species_vector[s], "_BSP_50_data.txt", sep=""), header=T, stringsAsFactors=F, sep="\t"))
#plot(bsp_raw$Time, log(bsp_raw$Median), main=Single_sp[s], ylim=c(4,16), lwd=3, type="l",frame.plot=F, lab=c(10,10,5), ylab="Population size", xlab="", col="#008B45", cex=0.5, xaxt="n")
#lines(bsp_raw$Time, log(bsp_raw$Mean), lty=2)
#par(new=TRUE)
#temp_time_bsp <- as.data.frame(matrix("numeric", nrow=37, ncol=3), stringsAsFactors=F)
#colnames(temp_time_bsp) <- c("Time_bin", "Average_popSize", "Slope")
#temp_time_bsp$Time_bin <- c(seq(50000, 22000, by=-2000), seq(21000, 0, by=-1000))
# Estimate the average pop size
#for (t in 1:(length(temp_time_bsp$Time_bin)-1)){
#  temp_time_bsp$Average_popSize[t] <- mean.default(bsp_raw[which(bsp_raw$Time <= temp_time_bsp$Time_bin[t] & bsp_raw$Time >= temp_time_bsp$Time_bin[t+1]),5])
#}
# Estimate the slope of the average pop size
#for (t in 1:(length(temp_time_bsp$Time_bin)-2)){ 
#  temp_time_bsp$Slope[t] <- ((log(as.numeric(temp_time_bsp$Average_popSize[t+1])) - log(as.numeric(temp_time_bsp$Average_popSize[t])))/ (as.numeric(temp_time_bsp$Time_bin[t+1]) - as.numeric(temp_time_bsp$Time_bin[t])))
#}
bp <- boxplot(points_plot$Biome ~ as.numeric(points_plot$Time_bin), border="#EE760095", boxwex=.7, col="#EE760095", at=tick/1000, las=2, yaxt="n", frame=F, xaxt="n",  cex=0.5, xlim=c(0,50))
tick_labels <- c(0, rep(NA, 4), 5000, rep(NA, 4), 10000, rep(NA, 4), 15000, rep(NA, 4), 20000, rep(NA, 4), 25000, rep(NA, 4), 30000, rep(NA, 4), 35000, rep(NA, 4), 40000, rep(NA, 4), 45000, rep(NA, 4), 50000)
axis(side=1, at=seq(0, 50, by=1), labels=tick_labels, las=2, cex=0.5)
axis(side=4, at=seq(0, 1, by=0.2), cex=0.5)
mtext("Biome", side=2, line=3, cex.lab=2, col="#EE760095")
#text(tick/1000, rep(-0.001, length(tick)) ,labels=as.character(bp$n), cex=0.4)
lines(as.numeric(bp$names)/1000, bp$stats[3,], col="#EE7600", lwd=3)
#par(new=TRUE)
#plot(temp_time_bsp$Time_bin, temp_time_bsp$Slope, type="l", lwd=3, col="#595959",axes=F, ylab="", xlab="")
#axis(side=2, cex=0.5)
#mtext("Slope population size", side=4, line=3, cex.lab=1, col="red")
}
# plotting all
clim <- read.delim("~/Desktop/PhD/Thesis/1stChapter/clim_greenland", header=T, stringsAsFactors=F)
#### Extract the first 50000 years
clim_temp <- clim[clim$Years <= 50000,]
clim_50 <- clim_temp[c(0,seq(1,2000,2)),]
#### Delete the repeated years
X0<- matrix(ncol=37, nrow=21)
for (x in 1:21){
  X0[x,] <- bins_image
}
X1<- matrix(ncol=23, nrow=21)
for (x in 1:21){
  X1[x,] <- bins_image[1:23]
}
X2<- matrix(ncol=14, nrow=21)
for (x in 1:21){
  X2[x,] <- bins_image[24:37]
}
Y0<- matrix(ncol=37, nrow=21)
for (y in 1:37){
  Y0[,y] <- seq(-32,-12, by=1)
}
Y1<- matrix(ncol=23, nrow=21)
for (y in 1:23){
  Y1[,y] <- seq(-33,-13, by=1)
}
Y2<- matrix(ncol=14, nrow=21)
for (y in 1:14){
  Y2[,y] <- seq(-31,-11, by=1)
}
#poly.image(X1, Y1, z=image_climate_biome[,6:28])
#layout.show(ln)
#plot.new()
#### Plot image using all the values NO LAYOUT, add=T#####
colores_biomes <- terrain.colors(33, alpha=0.7)
ln<- layout(as.matrix(cbind(c(1,1,1), c(1,1,1),c(2,2,1))), heights=c(70,10,35),widths=(c(30,70,25)))
#plot.new()
par(mar=c(5,14,5,6))
plot(clim_50$Years, clim_50$Del_18O, type="l", xaxs="i", yaxs="i",xaxt="n", yaxt="n", ylab="",xlab="",frame.plot=F,xlim=c(0, 51000),ylim=c(-46, -10) )
mtext("Biome per species per time bin", side=3,line=2,cex=1.7)
#mtext("ranking relative to all the records", side=3,line=0)
axis(1, at=bins, labels=bins, las=2)
axis(4, at=(c(-46, -34)), labels=c(-46, -34), cex=0.7)
mtext("Temperature", side=4, line=2,at=-40)
#par(mar=c(0,0,0,0))
#poly.image(add=T,X1, Y1, z=image_climate[,6:28],col=paste(colores, 90, sep=""), ylab="", xlab="", xaxs="i", yaxs="i", xlim=c(-500, 22000), border="white")
#par(mar=c(0,0,0,0))
#poly.image(add=T,X2, Y2, z=image_climate[,29:42],col=paste(colores, 90, sep=""), ylab="", xlab="", xaxs="i", yaxs="i", xlim=c(22000, 50000), border="white")
#poly.image(add=T,X0, Y0, z=image_climate_biome[,6:42],col=colores_biomes, ylab="", xlab="", xaxs="i", yaxs="i", xlim=c(0, 50000), border="white")
poly.image(add=T,X0, Y0, z=image_climate_biome[,6:42],col=colores_biomes, ylab="", xlab="", xaxs="i", yaxs="i", xlim=c(0, 50000), border="white")


axis(1, at=bins,pos=-32.5, labels=NA)
axis(2, at=Y0[,1], labels=rownames(image_climate_biome), las=2)
#image(t(image_climate[,6:42]), col=paste(colores, 90, sep=""), axes=F,)
#axis(1, at=bins_image, labels=colnames(image_climate[,6:42]), las=2)
#text((x-1)/36, ((y-1)/20)-0.01, t(image_climate_n_biome[,6:42])[x,y], cex=0.5)
for (line in 1:dim(image_climate_n_biome[,6:42])[1]){
text(X0[1,], Y0[line], t(image_climate_n_biome[line,6:42]), cex=0.8)
}
par(mar=c(0,0,0,0))
#plot(1,1)
plot(1,1, col="white", frame.plot=F, axes=F, xlim=c(0,2), ylab=NA,xlab=NA)
legend.gradient(cbind(c(1.5,1.7,1.7,1.5), c(1.25,1.25,0.7,0.7)), col=colores_biomes,limits=c(min(na.omit(as.vector(image_climate_biome[,6:42]))), max(na.omit(as.vector(image_climate_biome[,6:42])))), title="")
segments(1.07,0.593,1.62,0.593)
segments(1.07,0.593,1.07,1.375)
segments(1.62,0.593,1.62,1.375)
segments(1.62,1.375,1.07,1.375)
#text(1.80,1.344, labels="Fast", srt=90)
#text(1.80,0.625, labels="Slow", srt=90)
















#########################################################################################################################
# Climate velocity
#########################################################################################################################





#### Plot image using all the colors
library(SDMTools)
colores <- colorRampPalette(c("#C6E2FF", "#63B8FF", "#1E90FF", "#7B68EE", "#D15FEE", "#FF69B4", "#EE3A8C", "#FF0000"))(6)
#colores <- colorRampPalette(c("#C6E2FF", "#FF0000"))(6)



clim <- read.delim("~/Desktop/PhD/Thesis/1stChapter/clim_greenland", header=T, stringsAsFactors=F)
#### Extract the first 50000 years
clim_temp <- clim[clim$Years <= 50000,]
clim_50 <- clim_temp[c(0,seq(1,2000,2)),]
#### Delete the repeated years
X0<- matrix(ncol=37, nrow=21)
for (x in 1:21){
  X0[x,] <- bins_image
}
X1<- matrix(ncol=23, nrow=21)
for (x in 1:21){
X1[x,] <- bins_image[1:23]
}
X2<- matrix(ncol=14, nrow=21)
for (x in 1:21){
  X2[x,] <- bins_image[24:37]
}
Y0<- matrix(ncol=37, nrow=21)
for (y in 1:37){
  Y0[,y] <- seq(-32,-12, by=1)
}
Y1<- matrix(ncol=23, nrow=21)
for (y in 1:23){
  Y1[,y] <- seq(-33,-13, by=1)
}
Y2<- matrix(ncol=14, nrow=21)
for (y in 1:14){
  Y2[,y] <- seq(-31,-11, by=1)
}
poly.image(X1, Y1, z=image_climate[,6:28])
layout.show(ln)
plot.new()
#### Plot image using all the values NO LAYOUT, add=T#####
ln<- layout(as.matrix(cbind(c(1,1,1), c(1,1,1),c(2,2,1))), heights=c(70,10,35),widths=(c(30,70,25)))
#plot.new()
par(mar=c(5,14,5,6))
plot(clim_50$Years, clim_50$Del_18O, type="l", xaxs="i", yaxs="i",xaxt="n", yaxt="n", ylab="",xlab="",frame.plot=F,xlim=c(0, 51000),ylim=c(-46, -10) )
mtext("Climate velocity per species per time bin", side=3,line=2,cex=1.7)
mtext("ranking relative to all the records", side=3,line=0)
axis(1, at=bins, labels=bins, las=2)
axis(4, at=(c(-46, -34)), labels=c(-46, -34), cex=0.7)
mtext("Temperature", side=4, line=2,at=-40)
#par(mar=c(0,0,0,0))
#poly.image(add=T,X1, Y1, z=image_climate[,6:28],col=paste(colores, 90, sep=""), ylab="", xlab="", xaxs="i", yaxs="i", xlim=c(-500, 22000), border="white")
#par(mar=c(0,0,0,0))
#poly.image(add=T,X2, Y2, z=image_climate[,29:42],col=paste(colores, 90, sep=""), ylab="", xlab="", xaxs="i", yaxs="i", xlim=c(22000, 50000), border="white")
poly.image(add=T,X0, Y0, z=image_climate[,6:42],col=paste(colores, 90, sep=""), ylab="", xlab="", xaxs="i", yaxs="i", xlim=c(0, 50000), border="white")
axis(1, at=bins,pos=-32.5, labels=NA)
#image(t(image_climate[,6:42]), col=paste(colores, 90, sep=""), axes=F,)
#axis(1, at=bins_image, labels=colnames(image_climate[,6:42]), las=2)
axis(2, at=Y0[,1], labels=rownames(image_climate), las=2)
par(mar=c(4.5,7,5,2))
#plot(1,1)
plot(1,1, col="white", frame.plot=F, axes=F, xlim=c(0,2), ylab=NA,xlab=NA)
legend.gradient(cbind(c(1.07,1.62,1.62,1.07), c(1.375,1.375,0.593,0.593)), cols=paste(colores, 99, sep=""),limits=c("",""), title="")
segments(1.07,0.593,1.62,0.593)
segments(1.07,0.593,1.07,1.375)
segments(1.62,0.593,1.62,1.375)
segments(1.62,1.375,1.07,1.375)
text(1.80,1.344, labels="Fast", srt=90)
text(1.80,0.625, labels="Slow", srt=90)
#text("Slowest", x=(1.44-0.55), y=0.7)
#text("Fastest", x=(1.44-0.55), y=1.2)

plot.new()
#########################################################################################################################

#########################################################################################################################

ln<- layout(as.matrix(cbind(c(1,1,1), c(1,1,1),c(2,2,1))), heights=c(70,10,35),widths=(c(30,70,25)))
#plot.new()
par(mar=c(5,14,5,6))
plot(clim_50$Years, clim_50$Del_18O, type="l", xaxs="i", yaxs="i",xaxt="n", yaxt="n", ylab="",xlab="",frame.plot=F,xlim=c(0, 51000),ylim=c(-46, -10) )
mtext("Climate velocity per species per time bin", side=3,line=2,cex=1.7)
mtext("ranking relative to single species records (per row)", side=3,line=0)
axis(1, at=bins, labels=bins, las=2)
axis(4, at=(c(-46, -34)), labels=c(-46, -34), cex=0.7)
mtext("Temperature", side=4, line=2,at=-40)
#par(mar=c(0,0,0,0))
#poly.image(add=T,X1, Y1, z=image_climate[,6:28],col=paste(colores, 90, sep=""), ylab="", xlab="", xaxs="i", yaxs="i", xlim=c(-500, 22000), border="white")
#par(mar=c(0,0,0,0))
#poly.image(add=T,X2, Y2, z=image_climate[,29:42],col=paste(colores, 90, sep=""), ylab="", xlab="", xaxs="i", yaxs="i", xlim=c(22000, 50000), border="white")
poly.image(add=T,X0, Y0, z=image_quartiles,col=paste(colores_quartiles, 90, sep=""), ylab="", xlab="", xaxs="i", yaxs="i", xlim=c(0, 50000), border="white")
axis(1, at=bins,pos=-32.5, labels=NA)
#image(t(image_climate[,6:42]), col=paste(colores, 90, sep=""), axes=F,)
#axis(1, at=bins_image, labels=colnames(image_climate[,6:42]), las=2)
axis(2, at=Y0[,1], labels=rownames(image_climate), las=2)
par(mar=c(4.5,7,5,2))
#plot(1,1)
plot(1,1, col="white", frame.plot=F, axes=F, xlim=c(0,2), ylab=NA,xlab=NA)
legend.gradient(cbind(c(1.07,1.62,1.62,1.07), c(1.375,1.375,0.593,0.593)), cols=paste(colores_quartiles, 99, sep=""),limits=c("",""), title="")
segments(1.07,0.593,1.62,0.593)
segments(1.07,0.593,1.07,1.375)
segments(1.62,0.593,1.62,1.375)
segments(1.62,1.375,1.07,1.375)
text(1.80,1.344, labels="Fast", srt=90)
text(1.80,0.625, labels="Slow", srt=90)

#########################################################################################################################

#########################################################################################################################



#### Plot image using all the values #####
ln<- layout(as.matrix(cbind(c(1,2), c(4,3))), heights=c(20,80),widths=(c(80,20)))
par(mar=c(1,13,4,2))
plot(clim_50$Years, clim_50$Del_18O, type="l", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylab="",frame.plot=F,ylim=c(-46, -34), main="Climate velocity per bin (color relative to ranking in the full database)")
axis(2, at=(c(-46, -34)), labels=c(-46, -34))
mtext("Temperature", side=2, line=2)
par(mar=c(5,13,0,2))
image(t(image_climate_rank[,6:42]), axes=F, col=paste(colores, 99, sep=""))
axis(1, at=(seq(0,1,by=1/36)), labels=colnames(image_climate_rank[,6:42]), las=2)
axis(2, at=(seq(0,1,by=1/20)), labels=rownames(image_climate_rank), las=2)
#axis(4, at=(seq(0,1,by=1/20)), labels=rowsum(image_climate_rank, ), las=2)
#for (x in 1:ncol(image_climate_rank)){
 # counter <- 0
 # for (y in 1:nrow(image_climate_rank)){
   # extemes<- c(max(t(image_climate_rank[,6:42])[,y], na.rm=T), min(t(image_climate_rank[,6:42])[,y], na.rm=T))
    #text((x-1)/36, ((y-1)/20)+0.01, t(image_climate_rank[,6:42])[x,y], cex=ifelse(is.na(sum(match(t(image_climate_rank[,6:42])[x,y],extemes))),0.3, 0.9))
    #text((x-1)/36, ((y-1)/20)-0.01, t(image_climate_n[,6:42])[x,y], cex=0.5)
    #counter <- counter + 1
    #print(sum(match(t(image_climate_rank)[x,y],extemes)))
  #}
#}
par(mar=c(0,0,0,0))
plot.new()
legend.gradient(cbind(c(0,0.5,0.5,0), c(1,1,0.5,0.5)), cols=paste(colores, 99, sep=""), limits=c("Slowest", "Fastest"), title="Ranking")


#### Plot image using quartiles boundaries
image_quartiles <- matrix(NA,ncol=37, nrow=21)
for (q in 1:dim(image_climate[,6:42])[1]){
  for (r in 1:dim(image_climate[,6:42])[2]){
    if (is.na(image_climate[q,r+5])){
      image_quartiles[q,r] <- NA
    }
    else if (image_climate[q,r+5] >= image_climate[q,1] & image_climate[q,r+5] < image_climate[q,2]){
      image_quartiles[q,r] <- 25
    }
    else if (image_climate[q,r+5] >= image_climate[q,2] & image_climate[q,r+5] < image_climate[q,3]){
      image_quartiles[q,r] <- 50
    }
    else if (image_climate[q,r+5] >= image_climate[q,3] & image_climate[q,r+5] < image_climate[q,4]){
      image_quartiles[q,r] <- 75
    }
    else if (image_climate[q,r+5] >= image_climate[q,4] & image_climate[q,r+5] < image_climate[q,5]){
      image_quartiles[q,r] <- 100
    }
}
}
colores_quartiles <- c("#BFEFFF", "#87CEEB", "#6CA6CD", "#4A708B")

ln<- layout(as.matrix(cbind(c(1,2), c(4,3))), heights=c(20,80),widths=(c(80,20)))
par(mar=c(1,13,4,2))
plot(clim_50$Years, clim_50$Del_18O, type="l", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylab="",frame.plot=F,ylim=c(-46, -34), main="Climate velocity per bin (color relative to quartiles per species)")
axis(2, at=(c(-46, -34)), labels=c(-46, -34))
mtext("Temperature", side=2, line=2)
par(mar=c(5,13,0,2))
image(t(image_quartiles), axes=F, col=paste(colores_quartiles, 99, sep=""))
axis(1, at=(seq(0,1,by=1/36)), labels=colnames(image_climate_rank[,6:42]), las=2)
axis(2, at=(seq(0,1,by=1/20)), labels=rownames(image_climate_rank), las=2)

#for (x in 1:ncol(image_quartiles)){
  #counter <- 0
  #for (y in 1:nrow(image_quartiles)){
   # extemes<- c(max(t(image_quartiles)[,y], na.rm=T), min(t(image_quartiles)[,y], na.rm=T))
   # text((x-1)/36, ((y-1)/20)+0.01, t(image_quartiles)[x,y], cex=0.6)
    #text((x-1)/36, ((y-1)/20)+0.01, t(image_quartiles)[x,y], cex=ifelse(is.na(sum(match(t(image_quartiles)[x,y],extemes))),0.3, 0.9))
    
    #text((x-1)/36, ((y-1)/20)+0.01, t(image_climate_rank[,6:42])[x,y], cex=ifelse(is.na(sum(match(t(image_climate_rank[,6:42])[x,y],extemes))),0.3, 0.9))
    #text((x-1)/36, ((y-1)/20)-0.01, t(image_climate_n[,6:42])[x,y], cex=0.5)
    #counter <- counter + 1
    #print(sum(match(t(image_climate_rank)[x,y],extemes)))
 # }
#}
par(mar=c(0,0,0,0))
plot.new()
legend.gradient(cbind(c(0,0.5,0.5,0), c(1,1,0.5,0.5)), cols=paste(colores_quartiles, 99, sep=""), limits=c("Slowest", "Fastest"), title="Ranking")
###################################
### Alex: Using generation times ###
###################################
setwd("/Users/afr/Desktop/Temporal_test/")
ABC_table_species <-  "species_constant.txt" 
ABC_table <- read.delim(ABC_table_species ,header=T, sep="\t", stringsAsFactors=FALSE, na.strings="NA")
clim_gen <- image_climate[match(ABC_table$Species, row.names(image_climate)),]
for (s in seq_along(rownames(clim_gen))){
  clim_gen[s,] <- (clim_gen[s,] *ABC_table$Gen_time[which(ABC_table$Species == row.names(clim_gen)[s])])
}

image_quartiles_gen <- matrix(NA,ncol=37, nrow=17)
for (q in 1:dim(clim_gen[,6:42])[1]){
  for (r in 1:dim(clim_gen[,6:42])[2]){
    if (is.na(clim_gen[q,r+5])){
      image_quartiles[q,r] <- NA
    }
    else if (clim_gen[q,r+5] >= clim_gen[q,1] & clim_gen[q,r+5] < clim_gen[q,2]){
      image_quartiles_gen[q,r] <- 25
    }
    else if (clim_gen[q,r+5] >= clim_gen[q,2] & clim_gen[q,r+5] < clim_gen[q,3]){
      image_quartiles_gen[q,r] <- 50
    }
    else if (clim_gen[q,r+5] >= clim_gen[q,3] & clim_gen[q,r+5] < clim_gen[q,4]){
      image_quartiles_gen[q,r] <- 75
    }
    else if (clim_gen[q,r+5] >= clim_gen[q,4] & clim_gen[q,r+5] < clim_gen[q,5]){
      image_quartiles_gen[q,r] <- 100
    }
  }
}
colores_quartiles <- c("#1E90FF","#C6E2FF", "#FFC0CB", "#FF0000")

ln<- layout(as.matrix(cbind(c(1,2), c(4,3))), heights=c(20,80),widths=(c(80,20)))
par(mar=c(1,13,4,2))
plot(clim_50$Years, clim_50$Del_18O, type="l", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylab="",frame.plot=F,ylim=c(-46, -34), main="Climate velocity per bin (color relative to quartiles per species/ gen time)")
axis(2, at=(c(-46, -34)), labels=c(-46, -34))
mtext("Temperature", side=2, line=2)
par(mar=c(5,13,0,2))
image(t(image_quartiles_gen), axes=F, col=paste(colores_quartiles, 99, sep=""))
#for (x in 1:ncol(image_quartiles)){
  #counter <- 0
  #for (y in 1:nrow(image_quartiles)){
   # extemes<- c(max(t(image_quartiles)[,y], na.rm=T), min(t(image_quartiles)[,y], na.rm=T))
    #text((x-1)/36, ((y-1)/16)+0.01, t(image_quartiles)[x,y], cex=0.6)
    axis(1, at=(seq(0,1,by=1/36)), labels=colnames(clim_gen[,6:42]), las=2)
    axis(2, at=(seq(0,1,by=1/16)), labels=rownames(clim_gen[,6:42]), las=2)
    
    
    #text((x-1)/36, ((y-1)/20)+0.01, t(image_quartiles)[x,y], cex=ifelse(is.na(sum(match(t(image_quartiles)[x,y],extemes))),0.3, 0.9))
    
    #text((x-1)/36, ((y-1)/20)+0.01, t(image_climate_rank[,6:42])[x,y], cex=ifelse(is.na(sum(match(t(image_climate_rank[,6:42])[x,y],extemes))),0.3, 0.9))
    #text((x-1)/36, ((y-1)/20)-0.01, t(image_climate_n[,6:42])[x,y], cex=0.5)
    #counter <- counter + 1
    #print(sum(match(t(image_climate_rank)[x,y],extemes)))
  #}
#}
par(mar=c(0,0,0,0))
plot.new()
legend.gradient(cbind(c(0,0.5,0.5,0), c(1,1,0.5,0.5)), cols=paste(colores_quartiles, 99, sep=""), limits=c("Slowest", "Fastest"), title="Ranking")

#################
ln<- layout(as.matrix(cbind(c(1,2), c(4,3))), heights=c(20,80),widths=(c(80,20)))
par(mar=c(1,13,4,2))
plot(clim_50$Years, clim_50$Del_18O, type="l", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylab="",frame.plot=F,ylim=c(-46, -34), main="Climate velocity per bin (color relative to full database species/ gen time)")
axis(2, at=(c(-46, -34)), labels=c(-46, -34))
mtext("Temperature", side=2, line=2)
par(mar=c(5,13,0,2))


image(t(clim_gen[,6:42]),col=paste(colores, 90, sep=""), axes=F)
axis(1, at=(seq(0,1,by=1/36)), labels=colnames(clim_gen[,6:42]), las=2)
axis(2, at=(seq(0,1,by=1/16)), labels=rownames(clim_gen[,6:42]), las=2)
par(mar=c(0,0,0,0))
plot.new()
legend.gradient(cbind(c(0,0.5,0.5,0), c(1,1,0.5,0.5)), cols=paste(colores, 90, sep=""), limits=c("Slowest", "Fastest"), title="Ranking")
########## using log for the values
ln<- layout(as.matrix(cbind(c(1,2), c(4,3))), heights=c(20,80),widths=(c(80,20)))
par(mar=c(1,13,4,2))
plot(clim_50$Years, clim_50$Del_18O, type="l", xaxs="i", yaxs="i", xaxt="n", yaxt="n", ylab="",frame.plot=F,ylim=c(-46, -34), main="Climate velocity per bin (color relative to full database species/ gen time)")
axis(2, at=(c(-46, -34)), labels=c(-46, -34))
mtext("Temperature", side=2, line=2)
par(mar=c(5,13,0,2))
image(t(image_quartiles_gen_2),col=colores_quartiles, axes=F)
axis(1, at=(seq(0,1,by=1/36)), labels=colnames(clim_gen[,6:42]), las=2)
axis(2, at=(seq(0,1,by=1/16)), labels=rownames(clim_gen[,6:42]), las=2)
par(mar=c(0,0,0,0))
plot.new()
legend.gradient(cbind(c(0,0.5,0.5,0), c(1,1,0.5,0.5)), cols=colores_quartiles, limits=c("Slowest", "Fastest"), title="Ranking")
bp_all_points_gen <- boxplot(as.vector(clim_gen[,6:42]), plot=F)
image_quartiles_gen_2 <- matrix(NA,ncol=37, nrow=17)
for (q in 1:dim(clim_gen[,6:42])[1]){
  for (r in 1:dim(clim_gen[,6:42])[2]){
    if (is.na(clim_gen[q,r+5])){
      image_quartiles[q,r] <- NA
    }
    else if (clim_gen[q,r+5] >= bp_all_points_gen$stats[1] & clim_gen[q,r+5] < bp_all_points_gen$stats[2]){
      image_quartiles_gen_2[q,r] <- 25
    }
    else if (clim_gen[q,r+5] >= bp_all_points_gen$stats[2] & clim_gen[q,r+5] < bp_all_points_gen$stats[3]){
      image_quartiles_gen_2[q,r] <- 50
    }
    else if (clim_gen[q,r+5] >= bp_all_points_gen$stats[3] & clim_gen[q,r+5] < bp_all_points_gen$stats[4]){
      image_quartiles_gen_2[q,r] <- 75
    }
    else if (clim_gen[q,r+5] >= bp_all_points_gen$stats[4] & clim_gen[q,r+5] < bp_all_points_gen$stats[5]){
      image_quartiles_gen_2[q,r] <- 100
    }
  }
}
plot(as.vector(log(clim_gen[,6:42])))
plot(as.vector(clim_gen[,6:42]))
hist(as.vector(clim_gen[,6:42]))
hist(log(as.vector(clim_gen[,6:42]))
hist(log(as.vector(clim_gen[,6:42]))+2.5)

###################################
### Alex: dating error and bins ###
###################################
bins_5to0 <- c(seq(50000, 22000, by=-2000), seq(21000, 0, by=-1000))
bins_0to5 <- c(seq(0, 21000,by=1000),seq(22000,52000, by=2000))
for (s in seq_along(Single_sp)){
  temp_DB_climate <- Full_DB_LL[which(Full_DB_LL$Species==Single_sp[s]),]
  temp_DB_sort <- temp_DB_climate[order(temp_DB_climate$Mean_Age),]
  plot(x=temp_DB_sort$Mean_Age, y=seq(1,length(temp_DB_sort$Median_Age),by=1), xlim=c(0,55000), pch=19, frame.plot=F,cex=0.2, yaxt="n",xaxt="n", ylab="", xlab="", xaxs="i", yaxs="i", main=Single_sp[s])
  h<-hist(temp_DB_sort$Mean_Age, breaks=bins_0to5, freq=T, col="#63B8FF80", border="white", add=T)
  axis(1,at=bins_0to5, las=2)
  axis(2,at=c(0, max(h$counts)))
  segments(x0=(temp_DB_sort$Mean_Age)-(temp_DB_sort$Cal_Sigma),y0=seq(1,length(temp_DB_sort$Median_Age),by=1), x1=(temp_DB_sort$Mean_Age)+(temp_DB_sort$Cal_Sigma))
  points(x=temp_DB_sort$Median_Age, y=seq(1,length(temp_DB_sort$Median_Age),by=1),col="red", pch=19, cex=0.3)
}
  



time <-  c(seq(50000, 22000, by=-2000), seq(21000, -1000, by=-1000))
height <- rep(1,times=length(time))
frame <- cbind(height, time)
barplot(rep(2,times=38), height=height, space=0, col="#63B8FF60", border="white")


##############################
### Alex: Correlate to BSP ###
##############################
x <- stats::runif(12); y <- stats::rnorm(12)
i <- order(x, y); x <- x[i]; y <- y[i]
plot(x,y, main = "arrows(.) and segments(.)")
## draw arrows from point to point :
s <- seq(length(x)-1)  # one shorter than data
arrows(x[s], y[s], x[s+1], y[s+1], col = 1:3)
s <- s[-length(s)]
segments(x[s], y[s], x[s+2], y[s+2], col = "pink")

### Plot velocity profile per time bin
#setwd("/Users/afr/Desktop/Evolution_ppt")
#bsp_raw <- na.omit(read.delim(paste(species_sp, "_BSP_50_data.txt", sep=""), header=T, stringsAsFactors=F, sep="\t"))
#bsp_raw <- na.omit(read.delim(paste(species_vector[s], "_BSP_50_data.txt", sep=""), header=T, stringsAsFactors=F, sep="\t"))
#plot(bsp_raw$Time, log(bsp_raw$Median), main=species_vector[s], ylim=c(4,16), lwd=3, type="l",frame.plot=F, lab=c(10,10,5), ylab="Population size", xlab="", col="#008B45", cex=0.5, xaxt="n")
#lines(bsp_raw$Time, log(bsp_raw$Mean), lty=2)
#par(new=TRUE)

length(bp$n)



temp_time_bsp <- as.data.frame(matrix("numeric", nrow=37, ncol=3), stringsAsFactors=F)
colnames(temp_time_bsp) <- c("Time_bin", "Average_popSize", "Slope")
temp_time_bsp$Time_bin <- c(seq(50000, 22000, by=-2000), seq(21000, 0, by=-1000))
# Estimate the average pop size
for (t in 1:(length(temp_time_bsp$Time_bin)-1)){
  temp_time_bsp$Average_popSize[t] <- mean.default(bsp_raw[which(bsp_raw$Time <= temp_time_bsp$Time_bin[t] & bsp_raw$Time >= temp_time_bsp$Time_bin[t+1]),5])
}
# Estimate the slope of the average pop size
for (t in 1:(length(temp_time_bsp$Time_bin)-2)){ 
  temp_time_bsp$Slope[t] <- ((log(as.numeric(temp_time_bsp$Average_popSize[t+1])) - log(as.numeric(temp_time_bsp$Average_popSize[t])))/ (as.numeric(temp_time_bsp$Time_bin[t+1]) - as.numeric(temp_time_bsp$Time_bin[t])))
}
bp <- boxplot(points_plot$Velocity ~ as.numeric(points_plot$Time_bin), border="#EE760095", boxwex=.7, col="#EE760095", at=tick/1000, las=2, yaxt="n", frame=F, xaxt="n",  cex=0.5, xlim=c(0,50), main=species_sp)
tick_labels <- c(0, rep(NA, 4), 5000, rep(NA, 4), 10000, rep(NA, 4), 15000, rep(NA, 4), 20000, rep(NA, 4), 25000, rep(NA, 4), 30000, rep(NA, 4), 35000, rep(NA, 4), 40000, rep(NA, 4), 45000, rep(NA, 4), 50000)
axis(side=1, at=seq(0, 50, by=1), labels=tick_labels, las=2, cex=0.5)
axis(side=4, at=seq(0, 1, by=0.2), cex=0.5)
mtext("Velocity", side=2, line=3, cex.lab=1, col="#EE760095")
#text(tick/1000, rep(-0.001, length(tick)) ,labels=as.character(bp$n), cex=0.4)
lines(as.numeric(bp$names)/1000, bp$stats[3,], col="#EE7600", lwd=3)
par(new=TRUE)
plot(temp_time_bsp$Time_bin, temp_time_bsp$Slope, type="l", lwd=3, col="#008B45",axes=F, ylab="", xlab="")
axis(side=2, cex=0.5)
mtext("Slope population size", side=4, line=3, cex.lab=1, col="red")




############## 
#### Correlation
temp_time_cli <- as.data.frame(matrix("numeric", nrow=length(bp$names), ncol=3), stringsAsFactors=F)
colnames(temp_time_cli) <- c("Time_bin", "Median_vel_clim", "slope_bsp")
temp_time_cli$Time_bin <- as.numeric(bp$names)
temp_time_cli$Median_vel_clim <- as.numeric(bp$stats[3,])
temp_time_cli$slope_bsp <- na.omit(temp_time_bsp$Slope[match(temp_time_bsp$Time_bin,temp_time_cli$Time_bin)])
plot(log(temp_time_cli$Median_vel_clim), log(as.numeric(temp_time_cli$slope_bsp)+1), main=paste(species_sp, "log", sep=" "))
cor.test(temp_time_cli$Median_vel_clim, sqrt(as.numeric(temp_time_cli$slope_bsp)))
temp_ccf <- ccf(temp_time_cli$Median_vel_clim, as.numeric(temp_time_cli$slope_bsp))
barplot(as.vector(temp_ccf$acf), names.arg=as.vector(temp_ccf$lag), col="#00B2EE95", border="white", ylim=c(min(as.vector(temp_ccf$acf))-0.2, max(as.vector(temp_ccf$acf)+0.2)))
temp_ccf <- ccf(temp_time_cli$Median_vel_clim, as.numeric(temp_time_cli$slope_bsp))
lag.plot(temp_time_cli$Median_vel_clim, layout=c(3,3))
#plot(temp_time_cli$Median_vel_clim, as.numeric(temp_time_cli$slope_bsp), main=species_vector[s])
}

temp_ccf <- ccf(temp_time_cli$Median_vel_clim, as.numeric(temp_time_cli$slope_bsp))
str(mp_ccf)
lag.plot(temp_time_cli$Median_vel_clim, 9)
lag.plot(as.numeric(mp_time_cli$slope_bsp))

