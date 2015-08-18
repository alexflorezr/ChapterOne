#####################
# Climate velocity ## WORKING, MAKE PLOT FOR ALL THE 4 SPECIES
#####################
library(raster)
setwd("/Users/afr/Desktop/Evolution_ppt/Paleo_cc_velocity/")
velocity_map_year <- stack("climate_change_velocity_perYear_tmp_prec_rescaled_everyKyrs_truncated_0_00005.grd")
# Define the directory where the files
setwd("/Users/afr/Desktop/Evolution_ppt/Corr_result/")
# Read the table with the species information 
##### Species_info <- read.delim("species_info.txt", header = T, sep="\t", stringsAsFactors=F)
#Create a grid for the graphs
##### par(mar=c(2,2,1,2), mfrow=c(2,5))
#### Loop throughout the species to extract the velocity values per time bin for every record####

species_vector <- c("Mammuthus_primigenius", "Bison_bison", "Coelodonta_antiquitatis", "Ovibos_moschatus")
#species_vector <- c("Mammuthus_primigenius")
for (s in seq_along(species_vector)){
  temp_db <- Full_DB_LL[Full_DB_LL$Species == species_vector[s],]
  temp_points <- as.data.frame(matrix(nrow=nrow(temp_db), ncol=7))
  colnames(temp_points) <- c("Longitude", "Latitude", "Time_sample", "Time_bin", "Layer", "cell","Velocity")
  matrix_time <- matrix("numeric", nrow=38, ncol=2 )
  matrix_time[,1] <- c(seq(50000, 22000, by=-2000), seq(21000, -1000, by=-1000))
  matrix_time[,2] <- seq(25, 62, by=1)
  for(i in seq_along(temp_db$Median_Age)){ # Takes longitude and latitude for data older than 5000 years and send it into a point matrix
    for (j in 1:(dim(matrix_time)[1]-1)){
      if(temp_db$Median_Age[i]<= as.numeric(matrix_time[j]) & temp_db$Median_Age[i]> as.numeric(matrix_time[j+1])){
        temp_points[i,c(1,2,3,4,5)] <- c(temp_db$Longitude[i], temp_db$Latitude[i],temp_db$Median_Age[i], matrix_time[j,1], matrix_time[j,2] )
      }
    }
  }
  for(k in seq_along(temp_points$Longitude)){
    temp_points[k,c(6,7)]<- extract(velocity_map_year, layer=as.numeric(temp_points[k,5]), nl=1, y=matrix(as.numeric(temp_points[k,c(1,2)]), nrow=1, ncol=2), cellnumbers=T)
  }
  points_plot <- temp_points[!is.na(temp_points$Velocity),]
  bp <- boxplot(points_plot$Velocity ~ as.numeric(points_plot$Time_bin), plot=F)
  tick <- as.vector(as.numeric(bp$names))
  # Check point
  sort(as.numeric(points_plot$Time_bin))
  # Save the velocity values for every species
  setwd("/Users/afr/Desktop/Evolution_ppt/Corr_result/")
  write.table(x=temp_points, file=paste(tolower(species_vector[s]), "_vel_bin.txt", sep=""), sep="\t",row.names=F)
  species_sp <- paste(strsplit(strsplit(tolower(species_vector[s]), split = "_")[[1]][1], split = "")[[1]][1], strsplit(strsplit(tolower(species_vector[s]), split = "_")[[1]][2], split = "")[[1]][1], sep="")
  ### Plot velocity profile per time bin
  setwd("/Users/afr/Desktop/Evolution_ppt")
  bsp_raw <- na.omit(read.delim(paste(species_sp, "_BSP_50_data.txt", sep=""), header=T, stringsAsFactors=F, sep="\t"))
  #bsp_raw <- na.omit(read.delim(paste(species_vector[s], "_BSP_50_data.txt", sep=""), header=T, stringsAsFactors=F, sep="\t"))
  plot(bsp_raw$Time, log(bsp_raw$Median), main=species_vector[s], ylim=c(4,16), lwd=3, type="l",frame.plot=F, lab=c(10,10,5), ylab="Population size", xlab="", col="#008B45", cex=0.5, xaxt="n")
  #lines(bsp_raw$Time, log(bsp_raw$Mean), lty=2)
  par(new=TRUE)
  bp <- boxplot(points_plot$Velocity ~ as.numeric(points_plot$Time_bin), border="#EE760095", boxwex=.7, col="#EE760095", at=tick/1000, las=2, yaxt="n", frame=F, xaxt="n",  cex=0.5, xlim=c(0,50))
  tick_labels <- c(0, rep(NA, 4), 5000, rep(NA, 4), 10000, rep(NA, 4), 15000, rep(NA, 4), 20000, rep(NA, 4), 25000, rep(NA, 4), 30000, rep(NA, 4), 35000, rep(NA, 4), 40000, rep(NA, 4), 45000, rep(NA, 4), 50000)
  axis(side=1, at=seq(0, 50, by=1), labels=tick_labels, las=2, cex=0.5)
  axis(side=4, at=seq(0, 1, by=0.2), cex=0.5)
  #mtext("Velocity", side=4, line=3, cex.lab=1, col="red")
  #text(tick/1000, rep(-0.001, length(tick)) ,labels=as.character(bp$n), cex=0.4)
  lines(as.numeric(bp$names)/1000, bp$stats[3,], col="#EE7600", lwd=2)
#}

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
bp <- boxplot(points_plot$Velocity ~ as.numeric(points_plot$Time_bin), border="#EE760095", boxwex=.7, col="#EE760095", at=tick/1000, las=2, yaxt="n", frame=F, xaxt="n",  cex=0.5, xlim=c(0,50), main=species_vector[s])
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
plot(log(temp_time_cli$Median_vel_clim), log(as.numeric(temp_time_cli$slope_bsp)+1), main=paste(species_vector[s], "log", sep=" "))
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

