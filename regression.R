rm(list=ls())
library(raster)
library(rgdal)
# Upload the climate velocity files
setwd("/Users/afr/Desktop/Regression/Cli_vel/")
cli_vel_Tmp <- stack("climate_change_velocity_perYear_tmp_rescaled_truncated_0_000005.grd")
cli_vel_Prc <- stack("climate_change_velocity_perYear_prec_rescaled_truncated_0_000005.grd")
cli_vel_TnP <- stack("climate_change_velocity_perYear_tmp_prec_rescaled_truncated_0_000005.grd")
# Upload the biomes files
setwd("/Users/afr/Desktop/Regression/Biomes/Koppen_full_cor/")
biome_kpp <- stack(sort(dir(),decreasing=T))
setwd("/Users/afr/Desktop/Regression/Biomes/Bio4_CO2/")
biome_bio <- stack(sort(dir(),decreasing=T))
# objects describing the variables
## Vector for the time bins from the present to the past
bins_0to50 <- c(seq(0,21000, by=1000), seq(22000, 50000, by=2000))
bins_50to0 <- c(seq(50000, 22000, by=-2000), seq(21000, 0, by=-1000))
layers <- seq(25, 61, by=1)
## Vector for the time bins middle point from the present to the past
bins_mid_point <- c(seq(-500, 21500, by=1000), seq(23000, 49000, by=2000))
# main loop for extracting values for every species
## Define the species to use in the loop
setwd("/Users/afr/Desktop/Regression")
Single_sp <- read.delim("Single_sp_regression", header=T, stringsAsFactors=F)
#########################################################################################################################
#########################################################################################################################
##### Extract the climatic values for every fossil record in the database
## TODO ## Convert this part of the script in a function
#########################################################################################################################
#########################################################################################################################
species_set_clim <- as.data.frame(matrix(nrow=0, ncol=12))
species_set_BSP <- as.data.frame(matrix(nrow=0, ncol=6))
for (s in seq_along(Single_sp$Species)){
  # Database to work with, contains all the species and data
  temp_DB_climate <- Full_DB_LL[which(Full_DB_LL$Species==Single_sp$Species[s]),]
  # Empty data frame to save the data for every species
  temp_points <- as.data.frame(matrix(nrow=nrow(temp_DB_climate), ncol=12))
  colnames(temp_points) <- c("Longitude", "Latitude", "Mean_date_record", "Time_bin", "Layer", "Rec_type", "Cell", "Cli_vel_tmp", "Cli_vel_prc", "Cli_vel_tnp", "Biome_kpp","Biome_bio")
  counter1 <- 1
  for (bin in seq_along(bins_50to0)[-37]){
    records <- which(temp_DB_climate$Mean_Age <= bins_50to0[bin] & temp_DB_climate$Mean_Age > bins_50to0[bin+1])
    points <- as.matrix(cbind.data.frame(as.numeric(temp_DB_climate$Longitude[records]), as.numeric(temp_DB_climate$Latitude[records])))
    colnames(points) <- c("longitude", "Latitude")
      if(length(records > 0)){
        temp_points[(counter1:(counter1 + length(records)-1)),1:6] <- cbind.data.frame(as.numeric(temp_DB_climate$Longitude[records]), as.numeric(temp_DB_climate$Latitude[records]),as.numeric(temp_DB_climate$Mean_Age[records]), rep(bins_50to0[bin], times=length(records)), rep(layers[bin], times=length(records)),ifelse(nchar(temp_DB_climate$Sequence[records]) > 1, "Seq", "Fossil") )
        temp_points[(counter1:(counter1 + length(records)-1)),which(colnames(temp_points) == c("Cell", "Cli_vel_tmp"))]<- extract(cli_vel_Tmp, layer=layers[bin], nl=1, y=points, cellnumbers=T)
        temp_points[(counter1:(counter1 + length(records)-1)),which(colnames(temp_points) == "Cli_vel_prc")]<- extract(cli_vel_Prc, layer=layers[bin], nl=1, y=points)
        temp_points[(counter1:(counter1 + length(records)-1)),which(colnames(temp_points) == "Cli_vel_tnp")]<- extract(cli_vel_TnP, layer=layers[bin], nl=1, y=points)
        temp_points[(counter1:(counter1 + length(records)-1)),which(colnames(temp_points) == "Biome_kpp")]<- extract(biome_kpp, layer=(layers[bin])+1, nl=1, y=points)
        temp_points[(counter1:(counter1 + length(records)-1)),which(colnames(temp_points) == "Biome_bio")]<- extract(biome_bio, layer=(layers[bin])+1, nl=1, y=points)
        counter1 <- counter1+length(records)
      }
  } 
  assign(paste(Single_sp$Sp[s], "points", sep="_"), temp_points)
  if (Species_set_clim ==TRUE){
    species_set_clim <- rbind(species_set_clim, get(paste(Single_sp$Sp[s], "points", sep="_")))
  }
#########################################################################################################################
#########################################################################################################################
##### Estimate the values of BSP slope for all the species in the species vector
## TODO ## Convert this part of the script in a function
#########################################################################################################################
#########################################################################################################################
  setwd("/Users/afr/Desktop/Regression/BSP/")
  temp_BSP <- read.delim(paste(tolower(Single_sp$Sp[s]), "_BSP_50_data.txt", sep=""), header=T, stringsAsFactors=F, sep="\t")
  temp_BSP_bin <- as.data.frame(matrix(nrow=length(temp_BSP$Median), ncol=6), stringsAsFactors=F)
  colnames(temp_BSP_bin) <- c("Time_bin", "Time_pop", "Pop_size_mean", "Pop_size_median", "Slope_pop_mean", "Slope_pop_median")
  # Estimate the average pop size
  counter2 <- 1
  for (pop_time in seq_along(bins_50to0)[-37]){
    pop_per_bin <- which(temp_BSP$Time <= bins_50to0[pop_time] & temp_BSP$Time > bins_50to0[pop_time+1])
    temp_BSP_bin[counter2:(counter2+length(pop_per_bin)-1), c(-5,-6)] <- cbind(rep(bins_50to0[pop_time], times=length(pop_per_bin)),temp_BSP[pop_per_bin,c(1,2,3)])
    counter2 <- counter2 + length(pop_per_bin)
  }
  for (time_pop in seq_along(temp_BSP_bin$Pop_size_mean)[-length(temp_BSP_bin$Pop_size_mean)]){
    temp_BSP_bin[time_pop, 5] <- (log(temp_BSP_bin$Pop_size_mean[time_pop+1]) - log(temp_BSP_bin$Pop_size_mean[time_pop]))/ (temp_BSP_bin$Time_pop[time_pop+1] - temp_BSP_bin$Time_pop[time_pop])
    temp_BSP_bin[time_pop, 6] <- (log(temp_BSP_bin$Pop_size_median[time_pop+1]) - log(temp_BSP_bin$Pop_size_median[time_pop]))/ (temp_BSP_bin$Time_pop[time_pop+1] - temp_BSP_bin$Time_pop[time_pop])
  }
  assign(paste(Single_sp$Sp[s], "BSP", sep="_"), temp_BSP_bin)
  if (Species_set_BSP ==TRUE){
  species_set_BSP <- rbind(species_set_BSP, get(paste(Single_sp$Sp[s], "BSP", sep="_")))
}
}
mean_clim_bsp <- as.data.frame(matrix(nrow=length(bins_50to0), ncol=5))
colnames(mean_clim_bsp) <- c("Time_bin", "Mean_clim_tmp", "Mean_clim_prc","Mean_clim_tnp", "Mean_slope_BSP")
mean_clim_bsp$Time_bin <- bins_50to0
for (bin in seq_along(bins_50to0)){
  mean_clim_bsp$Mean_clim_tmp[which(mean_clim_bsp$Time_bin == bins_0to50[bin])] <- mean.default(na.omit(species_set_clim[which(species_set_clim$Time_bin == bins_0to50[bin]),8]))
  mean_clim_bsp$Mean_clim_prc[which(mean_clim_bsp$Time_bin == bins_0to50[bin])] <- mean.default(na.omit(species_set_clim[which(species_set_clim$Time_bin == bins_0to50[bin]),9]))
  mean_clim_bsp$Mean_clim_tnp[which(mean_clim_bsp$Time_bin == bins_0to50[bin])] <- mean.default(na.omit(species_set_clim[which(species_set_clim$Time_bin == bins_0to50[bin]),10]))
  mean_clim_bsp$Mean_slope_BSP[which(mean_clim_bsp$Time_bin == bins_0to50[bin])] <- mean.default(na.omit(species_set_BSP[which(species_set_BSP$Time_bin == bins_0to50[bin]),5])) 
}


#### some plots 
bp_BSP <- boxplot(log(species_set_BSP$Slope_pop_mean) ~ species_set_BSP$Time_bin, col="green")
par(new=T)
bp_tmp <- boxplot(log(species_set_clim$Cli_vel_tmp) ~ species_set_clim$Time_bin, col="red")
bp_prc <- boxplot(log(species_set_clim$Cli_vel_prc) ~ species_set_clim$Time_bin, col="blue", add=T)
bp_tnp <- boxplot(log(species_set_clim$Cli_vel_tnp) ~ species_set_clim$Time_bin, col="grey", add=T)

plot(bp_tmp$stats[3,],bp_BSP$stats[3,], col="#CD000098", pch=19)
plot(bp_tmp$stats[3,],bp_BSP$stats[3,], col="#CD000098", pch=19
abline(lm(bp_BSP$stats[3,] ~ bp_tmp$stats[3,]), col="#CD0000")
cor_temp <- cor(bp_tmp$stats[3,], bp_BSP$stats[3,])
text(x=-2,y=-10, round(cor_temp, 3), col="#CD0000")
points(bp_prc$stats[3,], bp_BSP$stats[3,],col="#5CACEE98", pch=19)
plot(bp_prc$stats[3,], bp_BSP$stats[3,],col="#5CACEE98", pch=19)
abline(lm(bp_BSP$stats[3,] ~ bp_prc$stats[3,] ), col="#5CACEE")
cor_prc <- cor(bp_prc$stats[3,], bp_BSP$stats[3,])
text(x=-2,y=-9.5, round(cor_prc, 3), col="#5CACEE")
points(bp_tnp$stats[3,], bp_BSP$stats[3,],col="#8B8989", pch=19)
abline(lm(bp_BSP$stats[3,] ~ bp_tnp$stats[3,]), col="#8B8989")
cor_tnp <-cor(bp_tnp$stats[3,], bp_BSP$stats[3,])
cor.test(bp_prc$stats[3,], bp_BSP$stats[3,])
cor.test(bp_tmp$stats[3,], residuals(object=lm(bp_BSP$stats[3,] ~ bp_prc$stats[3,])))

cov(bp_prc$stats[3,], bp_BSP$stats[3,])

###### Residuals
lag2.plot(bp_tmp$stats[3,], residuals(object=lm(bp_BSP$stats[3,] ~ bp_prc$stats[3,])), 9)
lag2.plot(bp_prc$stats[3,], bp_BSP$stats[3,], 9)
############## stepwise regresion 
data_fit <- as.data.frame(cbind(bp_BSP$stats[3,],bp_tmp$stats[3,], bp_prc$stats[3,]))

fit_set <- lm(data_fit[,1] ~ ., data=data_fit)
step_set <- step(fit_set, direction="both", )

install.packages("MASS")
library(MASS)



ccf(bp_tmp$stats[3,],bp_BSP$stats[3,])
ccf(bp_prc$stats[3,],bp_BSP$stats[3,])
lag2.plot(bp_tmp$stats[3,],bp_BSP$stats[3,], 9)
lag2.plot(bp_prc$stats[3,],bp_BSP$stats[3,], 9)
######## normal plot points
temp_points_NA <- na.omit(temp_points[,1:12])
par(new=T)
plot(temp_points_NA$Time_bin, temp_points_NA$Cli_vel_tmp, col="#CD000098", pch=19, cex=1.2)
abline(lm(temp_points$Cli_vel_tmp ~ temp_points$Time_bin), col="#EE2C2C98", lwd=2)
points(temp_points_NA$Time_bin, temp_points_NA$Cli_vel_prc, col="#5CACEE98", pch=19, cex=1.2)
abline(lm(temp_points$Cli_vel_prc ~ temp_points$Time_bin), col="#5CACEE98", lwd=2)
points(temp_points_NA$Time_bin, temp_points_NA$Cli_vel_tnp, col="#8B898998", pch=19, cex=1.2)
abline(lm(temp_points$Cli_vel_tnp ~ temp_points$Time_bin), col="#8B898998", lwd=2)
######## normal boxplot
par(mar=c(4,4,4,4))
bp_tmp <- boxplot(temp_points$Cli_vel_tmp ~ as.numeric(temp_points$Time_bin), plot=F)
bp_tmp <- boxplot(temp_points$Cli_vel_tmp ~ as.numeric(temp_points$Time_bin), col="#CD000098", las=2, frame=F,  cex=0.5, xlim=c(0,50), at=as.numeric(bp_tmp$names)/1000, xaxt="n",  cex=0.5,)
abline(lm(bp_tmp$stats[3,] ~ as.numeric(bp$names)), col="#EE2C2C98", lwd=2)
bp_prc <- boxplot(temp_points$Cli_vel_prc ~ as.numeric(temp_points$Time_bin), plot=F)
bp_prc <- boxplot(temp_points$Cli_vel_prc ~ as.numeric(temp_points$Time_bin), col="#5CACEE98", las=2, frame=F,  cex=0.5, xlim=c(0,50), at=as.numeric(bp_prc$names)/1000, xaxt="n", yaxt="n", cex=0.5, add=T)
abline(lm(bp_prc$stats[3,] ~ as.numeric(bp$names)), col="#5CACEE98", lwd=2)
bp_tnp <- boxplot(temp_points$Cli_vel_tnp ~ as.numeric(temp_points$Time_bin), plot=F)
bp_tnp <- boxplot(temp_points$Cli_vel_tnp ~ as.numeric(temp_points$Time_bin), col="#8B898998", las=2, frame=F,  cex=0.5, xlim=c(0,50), at=as.numeric(bp_tnp$names)/1000, xaxt="n", yaxt="n", cex=0.5, add=T)
abline(lm(bp_tnp$stats[3,] ~ as.numeric(bp$names)), col="#8B898998", lwd=2)
tick_labels <- c(0, rep(NA, 4), 5000, rep(NA, 4), 10000, rep(NA, 4), 15000, rep(NA, 4), 20000, rep(NA, 4), 25000, rep(NA, 4), 30000, rep(NA, 4), 35000, rep(NA, 4), 40000, rep(NA, 4), 45000, rep(NA, 4), 50000)
axis(side=1, at=seq(0, 50, by=1), labels=tick_labels, las=2, cex=0.5)
######## correlation among climate velocity estimates
par(mar=c(6,6,8,6))
plot(temp_points$Cli_vel_tmp, temp_points$Cli_vel_prc, pch=19, cex=1.2, col="#EE762190", frame=F, xlim=c(-0.05,6), ylim=c(-0.1,6), xaxs="i",yaxs="i", xlab=NA, ylab=NA)
mtext("No correlation amongst climate velocity of different climatic variables ", side=3, line=5, cex=1.5)
mtext("Var 1", side=1, line=3)
mtext("Var 2", side=2, line=3)
abline(lm(temp_points$Cli_vel_prc ~ temp_points$Cli_vel_tmp), col="#EE7621", lwd=2)
points(temp_points$Cli_vel_tmp, temp_points$Cli_vel_tnp, pch=19, cex=1.2, col="#8B450090")
abline(lm(temp_points$Cli_vel_tnp ~ temp_points$Cli_vel_tmp), col="#8B4500", lwd=2)
points(temp_points$Cli_vel_prc, temp_points$Cli_vel_tnp, pch=19, cex=1.2, col="#00640090")
abline(lm(temp_points$Cli_vel_tnp ~ temp_points$Cli_vel_prc), col="#006400", lwd=2)
cor_all_cli_vel <- cor(na.omit(temp_points[,c(8,9,10)]))
text(round(cor_all_cli_vel[which(row.names(cor_all_cli_vel) == "Cli_vel_tmp"), which(colnames(cor_all_cli_vel) == "Cli_vel_prc")], 3), x=5,y=1.4, cex=0.5)
segments(x0=4.5,y0=1.4,x1=4.8,y1=1.4, col="#EE7621", lwd=2)
text(round(cor_all_cli_vel[which(row.names(cor_all_cli_vel) == "Cli_vel_tmp"), which(colnames(cor_all_cli_vel) == "Cli_vel_tnp")], 3), x=5,y=1.2, cex=0.5)
segments(x0=4.5,y0=1.2,x1=4.8,y1=1.2, col="#8B4500", lwd=2)
text(round(cor_all_cli_vel[which(row.names(cor_all_cli_vel) == "Cli_vel_prc"), which(colnames(cor_all_cli_vel) == "Cli_vel_tnp")], 3), x=5,y=1.0, cex=0.5)
segments(x0=4.5,y0=1.0,x1=4.8,y1=1.0, col="#006400", lwd=2)
######## correlation among logarithmic climate velocity estimates
par(new=T)
plot(log(temp_points$Cli_vel_tmp+0.001), log(temp_points$Cli_vel_prc+0.001), pch=17, cex=1.2, col="#EE762190", frame=F, xaxs="i",yaxs="i", lwd=0.5, xlim=c(-8,2), ylim=c(-8,2), xaxt="n", yaxt="n", xlab=NA, ylab=NA)
abline(lm(log(temp_points$Cli_vel_prc+0.01) ~ log(temp_points$Cli_vel_tmp+0.01)), col="#EE7621", lwd=2, lty=2)
axis(side=3, at=seq(-8, 2, by=1), cex=0.5)
axis(side=4, at=seq(-8, 2, by=1), cex=0.5)
mtext("Log(Var 1) + 0.01", side=3, line=3)
mtext("Log(Var 2) + 0.01", side=4, line=3)
points(log(temp_points$Cli_vel_tmp+0.001), log(temp_points$Cli_vel_tnp+0.001), pch=17, cex=1.2, col="#8B450090",frame=F, xaxs="i",yaxs="i",lwd=0.5)
abline(lm(log(temp_points$Cli_vel_tnp+0.01) ~ log(temp_points$Cli_vel_tmp+0.01)), col="#8B4500", lwd=2, lty=2)
points(log(temp_points$Cli_vel_prc+0.001), log(temp_points$Cli_vel_tnp+0.001), pch=17, cex=1.2, col="#00640090", frame=F, xaxs="i",yaxs="i",lwd=0.5)
abline(lm(log(temp_points$Cli_vel_tnp+0.01) ~ log(temp_points$Cli_vel_prc+0.01)), col="#006400", lwd=2, lty=2)
cor_all_cli_vel_log <- cor(na.omit(log(temp_points[,c(8,9,10)]+0.01)))
text(round(cor_all_cli_vel_log[which(row.names(cor_all_cli_vel_log) == "Cli_vel_tmp"), which(colnames(cor_all_cli_vel_log) == "Cli_vel_prc")], 3), x=1,y=0.2, cex=0.5)
segments(x0=0.2,y0=0.2,x1=0.8,y1=0.2, col="#EE7621", lwd=2, lty=2)
text(round(cor_all_cli_vel_log[which(row.names(cor_all_cli_vel_log) == "Cli_vel_tmp"), which(colnames(cor_all_cli_vel_log) == "Cli_vel_tnp")], 3), x=1,y=0, cex=0.5)
segments(x0=0.2,y0=0,x1=0.8,y1=0, col="#8B4500", lwd=2, lty=2)
text(round(cor_all_cli_vel_log[which(row.names(cor_all_cli_vel_log) == "Cli_vel_prc"), which(colnames(cor_all_cli_vel_log) == "Cli_vel_tnp")], 3), x=1,y=-0.2, cex=0.5)
segments(x0=0.2,y0=-0.2,x1=0.8,y1=-0.2, col="#006400", lwd=2, lty=2)
# box plots for the biomes 
boxplot(as.numeric(temp_points$Biome_kpp) ~ temp_points$Time_bin, border="#EE762190", col="#8B450090", frame=F, ylim=c(0,33), xaxs="i")
par(new=T)
Biome_bio_thousand <- temp_points$Biome_bio
Biome_bio_thousand[which(Biome_bio_thousand == -1000)] <- 0
boxplot(Biome_bio_thousand ~ temp_points$Time_bin, border="#00640090", col="#00640090", yaxt="n", frame=F, ylim=c(0,33), xaxs="i")
axis(side=4)
########### Correlation BSP and climate velocity
setwd("/Users/afr/Desktop/Regression/BSP/")
temp_BSP <- read.delim(paste(tolower(Single_sp$Sp[s]), "_BSP_50_data.txt", sep=""), header=T, stringsAsFactors=F, sep="\t")
#bsp_raw <- na.omit(read.delim(paste(species_vector[s], "_BSP_50_data.txt", sep=""), header=T, stringsAsFactors=F, sep="\t"))
plot(temp_BSP$Time, log(temp_BSP$Median), main=Single_sp[s], ylim=c(4,16), lwd=3, type="l",frame.plot=F, lab=c(10,10,5), ylab="Population size", xlab="", col="#008B45", cex=0.5, xaxt="n")
points(temp_BSP$Time, log(temp_BSP$Mean), lty=2, col="black")
#plot.new()
#par(new=TRUE)
temp_BSP_bin <- as.data.frame(matrix(nrow=length(temp_BSP$Median), ncol=6), stringsAsFactors=F)
colnames(temp_BSP_bin) <- c("Time_bin", "Time_pop", "Pop_size_mean", "Pop_size_median", "Slope_pop_mean", "Slope_pop_median")
# Estimate the average pop size
counter2 <- 1
for (pop_time in seq_along(bins_50to0)[-37]){
  pop_per_bin <- which(temp_BSP$Time <= bins_50to0[pop_time] & temp_BSP$Time > bins_50to0[pop_time+1])
  temp_BSP_bin[counter2:(counter2+length(pop_per_bin)-1), c(-5,-6)] <- cbind(rep(bins_50to0[pop_time], times=length(pop_per_bin)),temp_BSP[pop_per_bin,c(1,2,3)])
  counter2 <- counter2 + length(pop_per_bin)
}
for (time_pop in seq_along(temp_BSP_bin$Pop_size_mean)[-length(temp_BSP_bin$Pop_size_mean)]){
temp_BSP_bin[time_pop, 5] <- (log(temp_BSP_bin$Pop_size_mean[time_pop+1]) - log(temp_BSP_bin$Pop_size_mean[time_pop]))/ (temp_BSP_bin$Time_pop[time_pop+1] - temp_BSP_bin$Time_pop[time_pop])
temp_BSP_bin[time_pop, 6] <- (log(temp_BSP_bin$Pop_size_median[time_pop+1]) - log(temp_BSP_bin$Pop_size_median[time_pop]))/ (temp_BSP_bin$Time_pop[time_pop+1] - temp_BSP_bin$Time_pop[time_pop])
}

boxplot(log(temp_BSP_bin$Pop_size_median) ~ temp_BSP_bin$Time_bin)
points(temp_BSP_bin$Time_bin, log(temp_BSP_bin$Pop_size_median))
plot(temp_BSP_bin$Time_bin, log(temp_BSP_bin$Pop_size_median), add=T)
plot()
####### normal plot BSP slope
plot(temp_BSP_bin$Time_pop, temp_BSP_bin$Slope_pop_mean, pch=19, col="#00640090")
points(temp_BSP_bin$Time_pop, temp_BSP_bin$Slope_pop_median, pch=17, col="#8B450090")
###### LOG plot BSP slope
plot(temp_BSP_bin$Time_pop, log(temp_BSP_bin$Slope_pop_mean), pch=19, col="#00640090", xlim=c(0,51000),frame=F)
points(temp_BSP_bin$Time_pop, log(temp_BSP_bin$Slope_pop_median), pch=19, col="#8B450090")
#points(temp_points_NA$Time_bin, temp_points_NA$Cli_vel_tmp -8, col="#CD000098", pch=17, cex=1.2)
#bp_tmp <- boxplot(temp_points$Cli_vel_tmp ~ as.numeric(temp_points$Time_bin), col="#CD000098", las=2, frame=F,  cex=0.5, xlim=c(0,50), at=as.numeric(bp_tmp$names)/1000, xaxt="n",  cex=0.5, add=T)
box_slope_median <- boxplot(temp_BSP_bin$Slope_pop_median ~ temp_BSP_bin$Time_bin, las=2)
box_slope_mean <- boxplot(temp_BSP_bin$Slope_pop_mean ~ temp_BSP_bin$Time_bin, add=T, col="red")
par(new=T)
plot(as.numeric(bp_tmp$names), bp_tmp$stats[3,], xlim=c(0,50000 ))

plot.new()
plot(as.numeric(box_slope_mean$names),box_slope_mean$stats[3,], col="red",type="l")
points(as.numeric(box_slope_median$names),box_slope_median$stats[3,], col="blue",type="l")
par(new=T)
# plot the median values of climate velocity for every bin
par(mar=c(6,6,6,6))
plot(as.numeric(bp_tnp$names), bp_tnp$stats[3,], pch=17, las=2, xlim=c(0,51000), col="#8B898998", ylab=NA, yaxt="n", xaxt="n", frame=F, xlab=NA, xaxs="i", yaxs="i", ylim=c(0,0.65),las=2)
lines(as.numeric(bp_tnp$names), bp_tnp$stats[3,], col="#8B898998", lwd=3)
axis(1, las=2)
mtext("Time (years before present)",side=1, line=3)
axis(2)
mtext("Climate velocity (kilometers per year)",side=2, line=3)
points(as.numeric(bp_tmp$names), bp_tmp$stats[3,], pch=17, col="#CD000098")
lines(as.numeric(bp_tmp$names), bp_tmp$stats[3,], col="#CD000098", lwd=3)
points(as.numeric(bp_prc$names), bp_prc$stats[3,], pch=17, col="#5CACEE98")
lines(as.numeric(bp_prc$names), bp_prc$stats[3,], col="#5CACEE98", lwd=3)
# add the slope values for population size
par(new=T)
plot(as.numeric(box_slope_mean$names), box_slope_mean$stats[3,], pch=19, col="#2E8B5798", xlim=c(0,51000), xaxs="i", yaxs="i", xaxt="n", yaxt="n", frame=F, ylim=c(-0.00005,0.0009), ylab=NA, xlab=NA)
lines(as.numeric(box_slope_mean$names), box_slope_mean$stats[3,], col="#2E8B5798", lwd=3, lty=2)
points(as.numeric(box_slope_median$names), box_slope_median$stats[3,], pch=19, col="#228B2298")
lines(as.numeric(box_slope_median$names), box_slope_median$stats[3,], col="#228B2298", lwd=3, lty=2)
axis(4)
mtext("Change in population size (slope)",side=4, line=3)
abline(h=0)


### correlation test
cor_bins <- as.vector(as.numeric(bp_tmp$names))
cor_matrix <- as.data.frame(matrix(ncol=6, nrow=length(cor_bins)))
colnames(cor_matrix) <- c("Time_bin", "Slope_BSP_mean","Slope_BSP_median", "Cli_vel_tmp", "Cli_vel_prc", "Cli_vel_tnp")
cor_matrix[,1] <- cor_bins
for (cor_bin in seq_along(cor_matrix$Time_bin)){
  cor_matrix[,2] <- na.omit(box_slope_mean$stats[3,match(as.numeric(box_slope_mean$names),cor_matrix$Time_bin)])
  cor_matrix[,3] <- na.omit(box_slope_median$stats[3,match(as.numeric(box_slope_median$names),cor_matrix$Time_bin)])
  cor_matrix[,4] <- na.omit(bp_tmp$stats[3,match(as.numeric(bp_tmp$names),cor_matrix$Time_bin)])
  cor_matrix[,5] <- na.omit(bp_prc$stats[3,match(as.numeric(bp_prc$names),cor_matrix$Time_bin)])
  cor_matrix[,6] <- na.omit(bp_tnp$stats[3,match(as.numeric(bp_tnp$names),cor_matrix$Time_bin)])
}
cor(log(cor_matrix[2:length(cor_matrix[,2]),2]+0.5), log(cor_matrix[1:length(cor_matrix[,4])-1,6]+0.5))
# normal plot Slopes vs climate velocity
plot(cor_matrix[,2], cor_matrix[,2])
points(cor_matrix[,3], cor_matrix[,3], pch=19, col="#CD000098")
plot(cor_matrix[,4], cor_matrix[,2], pch=19, col="#CD000098" )
points(cor_matrix[,5], cor_matrix[,2], pch=19,col="#5CACEE98")
points(cor_matrix[,6], cor_matrix[,2], pch=19, col="#8B898998")
points(cor_matrix[,4], cor_matrix[,3], pch=17, col="#CD000098" )
points(cor_matrix[,5], cor_matrix[,3], pch=17,col="#5CACEE98")
points(cor_matrix[,6], cor_matrix[,3], pch=17, col="#8B898998")
# LOG plot Slopes vs climate velocity
plot(sqrt(cor_matrix[,4]), cor_matrix[,2], pch=19, col="#CD000098" )
abline(lm(cor_matrix[,2] ~ log(cor_matrix[,4])), col="#CD000098", lwd=3)
points(sqrt(cor_matrix[,5]), cor_matrix[,2], pch=19,col="#5CACEE98")
abline(lm(cor_matrix[,2] ~ log(cor_matrix[,5])),col="#5CACEE98", lwd=3)
points(sqrt(cor_matrix[,6]), cor_matrix[,2], pch=19, col="#8B898998")
abline(lm(cor_matrix[,2] ~ log(cor_matrix[,6])), col="#8B898998", lwd=3)
 points(cor_matrix[,4], cor_matrix[,3], pch=17, col="#CD000098" )
points(cor_matrix[,5], cor_matrix[,3], pch=17,col="#5CACEE98")
points(cor_matrix[,6], cor_matrix[,3], pch=17, col="#8B898998")
lag2.plot(cor_matrix[,4], cor_matrix[,2], max.lag=9) 

lag2.plot(cor_matrix[,4], cor_matrix[,2], max.lag=9, smooth=F)
lag2.plot(log(cor_matrix[,4]), cor_matrix[,2], max.lag=9, smooth=F)

lag2.plot(cor_matrix[,5], cor_matrix[,2], max.lag=9, smooth=F)
lag2.plot(log(cor_matrix[,5]), cor_matrix[,2], max.lag=9, smooth=F)

lag2.plot(cor_matrix[,6], cor_matrix[,2], max.lag=9, smooth=F)
lag2.plot(log(cor_matrix[,6]), cor_matrix[,2], max.lag=9, smooth=F)

lag4 <- lag(as.ts(cor_matrix[,3]), 4)
start <- attributes(lag4)$tsp[1]
end <- attributes(lag4)$tsp[2]
lag2.plot(log(cor_matrix[,5]), cor_matrix[,3], 9)
plot(log(cor_matrix[1:21,5]),cor_matrix[-start:end,3])
abline(lm(cor_matrix[start:end,3] ~ log(cor_matrix[1:21,5])))
lag2.plot
str(lag)
lag <- lag(cor_matrix[,3], k=5)
par(mfrow=c(3,2))
for (lags in 1:5){
plot(lag(as.ts(cor_matrix[,3]), lags), as.ts(cor_matrix[,4]),xy.labels = FALSE )
}
cor(lag(as.ts(cor_matrix[,3]), 5),as.ts(cor_matrix[,4]))
lag2.plot(cor_matrix[,3],cor_matrix[,4], 9)
?lag
par(mfrow=c(2,2))
plot(log(abs(cor_matrix[,2])))

plot(log10((cor_matrix[,2]+1)))
plot(abs(cor_matrix[,2]))
kk <- ccf(log(cor_matrix[,4]), cor_matrix[,2])
plot(cor_matrix[,1], cor_matrix[,2])
lag2.plot(log(cor_matrix[,6]), cor_matrix[,2],9)
axis(4)
plot(cor_matrix[,1], cor_matrix[,2])
points(cor_matrix[,1], cor_matrix[,3], col="red")
plot(cor_matrix[,1], cor_matrix[,4], col="green", pch=17)
points(cor_matrix[,1], cor_matrix[,5], col="black", pch=17)
points(cor_matrix[,1], cor_matrix[,6], col="yellow", pch=17)

plot(log(cor_matrix[,4]),cor_matrix[,3])
lm_slope_cli <- lm(cor_matrix[,3] ~ log(cor_matrix[,4]))
abline(lm_slope_cli)
resid <- residuals(lm_slope_cli)
plot(resid)
lag.plot(cor_matrix[,2], 9)
# Estimate the slope of the average pop size
for (t in 1:(length(temp_BSP_bin$Time_bin)-2)){ 
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
points(temp_time_bsp$Time_bin, temp_time_bsp$Slope, col="#595959", ylab="", xlab="")

####### all species correlation
plot(St_points$Time_bin, St_points$Cli_vel_tmp, , pch=19, col="#CD000098", ylim=c(0,2))
points(Bs_points$Time_bin, Bs_points$Cli_vel_tmp, , pch=19, col="#CD000098")
points(Mp_points$Time_bin, Mp_points$Cli_vel_tmp, , pch=19, col="#CD000098")
points(Om_points$Time_bin, Om_points$Cli_vel_tmp, , pch=19, col="#CD000098")
points(Rt_points$Time_bin, Rt_points$Cli_vel_tmp, , pch=19, col="#CD000098")



