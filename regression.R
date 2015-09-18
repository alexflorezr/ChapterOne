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
for (s in seq_along(Single_sp$Species)){
  # Database to work with, contains all the species and data
  temp_DB_climate <- Full_DB_LL[which(Full_DB_LL$Species==Single_sp$Species[s]),]
  # Empty data frame to save the data for every species
  temp_points <- as.data.frame(matrix(nrow=nrow(temp_DB_climate), ncol=14))
  colnames(temp_points) <- c("Longitude", "Latitude", "Mean_date_record", "Time_bin", "Layer", "Rec_type", "Cell", "Cli_vel_tmp", "Cli_vel_prc", "Cli_vel_tnp", "Biome_kpp","Biome_bio")
  #for(record in seq_along(temp_DB_climate$Mean_Age)){
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
   assign(paste(Single_sp$Sp[s], "points", sep="_"), temp_df)
}

######## normal plot points
temp_points_NA <- na.omit(temp_points)
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
# correlation pop size and climate velocity
setwd("/Users/afr/Desktop/Regression/BSP/")
temp_BSP <- na.omit(read.delim(paste(tolower(Single_sp$Sp[s]), "_BSP_50_data.txt", sep=""), header=T, stringsAsFactors=F, sep="\t"))
plot(temp_BSP$Time, log(temp_BSP$Median))
