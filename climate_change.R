#####################
# Climate velocity ## WORKING, MAKE PLOT FOR ALL THE 4 SPECIES
#####################
library(raster)
# Upload the climate velocity files
setwd("/Users/afr/Desktop/CHECK/Ditte_files/Velocity/paleo_cc_velocity/")
velocity_map_year <- stack("climate_change_velocity_perYear_tmp_prec_rescaled_everyKyrs_truncated_0_00005.grd")
# Define the directory where the files
setwd("/Users/afr/Desktop/Evolution_ppt/Corr_result/")
image_climate <- matrix(NA,ncol=37, nrow=21)
colnames(image_climate) <- c( seq(0, 21000, by=1000), seq(22000, 50000, by=2000))
rownames(image_climate) <- Single_sp
image_climate_rank <- image_climate
image_climate_n <- image_climate
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
  bp <- boxplot(points_plot$Velocity ~ as.numeric(points_plot$Time_bin), plot=F)
  tick <- as.vector(as.numeric(bp$names))
  # Check point
  sort(as.numeric(points_plot$Time_bin))
  # Save the velocity values for every species
### setwd("/Users/afr/Desktop/kk_temp/")
### write.table(x=temp_points, file=paste(tolower(Single_sp[s]), "_vel_bin.txt", sep=""), sep="\t",row.names=F)
### species_sp <- paste(strsplit(strsplit(tolower(Single_sp[s]), split = "_")[[1]][1], split = "")[[1]][1], strsplit(strsplit(tolower(Single_sp[s]), split = "_")[[1]][2], split = "")[[1]][1], sep="")
bp <- boxplot(points_plot$Velocity ~ as.numeric(points_plot$Time_bin), main=Single_sp[s], border="#EE760095", boxwex=.7, col="#EE760095", at=tick/1000, las=2, yaxt="n", frame=F, xaxt="n",  cex=0.5, xlim=c(0,50))
tick_labels <- c(0, rep(NA, 4), 5000, rep(NA, 4), 10000, rep(NA, 4), 15000, rep(NA, 4), 20000, rep(NA, 4), 25000, rep(NA, 4), 30000, rep(NA, 4), 35000, rep(NA, 4), 40000, rep(NA, 4), 45000, rep(NA, 4), 50000)
 axis(side=1, at=seq(0, 50, by=1), labels=tick_labels, las=2, cex=0.5)
 axis(side=4, at=seq(0, 1, by=0.2), cex=0.5)
lines(as.numeric(bp$names)/1000, bp$stats[1,], col="red", lwd=2)
lines(as.numeric(bp$names)/1000, bp$stats[2,], col="pink", lwd=2)
 lines(as.numeric(bp$names)/1000, bp$stats[3,], col="grey", lwd=2)
lines(as.numeric(bp$names)/1000, bp$stats[4,], col="lightblue", lwd=2)
lines(as.numeric(bp$names)/1000, bp$stats[5,], col="green", lwd=2)
  image_climate[s,as.vector(na.omit(match(bp$names, colnames(image_climate))))] <- bp$stats[3,]
  image_climate_rank[s,as.vector(na.omit(match(bp$names, colnames(image_climate))))] <- rank(bp$stats[3,])
  image_climate_n[s,as.vector(na.omit(match(bp$names, colnames(image_climate))))] <- bp$n                                                                                    
}

library(SDMTools)
colores <- colorRampPalette(c("#C6E2FF", "#63B8FF", "#1E90FF", "#7B68EE", "#D15FEE", "#FF69B4", "#EE3A8C", "#FF0000"))(6)
colores <- colorRampPalette(c("#C6E2FF", "#FF0000"))(6)


par(mar=c(5,16,4,4))
image(t(image_climate), col=paste(colores, 90, sep=""), axes=F)
image(t(image_climate_rank), axes=F, col=paste(colores, 99, sep=""), xlim=c(0,1.3))
axis(1, at=(seq(0,1,by=1/36)), labels=colnames(image_climate_rank), las=2)
axis(2, at=(seq(0,1,by=1/20)), labels=rownames(image_climate_rank), las=2)
#axis(4, at=(seq(0,1,by=1/20)), labels=rowsum(image_climate_rank, ), las=2)
for (x in 1:ncol(image_climate_rank)){
  counter <- 0
  for (y in 1:nrow(image_climate_rank)){
    extemes<- c(max(t(image_climate_rank)[,y], na.rm=T), min(t(image_climate_rank)[,y], na.rm=T))
    text((x-1)/36, ((y-1)/20)+0.01, t(image_climate_rank)[x,y], cex=ifelse(is.na(sum(match(t(image_climate_rank)[x,y],extemes))),0.3, 0.9))
    text((x-1)/36, ((y-1)/20)-0.01, t(image_climate_n)[x,y], cex=0.5)
    #counter <- counter + 1
    #print(sum(match(t(image_climate_rank)[x,y],extemes)))
  }
}
legend.gradient(cbind(c(1.05,1.2,1.2,1.05), c(0.95,0.95,0.7,0.7)), cols=paste(colores, 99, sep=""), limits=c("Slowest", "Fastest"), title="Ranking")

min(as.vector(na.omit(image_climate)))
which(image_climate == min(image_climate), arr.ind = TRUE)
min(image_climate,na.rm=T)
max(image_climate,na.rm=T)

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

