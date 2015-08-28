install.packages("fields")
library(fields)
dat<- data(RCMexample)
data
set.panel( 1,2)
par(pty="s")
# plot with grid modified
poly.image( RCMexample$x, RCMexample$y, RCMexample$z[,,1])
x<-RCMexample$x
z<-RCMexample$z[,,1]
# use midpoints of z
poly.image( RCMexample$x, RCMexample$y, RCMexample$z[,,1],midpoint=TRUE)

# images are very similar. 
set.panel()
# Regridding of x and y
l1<- poly.image.regrid( RCMexample$x)
l2<- poly.image.regrid( RCMexample$y)

# test that this works
i<- 1:10
plot( l1[i,i], l2[i,i])
points( RCMexample$x[i,i], RCMexample$y[i,i],col="red")


for(h in seq_along(bins)){
hist((temp_points$Biome[temp_points$Time_bin == bins[h]]), breaks=seq(0,33, by=1), main=bins[h], ylim=c(0,33))
}

pie_biome <- matrix(ncol=37, nrow=1)
colnames(pie_biome) <- c(seq(0, 21000, by=1000), seq(22000, 50000, by=2000))
for(h in seq_along(bins)){
  temp_hist <- hist((temp_points$Biome[temp_points$Time_bin == bins[h]]), breaks=seq(0,33, by=1), main=bins[h], ylim=c(0,33), plot=F)
  if (sum(temp_hist$counts != 0) > 0){
     temp_pie <- pie(temp_hist$counts, main=paste(bins[h]," n=",sum(temp_hist$counts), sep=""), col=colores_biomes, border="white", labels=NA, lty=0.5)
  }
  else {
    plot(1, type="n", axes=F, xlab="", ylab="")
  }
}
###########
par(mfrow=c(21,37), mar=c(0,0,0,0))
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
  for(h in seq_along(bins)){
    temp_hist <- hist((temp_points$Biome[temp_points$Time_bin == bins[h]]), breaks=seq(0,33, by=1), main=bins[h], ylim=c(0,33), plot=F)
    if (sum(temp_hist$counts != 0) > 0){
      temp_pie <- pie(temp_hist$counts, , col=colores_biomes, border="#242424", labels=NA,lwd=1)
      text(1,1,sum(temp_hist$counts), cex=1.2)
    }
    else {
      plot(1, type="n", axes=F, xlab="", ylab="")
    }
  }
}

main=sum(temp_hist$counts)
pie()