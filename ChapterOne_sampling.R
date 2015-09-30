
rm(list=ls())
Full_DB <- read.delim(file.choose(), header=T, stringsAsFactors=F)
Full_DB$Longitude <- as.numeric(Full_DB$Longitude)
Full_DB$Latitude <- as.numeric(Full_DB$Latitude)
full_vector <- as.vector(NULL)
for (i in seq_along(Full_DB$Latitude)){
  if (is.na(Full_DB$Longitude[i]) || is.na(Full_DB$Latitude[i])){
    full_vector <- c(full_vector,i)
  }
}
# remove the records without longitude and or latitude
Full_DB_LL <- Full_DB[-full_vector,]
# remove the records older than 50000 years
DATABASE <- Full_DB_LL[-which(Full_DB_LL$Median_Age > 50000),]

setwd("~/Desktop/Ch1_Figures/Figure_one/")
par(mfrow=c(4,4), mar=c(2,3,3,1), yaxs="i",xaxs="i",lwd=0.5)
Sp_LP_H <- read.delim("sp_LP-H.txt", stringsAsFactors=F, header=T)
#Sampling_aDNA <- function(Sp_LP_H, database, ind_color){
require(ape)
require(pegas)
require(adegenet)
plot.new()
Sp_bef_aft <- as.data.frame(matrix(nrow=dim(Sp_LP_H)[1], ncol=3))
colnames(Sp_bef_aft) <- c("Species", "Bef", "Aft")
Event_name <- "Late Pleistoce - Holocene"
setwd("~/Desktop/Ch1_Figures/Figure_one/Same_length_seqs/")
for (sp in seq_along(Sp_LP_H$Species)){
  #par(mar=c(5,5,5,5),yaxs="i",xaxs="i",lwd=0.5)
  if(is.element(Sp_LP_H$Species[sp], unique(DATABASE$Species))){
## STARTS: histogram for the fossil record and sequences ##
    temp_sp_rec_db <- Full_DB[which(Full_DB$Species == Sp_LP_H$Species[sp] & Full_DB$Mean_Age <= 50000),]
    temp_hist_rec <- hist(temp_sp_rec_db$Mean_Age, breaks=seq(0, 50000, 2000), plot=F)
    temp_hist_rec$counts[which(temp_hist_rec$counts == 0)] <- NA
    temp_max <- max(na.omit(temp_hist_rec$counts))
    temp_genus <- strsplit(Sp_LP_H$Species[sp], split="_")[[1]][1]
    temp_specific <- strsplit(Sp_LP_H$Species[sp], split="_")[[1]][2]
    main_name <- paste(temp_genus, temp_specific, sep=" ")
    temp_sp_seq_db <- temp_sp_rec_db[nchar(temp_sp_rec_db$Sequence) > 1,]
    temp_hist_seq <- hist(temp_sp_seq_db$Mean_Age,breaks=seq(0, 71000, 2000), plot=F )
    temp_hist_seq$counts[which(temp_hist_seq$counts == 0)] <- NA
    sp_color <- "#00BFFF"
    #Events <- as.numeric(Sp_LP_H[sp,c(2,3)])
    Events <- 11700
    temp_sp <-  paste(strsplit(temp_genus, split="")[[1]][1], strsplit(temp_specific, split="")[[1]][1], sep="")
    if( hist_rec_seq == TRUE){
      hist(temp_sp_rec_db$Mean_Age, breaks=seq(0, 50000, 2000), main=NULL, xaxt="n", yaxt="n", ylab=NULL, xlab=NULL, col="#00BFFF50", border=NA)
      axis(side=1)
      axis(side=2)
      SS_name_vector <- c("n","nuc div", "hap", "seg", "theta")
      for (name in seq_along(SS_name_vector)){
        text(paste(SS_name_vector[name], " =", sep=""), x=7000, y=(temp_max/2)+sum(rep(temp_max/12, times=name)), col="#363636", adj=c(1,0), cex=0.7)
      }
      mtext(side=3, main_name, line=1)
      hist(temp_sp_seq_db$Mean_Age,breaks=seq(0, 50000, 2000), add=T, col=paste(sp_color, 90, sep=""), border=NA)
      abline(v=Events, lwd=5, col=c("#FF7F00", "#EEC591"))
      #labels=as.character(temp_hist_rec$counts)
      #if (ind_color == T) {sp_color <= paste(Sp_LP_H$Color_color[sp], 90, sep="")}
      #text(x=temp_hist_seq$mids, y=0.5, labels=as.character(temp_hist_seq$counts), col="white", cex=0.7)
    }
    for (event in seq_along(Events)){  
      temp_db_seqs_bef <- temp_sp_seq_db[temp_sp_seq_db$Mean_Age > Events[event], ]
      temp_db_seqs_aft <- temp_sp_seq_db[temp_sp_seq_db$Mean_Age <= Events[event], ]
      temp_fasta <- read.FASTA(paste(tolower(temp_sp), "_fasta_new.fasta", sep=""))
      temp_db_Temp_Net <- temp_fasta[NULL]
      temp_vec <- c("temp_db_seqs_bef","temp_db_seqs_aft")
      start <- 1
      end <- dim(temp_db_seqs_bef)[1]
        for (bef_aft in seq_along(temp_vec)){
          temp_count_seq <- dim(get(temp_vec[bef_aft]))[1]
          temp_acc <- get(temp_vec[bef_aft])$Accession_GB
          temp_age <- get(temp_vec[bef_aft])$Mean_Age
          temp_err <- get(temp_vec[bef_aft])$Cal_Sigma
          temp_rep <- get(temp_vec[bef_aft])$Repeat
          Temp_seq_name_fasta <- paste(temp_sp,temp_acc, temp_rep, temp_age,temp_err, sep="_") 
          temp_seqs  <- temp_fasta[as.vector(na.exclude(match(Temp_seq_name_fasta, labels(temp_fasta))))]
          temp_db_Temp_Net <- c(temp_db_Temp_Net,  temp_seqs)
          names(temp_db_Temp_Net)[start:end] <- paste(names(temp_db_Temp_Net), "$", bef_aft, sep="")
          start <- end+1
          end <- dim(temp_db_seqs_bef)[1]+dim(temp_db_seqs_aft)[1]
          ## pop gen SS ##
            if (length(temp_seqs) > 0){
              temp_nuc_div <- nuc.div(temp_seqs)
              temp_num_hap <- haplotype(temp_seqs)
              temp_segSites <- seg.sites(temp_seqs)
              temp_theta  <- theta.s(length(temp_segSites), length(temp_seqs))
            } else{
              temp_nuc_div <- -999
              temp_num_hap <- as.matrix(-999, -999)
              temp_segSites <- -999
              temp_theta  <- -999
            }
          if (temp_vec[bef_aft] == "temp_db_seqs_bef"){
            temp_seqs_bef <- temp_seqs
            SS_bef <- c(temp_count_seq, temp_nuc_div, round(dim(temp_num_hap)[1]), round(length(temp_segSites)), temp_theta)
          }
          if (temp_vec[bef_aft] == "temp_db_seqs_aft"){
            temp_seqs_aft <- temp_seqs
            SS_aft <- c(temp_count_seq, temp_nuc_div, round(dim(temp_num_hap)[1]), round(length(temp_segSites)), temp_theta)
          }
        }
      write.dna(temp_db_Temp_Net, file=paste(temp_sp, "_Temp_Net.fasta", sep=""), format="fasta")
      assign(paste(temp_sp, "_Temp_Net", sep=""), TempNet(file=paste(temp_sp, "_Temp_Net.fasta", sep=""), save=T))
    }
      }
      if( hist_rec_seq == TRUE){
        for (ss_bef in seq_along(SS_bef)){
        text(round(SS_bef[ss_bef], digits=3), x=Events[event]+500, y=(temp_max/2)+sum(rep(temp_max/12, times=ss_bef)), col="#363636", adj=c(0,0), cex=0.7)
        }
        for (ss_aft in seq_along(SS_aft)){
          text(round(SS_aft[ss_aft], digits=3), x=Events[event]-500, y=(temp_max/2)+sum(rep(temp_max/12, times=ss_aft)), col="#363636", adj=c(1,0), cex=0.7)
        } 
      }
      if (subsampling == TRUE){
        nucdiv_sub <- as.numeric(vector(length=10000))
          if (length(temp_seqs_bef) > length(temp_seqs_aft)){ 
            for (sub in 1:10000){
              sub_seqs <- sample((1:length(temp_seqs_bef)), size=length(temp_seqs_aft), replace=T)
              nucdiv_sub[sub] <- nuc.div(temp_seqs_bef[sub_seqs])
            }
          }
          if (length(temp_seqs_bef) < length(temp_seqs_aft)){
              for (sub in 1:10000){
                sub_seqs <- sample((1:length(temp_seqs_aft)), size=length(temp_seqs_bef), replace=T)
                nucdiv_sub[sub] <- nuc.div(temp_seqs_aft[sub_seqs])
              }
          }
        if (length(temp_seqs_bef) == length(temp_seqs_aft)){
          temp_seqs_bef <- temp_seqs_bef[-sample((1:length(temp_seqs_bef)), 1)]
          for (sub in 1:10000){
            sub_seqs <- sample((1:length(temp_seqs_aft)), size=length(temp_seqs_bef), replace=T)
            nucdiv_sub[sub] <- nuc.div(temp_seqs_aft[sub_seqs])
          }
        }  
        }
      if (colnames(Sp_LP_H)[1+event] == "Fast"){
        col_event<- "#FFD39B"
      }
      if (colnames(Sp_LP_H)[1+event] == "Slow"){
        col_event<- "#FFD39B"
      }
      if (plot_sub == TRUE){
        if (event == 1){
          temp_hist_sub <- hist(nucdiv_sub, plot=F)
          temp_max_y_sub <- max(temp_hist_sub$counts)
          
          temp_max_x_sub <- max(temp_hist_sub$mids)
          temp_max_x_sub <- 0.035
          hist(nucdiv_sub, col=ifelse(length(temp_seqs_bef) > length(temp_seqs_aft),"#8B735580", "#FF8C0090" ),border=ifelse(length(temp_seqs_bef) > length(temp_seqs_aft),"#8B735580", "#FF8C0090" ), main=paste(main_name, " (", Event_name, ")", sep=""), xlab="Nucleotide diversity")
          abline(v=c(mean(nucdiv_sub), mean(nucdiv_sub)+sd(nucdiv_sub), mean(nucdiv_sub)-sd(nucdiv_sub)), lwd=c(3,3,3), lty=c(2,3,3), col=ifelse(length(temp_seqs_bef) > length(temp_seqs_aft),"#8B735580", "#FF8C0090" ))
          abline(v=SS_bef[2], lwd=6, col="#8B7355")
          text(paste("Before", " (", length(temp_seqs_bef),")", sep="") ,x=temp_max_x_sub*0.96, y=temp_max_y_sub*0.8, col="#8B7355", adj=c(0,0.5))
          #segments(x0=temp_max_x_sub*0.95, x1=temp_max_x_sub*0.91, y0=temp_max_y_sub*0.8, y1=temp_max_y_sub*0.8, col="#8B7355", lwd=6)
          abline(v=SS_aft[2], lwd=6, col="#FF8C00")
          text(paste("After", " (", length(temp_seqs_aft),")", sep=""),x=temp_max_x_sub*0.96, y=temp_max_y_sub*0.75, col="#FF8C00", adj=c(0,0.5))
          #segments(x0=temp_max_x_sub*0.95, x1=temp_max_x_sub*0.91, y0=temp_max_y_sub*0.75, y1=temp_max_y_sub*0.75, col="#FF8C00", lwd=6)
        }
      #if (event == 2){
        #hist(nucdiv_sub, col=paste(col_event, 50, sep=""), border=NA, xlim=c(0,0.07), main=main_name)
        #abline(v=c(mean(nucdiv_sub), mean(nucdiv_sub)+sd(nucdiv_sub), mean(nucdiv_sub)-sd(nucdiv_sub)), lwd=c(3,2,2), lty=c(2,3,3), col=col_event)
        #abline(v=SS_bef[2], lwd=4, col=col_event)
        #abline(v=SS_aft[2], lwd=4, col="yellow") 
        
      }
    }

Sp_bef_aft[sp,] <- c(Sp_LP_H$Species[sp], round(SS_bef[2], digits=4), round(SS_aft[2], digits=4))  
}

}
#}












temp_plot_by <- as.data.frame(matrix(nrow=length(Sp_bef_aft[,1]),ncol=4))
temp_plot_by[,c(1,2,3)] <- Sp_bef_aft[,c(1,2,3)]
temp_plot_by[,4] <- as.numeric(Sp_bef_aft[,3]) - as.numeric(Sp_bef_aft[,2])
colnames(temp_plot_by) <- c(colnames(Sp_bef_aft), "Diff")
# sort by Nuc div before
Plot_bef_aft_by_bef <- temp_plot_by[order(as.numeric(temp_plot_by$Bef), decreasing=T),]
barplot(as.numeric(Plot_bef_aft_by_bef$Bef), col="#8B735590", ylim=c(0,0.08))
barplot(as.numeric(Plot_bef_aft_by_bef$Aft), col="#FF8C0090", add=T)
# sort by Nuc div after
Plot_bef_aft_by_bef <- temp_plot_by[order(as.numeric(temp_plot_by$Bef), decreasing=T),]
barplot(as.numeric(Plot_bef_aft_by_bef$Bef), col="#8B735590", ylim=c(0,0.08))
barplot(as.numeric(Plot_bef_aft_by_bef$Aft), col="#FF8C0090", add=T)
# sort by Nuc div difference
Plot_bef_aft_by_diff <- temp_plot_by[order(as.numeric(temp_plot_by$Diff), decreasing=F),]
barplot(as.numeric(Plot_bef_aft_by_diff$Bef), col="#8B735590", ylim=c(0,0.08), plot=F)
barplot(as.numeric(Plot_bef_aft_by_diff$Aft), col="#FF8C0090", add=T, names.arg=Plot_bef_aft_by_diff$Species, las=2)

plot(x=seq(1,dim(Plot_bef_aft_by_diff)[1], by=1), ylim=c(-0.01,0.082), xlim=c(0,18), frame = F, xaxt="n", xlab=NA, ylab=NA, yaxt="n")
width <- 0.9
for (nuc in seq_along(Plot_bef_aft_by_diff[,1])){
  if (as.numeric(Plot_bef_aft_by_diff[nuc,3]) > 0){
    if (as.numeric(Plot_bef_aft_by_diff[nuc,2]) > as.numeric(Plot_bef_aft_by_diff[nuc,3])){
      polygon(x=c(nuc,nuc,nuc+width,nuc+width), y=c(0,Plot_bef_aft_by_diff[nuc,3], Plot_bef_aft_by_diff[nuc,3],0), col="#009ACD99", border=NA)
      polygon(x=c(nuc,nuc,nuc+width,nuc+width), y=c(Plot_bef_aft_by_diff[nuc,3],Plot_bef_aft_by_diff[nuc,2],Plot_bef_aft_by_diff[nuc,2],Plot_bef_aft_by_diff[nuc,3]), col="#CDAA7D", border=NA)
    }
    if (as.numeric(Plot_bef_aft_by_diff[nuc,2]) < as.numeric(Plot_bef_aft_by_diff[nuc,3])){
      polygon(x=c(nuc,nuc,nuc+width,nuc+width), y=c(Plot_bef_aft_by_diff[nuc,3],Plot_bef_aft_by_diff[nuc,2],Plot_bef_aft_by_diff[nuc,2],Plot_bef_aft_by_diff[nuc,3]), col="#009ACD99", border=NA)
      polygon(x=c(nuc,nuc,nuc+width,nuc+width), y=c(0,Plot_bef_aft_by_diff[nuc,2], Plot_bef_aft_by_diff[nuc,2],0), col="#CDAA7D", border=NA)
    }
  }
  if (as.numeric(Plot_bef_aft_by_diff[nuc,3]) < 0){
    polygon(x=c(nuc,nuc,nuc+width,nuc+width), y=c(0,Plot_bef_aft_by_diff[nuc,2], Plot_bef_aft_by_diff[nuc,2],0), col="#CDAA7D", border=NA)
    text(x=nuc+width/2, y=0 ,"T", col="#009ACD", cex=2, adj=c(0.5,0))
    
  }
  
  sp_letters <- paste(strsplit(strsplit(Plot_bef_aft_by_diff[nuc,1], split=c("_"))[[1]][1], split="")[[1]][1], strsplit(strsplit(Plot_bef_aft_by_diff[nuc,1], split=c("_"))[[1]][2], split="")[[1]][1], sep="")
  text(y=-0.005, x=nuc+(width/3), sp_letters, adj=c(0,1), cex=0.8)
}
polygon(x=c(2,2,3,3), y=c(0.07,0.073,0.073,0.07), col="#009ACD99",border=NA)
text("Holocene", x=3.1, y=0.07, adj=c(0, -0.5), cex=1.2)
polygon(x=c(2,2,3,3), y=c(0.066,0.069,0.069,0.066), col="#CDAA7D",border=NA)
text("Late Pleistocene", x=3.1, y=0.066, adj=c(0, -0.5), cex=1.2)
abline(h=c(0,0.01, 0.02, 0.04, 0.08))
text(x=0.7, y=c(0,0.01, 0.02, 0.04, 0.08),  c(0,0.01, 0.02, 0.04, 0.08), cex=0.7, adj=c(0.5,-0.5))



#################### delete bleow this line ###############

