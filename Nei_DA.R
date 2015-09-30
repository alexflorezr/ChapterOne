rm(list=ls())
setwd("C:/Katies Data/Amphibians/Popgentest/EastvWest")
#library(seqinr)  #might have to download .tgz file directly from CRAN site and install locally, not directly from CRAN repository
install.packages("ape")
install.packages("pegas")
install.packages("strataG")
install.packages("sidier")

library(ape)  
library(pegas)
library(strataG)
library(sidier)


#library(mmod)
library(adegenet)
#library(plyr)
library(strataG)
#library(iNEXT)
library(sidier)
library(reshape2)


######### What do I need?



##Alignment level:

#nucleotide diversity  ape, DNAbin(?), or strataG gtypes--define bases so only using letters
#haplotype diversity #may be able to do with pop genome, or modify dipnet? probably use code from below
#sidier to get unique haps, DNAbin in and out, writes fasta
#haplotypic.diversity in strataG? gtypes
#Tajima's D   pegas, DNAbin
#Fu's Fs  ##may be possible in PopGenome using paranthetical lists of populations? see instr. manual
#Global Fst
#Global Phist   hopefully in strataG?


##Sequence level:

#pairwise genetic distance matrix   ape
#geographic distance matrix   ape
#cost distance matrix

##Population level:

#Fst
#Phist
#Jost D for pop divergence?
#geographic distance matrix   ape
#cost distance matrix


####Next steps after that
#cluster/procrustes analysis of alignment level stats
#isolation by distance or replacement
#within versus between lineage cost distance, relative to geographic distance


###############################################


#Get data in
setwd("~/Desktop/Ch1_Figures/Figure_one/Same_length_seqs/")
temp_fasta <- read.dna(paste(tolower(temp_sp), "_fasta_new.fasta", sep=""), format = "fasta", as.character = FALSE, as.matrix=NULL)
# KATIES sp_align<-read.dna("Ensatina_eschscholtzii_align.fasta", format = "fasta", as.character = FALSE, as.matrix=NULL) #gets seq data in as DNAbin
sp_data<-read.csv("Ensatina_eschscholtzii.csv", stringsAsFactors = F)
na.omit(sp_data)
sp_data <- sp_data[match(rownames(sp_align), sp_data$UniqueInd),]#reorder to match alignment order
N_seqs<-length(sp_data$UniqueInd)
temp_data <- as.data.frame(matrix(nrow=191, ncol=3))
for (name in 1:dim(temp_data)[1]){
  uno <- strsplit(labels(temp_fasta), spli="_")[[name]][c(1,2,3)]
  temp_data[name,1] <- paste(uno[1], uno[2], uno[3], sep="_")
  temp_data[name,2] <- as.numeric(strsplit(labels(temp_fasta), spli="_")[[name]][4])
  if( as.numeric(strsplit(labels(temp_fasta), spli="_")[[name]][4]) > 11700){
    temp_data[name,3] <- "Pleistocene"
  }else{ temp_data[name,3] <- "Holocene"
  }
}
colnames(temp_data) <- c("name", "age", "era")
###get unique haps using sidier
##combine haplotype list with species csv info

haps<-FindHaplo(align=temp_fasta, saveFile=F)
as.data.frame(haps)
colnames(haps)<-c("UniqueInd", "haplotype")
full_data<-merge(temp_data, haps)

N_haps<-length(unique(full_data$haplotype))###add number of haplotypes as one of statistics to be output

sp_haps<-GetHaplo(align=sp_align, saveFile=T, outname="Ensatina_eschscholtzii_haps.fasta", format="fasta", seqsNames="Inf.Hap") #haps are now a DNAbin
sp_ghaps<-read.fasta("Ensatina_eschscholtzii_haps.fasta") #imports haps as gtypes
#str(sp_ghaps)
#View(full_data)

#create gtypes object for strataG
#this line borrowed from DIPnet
sp_gtype<-gtypes(gen.data=data.frame(full_data$UniqueInd,full_data$population,full_data$haplotype),id.col=1,strata.col=2,locus.col=3,dna.seq=sp_ghaps)




##Alignment level statistics


sp_nuc.div<-nuc.div(sp_align, variance = F, pairwise.deletion = FALSE) #nucleotide diversity, pegas

sp_hap.div<-haplotypic.diversity(sp_gtype) #haplotype diversity, strataG


###############getting very different estimates of TajD and signinficance!!!!
#DIPnet uses pegas version
sp_tajD<-tajima.test(sp_align) #Tajima's D, pegas, #Pval.beta - p-value assuming that D follows a beta distribution (Tajima, 1989)


TD<-tajimas.d(sp_gtype) #strataG

tst<-read.fasta("Ensatina_eschscholtzii_align.fasta")
TD2<-tajimas.d(tst)strataG
#####################

##Ramos-Onsins-Rozas Test of Neutrality
#png(paste(name,"_R2.png"))
#sp_R2<-R2.test(sp_align, B=10000, plot=TRUE, quiet=TRUE)
#dev.off()
sp_R2<-R2.test(sp_align, B=10000, plot=FALSE, quiet=TRUE) #Ramos-Onsins-Rozas Test of Neutrality, pegas


glbl_phist<-stat.phist(sp_gtype)   #Overall Phi_st, strataG

glbl_fst<-stat.fst(sp_gtype)

sp_result<-c(N_seqs, N_haps, sp_nuc.div, sp_hap.div, sp_tajD$D, sp_tajD$Pval.beta, sp_R2$R2, sp_R2$P.val, glbl_phist, glbl_fst)
sp_result


###########Sequence level statistics

sp_gendist<-dist.dna(sp_align, model = "K80", variance = FALSE, gamma = FALSE, pairwise.deletion = FALSE, base.freq = NULL, as.matrix = FALSE)
pairwise.gendist<-melt(as.matrix(sp_gendist), varnames = c("row", "col"), value.name="gendist") #using reshape2
sp_eucldist<-geod(sp_data$long, sp_data$lat, R=6371) #gives distance in kilometers
diffs<-as.dist(t(sp_eucldist), diag=F, upper=F)
pairwise.eucldist<-melt(as.matrix(diffs), varnames = c("row", "col"), value.name="eucldist") 

sp.distances<-c(pairwise.gendist, pairwise.eucldist)
View(sp.distances)




#need to add something to all eucldists to make all non-zero?
#need to add cost distance matrix



##this will probably have to change so that eucldist is on same dist matrix as cost 
##need to add something to eucldist and cost so that all values are non-zero? 
#or does that only matter when applying a mantel test or doing the IBD regression?

