################################################################################

setwd("~/Documents/MAC_DUI/SBD_MS")
library("radiator")

fm<-ape::read.dna(file="MACEU_cox1f_males_females.fasta",format="fasta")
pop<-c(rep("M",180),rep("F",166))
lab<-labels(fm)


gen<-mmod::as.genind.DNAbin(fm,pop)
adegenet::indNames(gen)<-lab
adegenet::locNames(gen)<-"cox1f"

radiator::detect_genomic_format(data = gen)

require(strataG) # for the gtypes format...
genepop_data <- genomic_converter(
  data = gen , 
  output =  "genepop")

#write.table(genepop_data, file="MACEU_cox1f_males_females.genepop")

str(genepop_data)

library(prabclus)

####### 2.1 - Grp1 = N-Atl Grp2 = Rif

geographic_coordinates <- read.csv("MACEU_fm_popcoord.csv", sep=",", header=T)
samples <- geographic_coordinates[,1]
popmap <- geographic_coordinates[,4]
geographic_coordinates <- geographic_coordinates[, 2:3]

geographic_distances <- coord2dist(coordmatrix=geographic_coordinates, file.format = "decimal2")
geographic_distances <- geographic_distances + quantile(geographic_distances)[2]
geographic_distances <- log(geographic_distances)

# Import genetic distances 1: a distances calculated in genepop

library(genepop)
library(tidyverse)

gene_ibd<-read.table("MACEU_fm_genepop_ibd.txt",header = F)
pop<-as.data.frame(rep('pop',346))
 rows_insert(gene_ibd,pop,by=2)

gene_idbformat<-gene_ibd[rep(1:nrow(gene_ibd),"pop",each=2),]

d<-ibd(inputFile = "MACEU_fm_genepop_ibd_format.txt",outputFile = "ibd_a_like", dataType = 'haploid',statistic = 'a-like',bootstrapMethod = 'BCa')

ncol <- max(count.fields("./1-Genet_dist_AtlN_Rif.txt", sep = " "))
genetic_distances_a <- as.matrix(read.table("./1-Genet_dist_a_AtlN_Rif.txt", sep=" ", fill=T, col.names=paste0('V', seq_len(ncol))))
for(i in 2:ncol(genetic_distances_a)){
  COL <- max(which(is.na(genetic_distances_a[,i])))
  genetic_distances_a[COL,i] <- 0
}
genetic_distances_a[upper.tri(genetic_distances_a)] <- t(genetic_distances_a)[upper.tri(genetic_distances_a)]
colnames(genetic_distances_a) <- samples
rownames(genetic_distances_a) <- samples

####analyses de rarefaction #####

library(spider)

fm<-ape::read.dna(file="MACEU_cox1f_males_females.fasta",format="fasta")
pop_sex<-c(rep("M",180),rep("F",166))
lab<-labels(fm)
pop_pop<-c(rep("AYT_M",33),rep("BRE_M",29), rep("BRI_M",30),rep("CRO_M",28), rep("FOU_M",30),rep("MAH_M",30),
           rep("AYT_F",31), rep("CRO_F",20), rep("FOU_F",30), rep("MAH_F",29), rep("BRI_F",30), rep("BRE_F",26))


gen_pop<-mmod::as.genind.DNAbin(fm,pop_pop)
adegenet::indNames(gen)<-lab
adegenet::locNames(gen)<-"cox1f"

gen_sex<-mmod::as.genind.DNAbin(fm,pop_sex)
adegenet::indNames(gen)<-lab
adegenet::locNames(gen)<-"cox1f"



