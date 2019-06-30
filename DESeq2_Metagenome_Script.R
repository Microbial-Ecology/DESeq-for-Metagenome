#Metagenome DESeq2 analysis

install.packages("BiocManager")
install.packages("ggplot2")
install.packages("plyr")
install.packages("ape")
install.packages("picante")
install.packages("tidyverse")

library(BiocManager)

BiocManager::install("phyloseq")
BiocManager::install("DESeq2")


library(ggplot2)
library(phyloseq)
library(DESeq2)
library(plyr)
library(ape)
library(picante)
library(tidyverse)
library(rjson)

otu<-import_biom("Taxamerge.biom")
map<-import_qiime_sample_data("RottnestFebMap.txt")
sample_names(otu)
sample_names(map)
J <- merge_phyloseq(otu, map)

save(J,file="Taxa.phyloseq")

ntaxa(J)
rank_names(J)
sample_variables(J)

#The biom files name the ranks as Rank1, Rank 2, Rank 3, Rank 4, Rank 5, Rank 6 and Rank 7. To change it to Kingdom, Phylum, ..., Species - run the command bellow
colnames(tax_table(J)) = c("Kingdom","Phylum", "Class", "Order", "Family", "Genus", "Species")

#DESeq2 transformation

packageVersion("DESeq2")

#cut alpha level at 0.05. We might want to change this later to higher level of significance if we get heaps of otu responding

alpha = 0.05

#alpha can be changed to 0.005  or 0.01 if necessary

alpha = 0.005
alpha = 0.01


#set up model Lithifying vs Nonlithifying
#sets it up so that the first factor will be on the bottom, so > 0 means it is higher in Nonlithifying < 0 means higher in Lithifying
sample_data(J)$Type <- factor(sample_data(J)$Type, levels = c("Lithifying", "Nonlithifying"))

#look at lith Vs nonlith
dds<-phyloseq_to_deseq2(J,~Type)
dds<-DESeq(dds,test="Wald",fitType = "parametric")
com<-results(dds)

#pull out model parameters for plotting
sigtab = com[which(com$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(J)[rownames(sigtab), ], "matrix"))

test=as.data.frame(com)
#we can write out the above as a csv if we want to join and save it with other data later on.

# order the graph (Genus). You can change this to other taxonomic levels too
x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))

#graph
DA<-ggplot(sigtab, aes(y=Genus, x=log2FoldChange)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 1)+geom_point(aes(fill=Phylum), size=3.5,color="black",shape=21) + 
  theme(legend.key = element_rect(fill="white"), panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"), panel.background = element_rect(fill="white",colour="grey50"), axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

DA<-DA+theme(panel.grid.minor.x=element_blank(),
             panel.grid.major.x=element_blank())
DA

#set up model pair-wise mat morphologies

##B1

B1<-subset_samples(J,Mat== "Blister"| Mat== "Flocculent")
sample_data(B1)

#set up model Blister vs Flocculent
sample_data(B1)$Mat <- factor(sample_data(B1)$Mat, levels = c("Blister", "Flocculent"))
levels(sample_data(B1)$Mat)

#look at Blister vs Flocculent
dds1<-phyloseq_to_deseq2(B1,~Mat)
dds1<-DESeq(dds1,test="Wald",fitType = "parametric")
com1<-results(dds1)

#pull out model parameters for plotting
sigtab = com1[which(com1$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(B1)[rownames(sigtab), ], "matrix"))

test=as.data.frame(com1)
#we can write out the above as a csv if we want to join and save it with other data later on.

# order the graph (Family). You can change this to other taxonomic levels too
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

#graph
DA1<-ggplot(sigtab, aes(y=Family, x=log2FoldChange)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 1)+geom_point(aes(fill=Phylum), size=3.5,color="black",shape=21) + 
  theme(legend.key = element_rect(fill="white"), panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"), panel.background = element_rect(fill="white",colour="grey50"), axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

DA1<-DA1+theme(panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank())
DA1


##B2

B2<-subset_samples(J,Mat== "Blister"| Mat== "Pustular")
sample_data(B2)

#set up model Blister vs Pustular
sample_data(B2)$Mat <- factor(sample_data(B2)$Mat, levels = c("Blister", "Pustular"))
levels(sample_data(B2)$Mat)

#look at Blister vs Pustular
dds2<-phyloseq_to_deseq2(B2,~Mat)
dds2<-DESeq(dds2,test="Wald",fitType = "parametric")
com2<-results(dds2)

#pull out model parameters for plotting
sigtab = com2[which(com2$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(B2)[rownames(sigtab), ], "matrix"))

test=as.data.frame(com2)
#we can write out the above as a csv if we want to join and save it with other data later on.

# order the graph (Family). You can change this to other taxonomic levels too
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

#graph
DA2<-ggplot(sigtab, aes(y=Family, x=log2FoldChange)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 1)+geom_point(aes(fill=Phylum), size=3.5,color="black",shape=21) + 
  theme(legend.key = element_rect(fill="white"), panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"), panel.background = element_rect(fill="white",colour="grey50"), axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

DA2<-DA2+theme(panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank())
DA2


##B3

B3<-subset_samples(J,Mat== "Blister"| Mat== "Nonlithifying cohesive")
sample_data(B3)

#set up model Blister vs NC
sample_data(B3)$Mat <- factor(sample_data(B3)$Mat, levels = c("Blister", "Nonlithifying cohesive"))
levels(sample_data(B3)$Mat)

#look at Blister vs NC
dds3<-phyloseq_to_deseq2(B3,~Mat)
dds3<-DESeq(dds3,test="Wald",fitType = "parametric")
com3<-results(dds3)

#pull out model parameters for plotting
sigtab = com3[which(com3$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(B3)[rownames(sigtab), ], "matrix"))

test=as.data.frame(com3)
#we can write out the above as a csv if we want to join and save it with other data later on.

# order the graph (Family). You can change this to other taxonomic levels too
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

#graph
DA3<-ggplot(sigtab, aes(y=Family, x=log2FoldChange)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 1)+geom_point(aes(fill=Phylum), size=3.5,color="black",shape=21) + 
  theme(legend.key = element_rect(fill="white"), panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"), panel.background = element_rect(fill="white",colour="grey50"), axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

DA3<-DA3+theme(panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank())
DA3


##B4

B4<-subset_samples(J,Mat== "Blister"| Mat== "Nonlithifying loosely cohesive")
sample_data(B4)

#set up model Blister vs NLC
sample_data(B4)$Mat <- factor(sample_data(B4)$Mat, levels = c("Blister", "Nonlithifying loosely cohesive"))
levels(sample_data(B4)$Mat)

#look at Blister vs NLC
dds4<-phyloseq_to_deseq2(B4,~Mat)
dds4<-DESeq(dds4,test="Wald",fitType = "parametric")
com4<-results(dds4)

#pull out model parameters for plotting
sigtab = com4[which(com4$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(B4)[rownames(sigtab), ], "matrix"))

test=as.data.frame(com4)
#we can write out the above as a csv if we want to join and save it with other data later on.

# order the graph (Family). You can change this to other taxonomic levels too
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

#graph
DA4<-ggplot(sigtab, aes(y=Family, x=log2FoldChange)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 1)+geom_point(aes(fill=Phylum), size=3.5,color="black",shape=21) + 
  theme(legend.key = element_rect(fill="white"), panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"), panel.background = element_rect(fill="white",colour="grey50"), axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

DA4<-DA4+theme(panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank())
DA4


##B5

B5<-subset_samples(J,Mat== "Flocculent"| Mat== "Pustular")
sample_data(B5)

#set up model Flocculent vs Pustular
sample_data(B5)$Mat <- factor(sample_data(B5)$Mat, levels = c("Flocculent", "Pustular"))
levels(sample_data(B5)$Mat)

#look at Flocculent vs Pustular
dds5<-phyloseq_to_deseq2(B5,~Mat)
dds5<-DESeq(dds5,test="Wald",fitType = "parametric")
com5<-results(dds5)

#pull out model parameters for plotting
sigtab = com5[which(com5$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(B5)[rownames(sigtab), ], "matrix"))

test=as.data.frame(com5)
#we can write out the above as a csv if we want to join and save it with other data later on.

# order the graph (Family). You can change this to other taxonomic levels too
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

#graph
DA5<-ggplot(sigtab, aes(y=Family, x=log2FoldChange)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 1)+geom_point(aes(fill=Phylum), size=3.5,color="black",shape=21) + 
  theme(legend.key = element_rect(fill="white"), panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"), panel.background = element_rect(fill="white",colour="grey50"), axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

DA5<-DA5+theme(panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank())
DA5

##B6

B6<-subset_samples(J,Mat== "Flocculent"| Mat== "Nonlithifying cohesive")
sample_data(B6)

#set up model Flocculent vs NC
sample_data(B6)$Mat <- factor(sample_data(B6)$Mat, levels = c("Flocculent", "Nonlithifying cohesive"))
levels(sample_data(B6)$Mat)

#look at Flocculent vs NC
dds6<-phyloseq_to_deseq2(B6,~Mat)
dds6<-DESeq(dds6,test="Wald",fitType = "parametric")
com6<-results(dds6)

#pull out model parameters for plotting
sigtab = com6[which(com6$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(B6)[rownames(sigtab), ], "matrix"))

test=as.data.frame(com6)
#we can write out the above as a csv if we want to join and save it with other data later on.

# order the graph (Family). You can change this to other taxonomic levels too
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

#graph
DA6<-ggplot(sigtab, aes(y=Family, x=log2FoldChange)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 1)+geom_point(aes(fill=Phylum), size=3.5,color="black",shape=21) + 
  theme(legend.key = element_rect(fill="white"), panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"), panel.background = element_rect(fill="white",colour="grey50"), axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

DA6<-DA6+theme(panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank())
DA6


##B7

B7<-subset_samples(J,Mat== "Flocculent"| Mat== "Nonlithifying loosely cohesive")
sample_data(B7)

#set up model Flocculent vs NLC
sample_data(B7)$Mat <- factor(sample_data(B7)$Mat, levels = c("Flocculent", "Nonlithifying loosely cohesive"))
levels(sample_data(B7)$Mat)

#look at Flocculent vs NLC
dds7<-phyloseq_to_deseq2(B7,~Mat)
dds7<-DESeq(dds7,test="Wald",fitType = "parametric")
com7<-results(dds7)

#pull out model parameters for plotting
sigtab = com7[which(com7$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(B7)[rownames(sigtab), ], "matrix"))

test=as.data.frame(com7)
#we can write out the above as a csv if we want to join and save it with other data later on.

# order the graph (Family). You can change this to other taxonomic levels too
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

#graph
DA7<-ggplot(sigtab, aes(y=Family, x=log2FoldChange)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 1)+geom_point(aes(fill=Phylum), size=3.5,color="black",shape=21) + 
  theme(legend.key = element_rect(fill="white"), panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"), panel.background = element_rect(fill="white",colour="grey50"), axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

DA7<-DA7+theme(panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank())
DA7


##B8

B8<-subset_samples(J,Mat== "Pustular"| Mat== "Nonlithifying cohesive")
sample_data(B8)

#set up model Pustular vs NC
sample_data(B8)$Mat <- factor(sample_data(B8)$Mat, levels = c("Pustular", "Nonlithifying cohesive"))
levels(sample_data(B8)$Mat)

#look at Pustular vs NLC
dds8<-phyloseq_to_deseq2(B8,~Mat)
dds8<-DESeq(dds8,test="Wald",fitType = "parametric")
com8<-results(dds8)

#pull out model parameters for plotting
sigtab = com8[which(com8$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(B8)[rownames(sigtab), ], "matrix"))

test=as.data.frame(com8)
#we can write out the above as a csv if we want to join and save it with other data later on.

# order the graph (Family). You can change this to other taxonomic levels too
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

#graph
DA8<-ggplot(sigtab, aes(y=Family, x=log2FoldChange)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 1)+geom_point(aes(fill=Phylum), size=3.5,color="black",shape=21) + 
  theme(legend.key = element_rect(fill="white"), panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"), panel.background = element_rect(fill="white",colour="grey50"), axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

DA8<-DA8+theme(panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank())
DA8


##B9

B9<-subset_samples(J,Mat== "Pustular"| Mat== "Nonlithifying loosely cohesive")
sample_data(B9)

#set up model Pustular vs NLC
sample_data(B9)$Mat <- factor(sample_data(B9)$Mat, levels = c("Pustular", "Nonlithifying loosely cohesive"))
levels(sample_data(B9)$Mat)

#look at Pustular vs NLC
dds9<-phyloseq_to_deseq2(B9,~Mat)
dds9<-DESeq(dds9,test="Wald",fitType = "parametric")
com9<-results(dds9)

#pull out model parameters for plotting
sigtab = com9[which(com9$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(B9)[rownames(sigtab), ], "matrix"))

test=as.data.frame(com9)
#we can write out the above as a csv if we want to join and save it with other data later on.

# order the graph (Family). You can change this to other taxonomic levels too
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

#graph
DA9<-ggplot(sigtab, aes(y=Family, x=log2FoldChange)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 1)+geom_point(aes(fill=Phylum), size=3.5,color="black",shape=21) + 
  theme(legend.key = element_rect(fill="white"), panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"), panel.background = element_rect(fill="white",colour="grey50"), axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

DA9<-DA9+theme(panel.grid.minor.x=element_blank(),
               panel.grid.major.x=element_blank())
DA9


##B10

B10<-subset_samples(J,Mat== "Nonlithifying cohesive"| Mat== "Nonlithifying loosely cohesive")
sample_data(B10)

#set up model NC vs NLC
sample_data(B10)$Mat <- factor(sample_data(B10)$Mat, levels = c("Nonlithifying cohesive", "Nonlithifying loosely cohesive"))
levels(sample_data(B10)$Mat)

#look at Pustular vs NLC
dds10<-phyloseq_to_deseq2(B10,~Mat)
dds10<-DESeq(dds10,test="Wald",fitType = "parametric")
com10<-results(dds10)

#pull out model parameters for plotting
sigtab = com10[which(com10$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(B10)[rownames(sigtab), ], "matrix"))

test=as.data.frame(com10)
#we can write out the above as a csv if we want to join and save it with other data later on.

# order the graph (Family). You can change this to other taxonomic levels too
x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))

#graph
DA10<-ggplot(sigtab, aes(y=Family, x=log2FoldChange)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 1)+geom_point(aes(fill=Phylum), size=3.5,color="black",shape=21) + 
  theme(legend.key = element_rect(fill="white"), panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey"), panel.background = element_rect(fill="white",colour="grey50"), axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

DA10<-DA10+theme(panel.grid.minor.x=element_blank(),
                 panel.grid.major.x=element_blank())
DA10