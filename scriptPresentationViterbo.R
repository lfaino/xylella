#install.packages("gtools")
library("xlsx")
library("ggplot2")
library(plyr)
require(graphics)
library(randomcoloR)
library(stringi)
library(gtools)
library("randomcoloR")
#install.packages("caTools")
library(RColorBrewer)
library(seqinr)
library(gplots)
setwd("~/lfainoData/Xylella")

##subspecies gDNA
data1<-read.delim("Xf_first_bast/Xf_first_bast.txt")
data2<-read.delim("regression/regression.txt")
data<-rbind(data1,data2)
data$Samples<-paste(data$Subset,data$Barcode,sep="-")
match<-read.delim("corrispondance.txt")
df4<-merge(data, match, by.x = "strain", by.y = "strain", all.x=TRUE)
col <- c("pauca"="#FF0000", "sandy"="#32CD32", "multiplex"="#808080", "fastidiosa"="#0000FF", "morus"="#808000", "NA" = "black") 
theme_update(plot.title = element_text(hjust = 0.5, size = 22), 
             axis.text.x = element_text(angle = 90, hjust = 1),
             legend.text =element_text(size = 15),legend.title=element_text(size=18),
             axis.text=element_text(size=8), axis.title=element_text(size=15))
pdf(file = "Xylella_gDNA_Xf_1.pdf", paper = "a4r")
ggplot(df4, aes(fill=Subspecies, y=Reads, x=Samples)) +  geom_bar(stat="identity", position = "fill") + scale_fill_manual(values = col) +
  ggtitle("Olives samples gDNA Xf") +scale_x_discrete(limits = levels(as.factor(df4$Samples)))
dev.off()
aggregate(Reads ~ Samples + Subspecies, data = df4, FUN = 'sum')

##Species LOD PCR
data3<-read.delim("LODpcr/LODpcr.txt")
data3 <- transform(data3, Samples = factor(Samples, levels = mixedsort(levels(data3$Samples))))
palette <- distinctColorPalette(100)
theme_update(plot.title = element_text(hjust = 0.5, size = 22), 
             axis.text.x = element_text(angle = 90, hjust = 1),
             legend.text =element_text(size = 6),legend.title=element_text(size=8),
             axis.text=element_text(size=8), axis.title=element_text(size=15))
pdf(file = "LOD_PCR_Xf_1.pdf", height = 8, width = 16)
ggplot(data3, aes(fill=Species, y=Reads, x=Samples)) +  geom_bar(stat="identity", position = "fill") + #scale_fill_manual(values = col) +
  ggtitle("LOD PCR Xf") +scale_x_discrete(limits = rev(levels(as.factor(data3$Samples)))) + #scale_fill_brewer(palette = "Dark2") # +scale_color_brewer(palette="Dark2")
  scale_fill_manual(values=palette) #wes.palette(n=3, name="GrandBudapest"))
dev.off()

##subspecies LOD consensus
data60<-read.delim("LODpcr/summary.MLST.txt")
match<-read.delim("corrispondance.txt")
df50<-merge(data60, match, by.x = "strain", by.y = "strain", all.x=TRUE)
df50 <- transform(df50, Samples = factor(Samples, levels = mixedsort(levels(df50$Samples))))
col <- c("pauca"="#FF0000", "sandy"="#32CD32", "multiplex"="#808080", "fastidiosa"="#0000FF", "morus"="#808000", "NA" = "black") 
theme_update(plot.title = element_text(hjust = 0.5, size = 22), 
             axis.text.x = element_text(angle = 90, hjust = 1),
             legend.text =element_text(size = 15),legend.title=element_text(size=18),
             axis.text=element_text(size=8), axis.title=element_text(size=15))
pdf(file = "LOD_PCR_Xfsubsp_1.pdf", height = 8, width = 16)
ggplot(df50, aes(fill=Subspecies, y=Reads, x=Samples)) +  geom_bar(stat="identity", position = "fill") + scale_fill_manual(values = col) +
  ggtitle("Olives samples gDNA Xf") +scale_x_discrete(limits = rev(levels(as.factor(df50$Samples))))
dev.off()
#aggregate(Reads ~ Samples + Subspecies, data = df4, FUN = 'sum')



### syntetic

data5<-read.delim(file = "subset.txt", header = FALSE)
ggplot(data5, aes(x=V1, y=V2)) + geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) + scale_y_log10()


##subspecies gDNA
data6<-read.delim("LODpcr/LODpcrXf.txt")
match<-read.delim("corrispondance.txt")
df5<-merge(data6, match, by.x = "strain", by.y = "strain", all.x=TRUE)
df5 <- transform(df5, Samples = factor(Samples, levels = mixedsort(levels(df5$Samples))))
col <- c("pauca"="#FF0000", "sandy"="#32CD32", "multiplex"="#808080", "fastidiosa"="#0000FF", "morus"="#808000", "NA" = "black") 
theme_update(plot.title = element_text(hjust = 0.5, size = 22), 
             axis.text.x = element_text(angle = 90, hjust = 1),
             legend.text =element_text(size = 15),legend.title=element_text(size=18),
             axis.text=element_text(size=8), axis.title=element_text(size=15))
pdf(file = "LOD_PCR_Xfsubsp_1.pdf", height = 8, width = 16)
ggplot(df5, aes(fill=Subspecies, y=Reads, x=Samples)) +  geom_bar(stat="identity", position = "fill") + scale_fill_manual(values = col) +
  ggtitle("Olives samples gDNA Xf") +scale_x_discrete(limits = rev(levels(as.factor(df5$Samples))))
dev.off()
#aggregate(Reads ~ Samples + Subspecies, data = df4, FUN = 'sum')

##subspecies unknown gDNA
data7<-read.delim("unknownS/gDNaUnk.txt")
match<-read.delim("corrispondance.txt")
df6<-merge(data7, match, by.x = "strain", by.y = "strain", all.x=TRUE)
df6 <- transform(df6, Samples = factor(Samples, levels = mixedsort(levels(df6$Samples))))
col <- c("pauca"="#FF0000", "sandy"="#32CD32", "multiplex"="#808080", "fastidiosa"="#0000FF", "morus"="#808000", "NA" = "black") 
theme_update(plot.title = element_text(hjust = 0.5, size = 22), 
             axis.text.x = element_text(angle = 90, hjust = 1),
             legend.text =element_text(size = 15),legend.title=element_text(size=18),
             axis.text=element_text(size=8), axis.title=element_text(size=15))
pdf(file = "UnknowndirectgDNA_1.pdf", height = 8, width = 16)
ggplot(df6, aes(fill=Subspecies, y=Reads, x=Samples)) +  geom_bar(stat="identity", position = "fill") + scale_fill_manual(values = col) +
  ggtitle("Olives samples gDNA Xf") +scale_x_discrete(limits = rev(levels(as.factor(df6$Samples))))
dev.off()
#aggregate(Reads ~ Samples + Subspecies, data = df4, FUN = 'sum')

##species unknown gDNA
data30<-read.delim("unknownS/summary.Bacteria.gDNA.txt")
data30 <- transform(data30, Samples = factor(Samples, levels = mixedsort(levels(data30$Samples))))
palette <- distinctColorPalette(100)
theme_update(plot.title = element_text(hjust = 0.5, size = 22), 
             axis.text.x = element_text(angle = 90, hjust = 1),
             legend.text =element_text(size = 6),legend.title=element_text(size=8),
             axis.text=element_text(size=8), axis.title=element_text(size=15))
pdf(file = "LOD_PCR_Xf_1.pdf", height = 8, width = 16)
ggplot(data30, aes(fill=Species, y=Reads, x=Samples)) +  geom_bar(stat="identity", position = "fill") + #scale_fill_manual(values = col) +
  ggtitle("LOD PCR Xf") + scale_x_discrete(limits = rev(levels(as.factor(data30$Samples)))) + #scale_fill_brewer(palette = "Dark2") # +scale_color_brewer(palette="Dark2")
  scale_fill_manual(values=palette) #wes.palette(n=3, name="GrandBudapest"))
dev.off()


##subspecies unknown PCR
data8<-read.delim("unknownPCRs/output.txt")
match<-read.delim("corrispondance.txt")
df7<-merge(data8, match, by.x = "strain", by.y = "strain", all.x=TRUE)
df7 <- transform(df7, Samples = factor(Samples, levels = mixedsort(levels(df7$Samples))))
col <- c("pauca"="#FF0000", "sandyi"="#32CD32", "multiplex"="#808080", "fastidiosa"="#0000FF", "morus"="#808000", "NA" = "black") 
theme_update(plot.title = element_text(hjust = 0.5, size = 22), 
             axis.text.x = element_text(angle = 90, hjust = 1),
             legend.text =element_text(size = 15),legend.title=element_text(size=18),
             axis.text=element_text(size=8), axis.title=element_text(size=15))
pdf(file = "UnknowndirectPCR_1.pdf", height = 8, width = 16)
ggplot(df7, aes(fill=Subspecies, y=Reads, x=Samples)) +  geom_bar(stat="identity", position = "fill") + scale_fill_manual(values = col) +
  ggtitle("Olives samples gDNA Xf") +scale_x_discrete(limits = rev(levels(as.factor(df7$Samples)))) #+ scale_y_log10()
dev.off()
#aggregate(Reads ~ Samples + Subspecies, data = df4, FUN = 'sum')


##subspecies unknown PCR polish
data10<-read.delim("unknownPCRs/MLST_2/summary.1.txt")
match<-read.delim("corrispondance.txt")
df10<-merge(data10, match, by.x = "strain", by.y = "strain", all.x=TRUE)
df10 <- transform(df10, Samples = factor(Samples, levels = mixedsort(levels(df10$Samples))))
col <- c("pauca"="#FF0000", "sandyi"="#32CD32", "multiplex"="#808080", "fastidiosa"="#0000FF", "morus"="#808000", "NA" = "black") 
theme_update(plot.title = element_text(hjust = 0.5, size = 22), 
             axis.text.x = element_text(angle = 90, hjust = 1),
             legend.text =element_text(size = 15),legend.title=element_text(size=18),
             axis.text=element_text(size=8), axis.title=element_text(size=15))
pdf(file = "UnknowndirectPCR_2.pdf", height = 8, width = 16)
ggplot(df10, aes(fill=Subspecies, y=Reads, x=Samples)) +  geom_bar(stat="identity", position = "fill") + scale_fill_manual(values = col) +
  ggtitle("Olives samples gDNA Xf") +scale_x_discrete(limits = rev(levels(as.factor(df10$Samples)))) #+ scale_y_log10()
dev.off()

##matrix align

myseqs <- read.alignment("consensus/fromDatabase_sanger_nanopore.aln.fas", format = "fasta")
m<-as.matrix(dist.alignment(myseqs, matrix = "identity" ))

pdf(file = "MatchConsensus.pdf", height = 16, width = 16)

heatmap.2(m)

dev.off()

write.csv(m, file = "test.txt")

##matrix align zoom


myseqs <- read.alignment("consensus/87_sanger_nanopore.aln", format = "fasta")
m<-as.matrix(dist.alignment(myseqs, matrix = "identity" ))

pdf(file = "MatchConsensus.zoom.pdf", height = 16, width = 16)

heatmap.2(m)

dev.off()
