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
library(dplyr)
library(stringr)


setwd("~/lfainoData/XylellaOct_2019/")
data2<-read.csv("MLSTVale.fastq.txt", header = T)
data3 <- select(data2, tribu, A)
data3 = data3[-1,]
colnames(data3) <- c("tribu", "Reads")
data3$strain <- str_replace(data3$tribu, "t__", "")
match<-read.delim("../Xylella/repository/xylella/corrispondance.txt")
df4<-merge(data3, match, by.x = "strain", by.y = "strain", all.x=TRUE)
df4$Subspecies <- as.character(df4$Subspecies)
df4$Subspecies[is.na(df4$Subspecies)] <- "unknown"
df4$Samples <- as.character("1")
col <- c("pauca"="#FF0000", "sandyi"="#DFAA32", "multiplex"="#808080", "fastidiosa"="#0000FF", "morus"="#808000", "unknown" = "black") 
theme_update(plot.title = element_text(hjust = 0.5, size = 22), 
             axis.text.x = element_text(angle = 90, hjust = 1),
             legend.text =element_text(size = 15),legend.title=element_text(size=18),
             axis.text=element_text(size=8), axis.title=element_text(size=15))
ggplot(df4, aes(fill=Subspecies, y=Reads, x=Samples)) +  geom_bar(stat="identity", position = "fill") + scale_fill_manual(values = col) +
  ggtitle("Olives samples gDNA Xf") +scale_x_discrete(limits = levels(as.factor(df4$Samples)))

