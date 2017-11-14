#############
## Figure1 ##
#############
load("Ptep_top_RNAseq_obj.Rdata")
# make new column for total siRNA abundance from each TE
Ptep_top_RNAseq[,8]<-Ptep_top_RNAseq[,3]+Ptep_top_RNAseq[,4]
# select 15th-highest siRNA abundance value and set as threshold
siRNA_threshold<-sort(Ptep_top_RNAseq[,8],decreasing=T)[15]
# mark whether each TEs is above (TRUE) or below (FALSE) the threshold
siRNAtop<-c()
for (TE in Ptep_top_RNAseq[,8]) {
  if (TE>=siRNA_threshold) {
    siRNAtop<-c(siRNAtop,"TRUE")
  } else {
    siRNAtop<-c(siRNAtop,"FALSE")
  }
}
# add this column to array
Ptep_top_RNAseq[,9]<-siRNAtop

load("Hmel_top_RNAseq_obj.Rdata")
# make new column for total siRNA abundance from each TE
Hmel_top_RNAseq[,8]<-Hmel_top_RNAseq[,3]+Hmel_top_RNAseq[,4]
# select 15th-highest siRNA abundance value and set as threshold
siRNA_threshold<-sort(Hmel_top_RNAseq[,8],decreasing=T)[15]
# mark whether each TEs is above (TRUE) or below (FALSE) the threshold
siRNAtop<-c()
for (TE in Hmel_top_RNAseq[,8]) {
  if (TE>=siRNA_threshold) {
    siRNAtop<-c(siRNAtop,"TRUE")
  } else {
    siRNAtop<-c(siRNAtop,"FALSE")
  }
}
# add this column to array
Hmel_top_RNAseq[,9]<-siRNAtop

# calculate Spearman Rank correlation between antisense enrichment and siRNA abundance for Hmel & Ptep
HmelCorr<-cor.test(x=(Hmel_top_RNAseq[,6]-Hmel_top_RNAseq[,5]),y=(Hmel_top_RNAseq[,3]+Hmel_top_RNAseq[,4]),method="spearman")
PtepCorr<-cor.test(x=(Ptep_top_RNAseq[,6]-Ptep_top_RNAseq[,5]),y=(Ptep_top_RNAseq[,3]+Ptep_top_RNAseq[,4]),method="spearman")

# plot sense & antisense expression of TEs, and highlight top siRNA targets
size=2
# plot Hmel
# NB: if TE in top 25 RNAi then col 9 ==T
plot(log2(Hmel_top_RNAseq[,5]),log2(Hmel_top_RNAseq[,6]), ylab="antisense",xlab="sense", pch="o",col="red",cex=size)
points(log2(Hmel_top_RNAseq[Hmel_top_RNAseq[,9]==T,5]),log2(Hmel_top_RNAseq[Hmel_top_RNAseq[,9]==T,6]),col="red",pch=19,cex=size)
#plot Ptep
points(log2(Ptep_top_RNAseq[,5]),log2(Ptep_top_RNAseq[,6]),pch="o", col="blue",cex=size)
points(log2(Ptep_top_RNAseq[Ptep_top_RNAseq[,9]==T,5]),log2(Ptep_top_RNAseq[Ptep_top_RNAseq[,9]==T,6]),col="blue",pch=19,cex=size)

#############
## Figure2 ##
#############
library(ggplot2)
# read in data
TEcounts<-read.table("updated_12_08.txt",header=T)
# make new column for difference between antisense and sense (both on log2 scale)
TEcounts[,3]<-log2(TEcounts[,2])-log2(TEcounts[,1])
# make new column for species
TEcounts[,4]<-c("Aaeg","Avul","Bter","Cexi","Dmel","Dvir","Hmel","Lmig","Lpol","Nves","Ofas","Ptep","Pxyl","Tcas","Tni","Smar","Amel","Apis","Mdom","Dvirilis")
# make new column for RdRP presence/absence
TEcounts[,5]<-c("NoRdRP","NoRdRP","NoRdRP","RdRP","NoRdRP","NoRdRP","NoRdRP","NoRdRP","RdRP","NoRdRP","NoRdRP","RdRP","NoRdRP","NoRdRP","NoRdRP","RdRP","NoRdRP","NoRdRP","NoRdRP","NoRdRP")
# make into dataframe and add headers
TEcountsDF<-as.data.frame(TEcounts)
names(TEcountsDF)<-c("antisense","sense","Enrichment","Species","RdRP")
# pre-emptively label Smar for highlighting
highlight.species <- "Smar"
# add a new column to dataframe for highlighting
TEcountsDF$highlight <- ifelse(TEcountsDF$Species == highlight.species, "highlight", "normal")
textdf <- TEcountsDF[TEcountsDF$GeneName == highlight.species, ]
mycolours <- c("highlight" = "red", "normal" = "black")
# plot
RdRPboxplot<-ggplot(data=TEcountsDF,aes(x=RdRP,y=Enrichment))
RdRPboxplot<-RdRPboxplot+theme_bw()
RdRPboxplot<-RdRPboxplot+geom_boxplot(outlier.colour="NA",fill=c("grey"))
RdRPboxplot<-RdRPboxplot+geom_point(size=2,position=position_jitter(width = 0.2),aes(colour=highlight))
RdRPboxplot<-RdRPboxplot+scale_color_manual("Status", values = mycolours)
RdRPboxplot<-RdRPboxplot+xlab("")
RdRPboxplot<-RdRPboxplot+ylab("Antisense enrichment at siRNA -producing loci")
RdRPboxplot<-RdRPboxplot+scale_x_discrete(limits=c("RdRP","NoRdRP"))
RdRPboxplot<-RdRPboxplot+theme(legend.position="none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.border=element_blank(),axis.line.y=element_line(),axis.ticks.x=element_blank())
RdRPboxplot<-RdRPboxplot+annotate("text",x=0.7,y=0,label="Smar",colour="red")
RdRPboxplot
ggsave("Labelled_plot_antisense_enrichment.pdf",plot=RdRPboxplot,device="pdf",width=2.1,height=3.8,units="in")

