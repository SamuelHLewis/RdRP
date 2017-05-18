#############
## Figure1 ##
#############
load("Ptep_top_RNAseq_obj.Rdata")
load("Hmel_top_RNAseq_obj.Rdata")

size=2

plot(log2(Hmel_top_RNAseq[,5]),log2(Hmel_top_RNAseq[,6]), ylab="antisense",xlab="sense", pch=18,col="black",cex=size)
points(log2(Hmel_top_RNAseq[Hmel_top_RNAseq[,7]==T,5]),log2(Hmel_top_RNAseq[Hmel_top_RNAseq[,7]==T,6]), col="red", pch=1,cex=size+1,lwd=size)
#If in top 25 RNAi then col 7 ==T

points(log2(Ptep_top_RNAseq[,5]),log2(Ptep_top_RNAseq[,6]),pch=18, col="gray",cex=size)
points(log2(Ptep_top_RNAseq[Ptep_top_RNAseq[,7]==T,5]),log2(Ptep_top_RNAseq[Ptep_top_RNAseq[,7]==T,6]), col="blue", pch=1,,cex=size+1,lwd=size)

#############
## Figure2 ##
#############
 
load("Logdiff_obj.Rdata")
load("Smar_obj.Rdata")
load("m_obj.Rdata")
load("texts_m_obj.Rdata")
load("texts_nom_obj.Rdata")

#this loads the germline log2(antisense/sense) for all species and the Smar fat body for comparison 
#m is the Rdrp species in the Logdiff object.  
#texts_m and texts_nom are labels for the species.

boxplot(Logdiff[m],Logdiff[-m], col=c("red","green"),boxwex=0.4,names=c("rdrp","no_rdrp"), ylab="log2_antisense_enrichment_siRNAs")
points(rep(2,length=length(Logdiff)-length(m)),Logdiff[-m], pch=18,col="grey",cex=1.4)
for(i in 1:length(texts_m)){text(1.1, Logdiff[m[i]]+0.1,texts_m[i],cex=0.5)}
points(c(1,1,1),Logdiff[m], pch=18,col="black",cex=1.4)
points(2,Smar,col="red", pch=18,cex=1.4)
text(2.1,Smar+0.1,"Smar", col="red", cex=0.5)
dev.copy(pdf, "Labelled_plot_antisense_enrichment.pdf")
dev.off()
 
