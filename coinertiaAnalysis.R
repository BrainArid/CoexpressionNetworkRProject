library("made4")

ciaAndPlot <- function(d1,d2, filename, title, sampleColors, geneColors)
{
  ciaObj <- cia(as.matrix(d1), as.matrix(d2), cia.nf=2, cia.scan=FALSE, nsc=TRUE);
  png(filename=filename);
  plot.cia(ciaObj, nlab = 0, genecol = geneColors,genelabels1 = rownames(Data$ma), genelabels2 = rownames(Data$rs_DESeq), col=sampleColors, clabel=0, title=title)
  dev.off();
  
  plot.cia(ciaObj, nlab = 0, genecol = geneColors,genelabels1 = rownames(Data$ma), genelabels2 = rownames(Data$rs_DESeq), col=sampleColors, clabel=0, title=title)
  rm(ciaObj);
}

Data <- readRDS("../normalization.RDS");
sampleColors<- c(rep("blue",6), rep("orange",91));
geneColors<- "red";#c(rep("blue",6), rep("orange",91));

ciaAndPlot(Data$ma, Data$ma, "cia_ma_ma.png", "Microarray Vs Microarray", sampleColors, geneColors);
ciaAndPlot(Data$ma, Data$rs_DESeq, "cia_ma_DESeq.png", "Microarray Vs DESeq", sampleColors, geneColors);
ciaAndPlot(Data$rs_RPKM, Data$rs_DESeq, "cia_RPKM_DESeq.png", "RPKM Vs DESeq",sampleColors, geneColors);
ciaAndPlot(Data$rs_DESeq, Data$rs_Ubi, "cia_DESeq_Ubi.png", "DESeq Vs Total Ubiquitous genes",sampleColors, geneColors);
ciaAndPlot(Data$rs_DESeq, Data$rs_quant, "cia_DESeq_quant.png", "DESeq Vs quantile",sampleColors, geneColors);
ciaAndPlot(Data$rs_DESeq, Data$rs_raw, "cia_DESeq_raw.png", "DESeq Vs raw",sampleColors, geneColors);