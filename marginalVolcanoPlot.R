CON <- 1:6
CAN <- 7:97
source("CoexpressionNetworkRProject/plot2Groups.R");

#DESeq
head(res[order(res$padj),])
name <- "COL10A1"
Data$rs_DESeq[name,]
con_Avg <- mean(Data$rs_DESeq[name,CON])
can_Avg <- mean(Data$rs_DESeq[name,CAN])
con_Std <- sqrt(var(Data$rs_DESeq[name,CON]))
can_Std <- sqrt(var(Data$rs_DESeq[name,CAN]))
log2(can_Avg/con_Avg)
res[name,]$log2FoldChange

#volcano plot
plot(y=-log(res$pvalue,base=10),x=-res$log2FoldChange,ylab="-log10(p-value)",xlab="Fold-change",main="RNASeq (DESeq) volcano")

#this doesn't look right
hist(res$pvalue[res$pvalue<1],breaks=500)
hist(res$pvalue[res$pvalue<0.001],breaks=500)
hist(res$pvalue[res$pvalue<0.0001],breaks=500)
hist(res$pvalue[res$pvalue<0.00001],breaks=500)
hist(res$pvalue[res$pvalue<0.000001],breaks=500)
hist(res$pvalue[res$pvalue<0.0000001],breaks=500)
hist(res$pvalue[res$pvalue<0.000000001],breaks=500)
hist(res$pvalue[res$pvalue<0.00000000001],breaks=500)
hist(res$pvalue[res$pvalue<0.000000000001],breaks=500)

hist(log(-res$log2FoldChange,base=10),breaks=500)

#limma
head(efit$p.value[order(efit$p.value[,2]),])
Data$ma[name,]
con_Avg <- mean(unlist(Data$ma[name,CON]))
can_Avg <- mean(unlist(Data$ma[name,CAN]))
con_Std <- sqrt(var(unlist(Data$ma[name,CON])))
can_Std <- sqrt(var(unlist(Data$ma[name,CAN])))
log2(can_Avg/con_Avg)

topTable(efit, coef=2, number=length(efit$p.value))[name,1];

#volcano plot
plot(y=-log(topTable(efit, coef=2, number=length(efit$p.value))[,5],base=10),x=topTable(efit, coef=2, number=length(efit$p.value))[,1],ylab="-log10(p-value)",xlab="Fold-change",main="Microarray (limma) volcano")

#this doesn't look right
hist <- hist(-log(topTable(efit, coef=2, number=length(efit$p.value))[topTable(efit, coef=2, number=length(efit$p.value))[,5]<1,5]),breaks=500,plot=FALSE)
plot(hist,ylab="Frequency(Count)",xlab="-log(p-value)",main="Microarray (limma) -log(p-value) distribution")

hist <- hist(log(topTable(efit, coef=2, number=length(efit$p.value))[,1],base=10),breaks=120)
plot(hist,ylab="Frequency(Count)",xlab="log10(Fold-change)",main="Microarray (limma) log10(Fold-change) distribution")

plot2Groups(GroupA=topTable(efit, coef=2, number=length(efit$p.value))[,1],GroupB=-log(topTable(efit, coef=2, number=length(efit$p.value))[,5],base=10),byRow=TRUE,ylab="-log(p-value) ",xlab="Fold-change",main="Microarray (limma) volcano",file="Microarray Volcano.png",histA=TRUE,histB=TRUE,breaksA=120,breaksB=120)

marginalVolcanoPlots <- function(pVals, fcVals, mainPrefix)
{
  #this doesn't look right
  png(filename=paste0(mainPrefix,"pValDist.png"));
  hist <- hist(-log(pVals,base=10),breaks=500,plot=FALSE)
  plot(hist,ylab="Frequency(Count)",xlab="-log10(p-value)",main=paste0(mainPrefix, " -log10(p-value) distribution"))
  dev.off();
  
  png(filename=paste0(mainPrefix,"fcDist.png"));
  hist <- hist(log(fcVals,base=10),breaks=120)
  plot(hist,ylab="Frequency(Count)",xlab="log10(Fold-change)",main=paste0(mainPrefix," log10(Fold-change) distribution"))
  lines(x=rep(x=mean(log10(x=fcVals)),times=2),y=c(0,max(hist$counts)),col="Red")
  
  dev.off();
  
  plot2Groups(GroupA=log(fcVals,base=10),GroupB=-log(pVals,base=10),byRow=TRUE,ylab="-log(p-value) ",xlab="log(Fold-change)",main=paste0(mainPrefix," volcano"),file=paste0(mainPrefix,"_Volcano.png"),histA=TRUE,histB=TRUE,breaksA=120,breaksB=120)
}

source("CoexpressionNetworkRProject/plot2Groups.R");

pVals <- topTable(efit, coef=2, number=length(efit$p.value))[,5];
fcVals <- 2^topTable(efit, coef=2, number=length(efit$p.value))[,1];
mainPrefix <- "Microarray (limma)";
marginalVolcanoPlots(pVals, fcVals, mainPrefix);

pVals <- res[!is.na(res$padj),]$padj
fcVals <- 2^res[!is.na(res$padj),]$log2FoldChange
mainPrefix <- "RNASeq (DESeq)";
marginalVolcanoPlots(pVals, fcVals, mainPrefix);