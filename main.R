#Coexpression Network Project main script
#By: Brian Arand
#September 2014

#workingDirectories <- c("/home/barand/microArray_v_RNASeq/","C:/Users/Student/My Research/microArray v RNA Seq/");
#for(wd in workingDirectories)
#  setwd(wd);

print(paste("Current working directory: ", getwd()));
print("Changing Current working directory.")

#setwd("C:/Users/Student/My Research/microArray v RNA Seq/");
setwd("/home/barand/microArray_v_RNASeq/")

print(paste("Current working directory: ", getwd()));

#source("http://bioconductor.org/biocLite.R");

#read in arguments
print("Reading in command line arguments.");
args <- commandArgs(trailingOnly = TRUE);
print(paste("commandArgs: ",args));

if(length(args) > 0)
{
  #Parse arguments (we expec the form --argName=argValue)
  parseArgs <- function (x) strsplit(sub("^--","",x), "=");
  argsDF <- as.data.frame(do.call("rbind", parseArgs(args)));
  args <- as.character(argsDF$V2)
  names(args) <- argsDF$V1
  rm(argsDF);
}
args<- as.list(args);

#initialize arguments if 
initializeBooleanArg <- function(arg, default){
  if(is.null(arg))
  {
    arg <- default;
  }
  else if(is.character(arg))
  {
    arg <- as.logical(arg);
  }
  return(arg);
}

args$diffCoexFlag <- initializeBooleanArg(arg=args$diffCoexFlag, default=TRUE);
args$diffExprsFlag <- initializeBooleanArg(arg=args$diffExprsFlag, default=FALSE);
args$normFlagRPKM <- initializeBooleanArg(arg=args$normFlagRPKM, default=TRUE);
args$normFlagUbi <- initializeBooleanArg(arg=args$normFlagUbi, default=TRUE);
args$normFlagDESeq <- initializeBooleanArg(arg=args$normFlagDESeq, default=TRUE);
args$normFlagQuant <- initializeBooleanArg(arg=args$normFlagQuant, default=TRUE);
args$QCFlag <- initializeBooleanArg(arg=args$QCFlag, default=FALSE);

#import data
print("Importing data files.");
maDir <- "Data/BRCA/Batch 47/Expression-Genes/UNC__AgilentG4502A_07_3/Level_3/";
rsDir <- "Data/BRCA/Batch 47/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/";

#Metadata
metaData <- read.table(file="Data/BRCA/Batch 47/file_manifest.txt",header=TRUE,sep="\t");
metaData <- cbind(metaData, control=(substr(x=metaData[,5],start=14,stop=16)=='11'));#addControl bool
metaData[,6] <- gsub(pattern="-",replacement=".", x=metaData[,6]);#replace '-' with '.' to make mapping easier later

maControlFiles <- paste(sep='', maDir ,as.character(metaData[metaData[,"control"] & metaData[,"Platform.Type"]=="Expression-Genes" ,"File.Name"]));
maCancerFiles <- paste(sep='', maDir ,as.character(metaData[!metaData[,"control"] & metaData[,"Platform.Type"]=="Expression-Genes" ,"File.Name"]));
rsControlFiles <- paste(sep='', rsDir ,as.character(metaData[metaData[,"control"] & metaData[,"Platform.Type"]=="RNASeqV2" & grepl(x=metaData[,"File.Name"], pattern="*.rsem.genes.results") ,"File.Name"]));
rsCancerFiles <- paste(sep='', rsDir ,as.character(metaData[!metaData[,"control"] & metaData[,"Platform.Type"]=="RNASeqV2" & grepl(x=metaData[,"File.Name"], pattern="*.rsem.genes.results") ,"File.Name"]));

source("CoexpressionNetworkRProject/constructCrossSampleFrame.R");
maConData <- constructCrossSampleFrame(inFiles=maControlFiles,rows2Ignore=c(1));
maCanData <- constructCrossSampleFrame(inFiles=maCancerFiles,rows2Ignore=c(1));
rsConData <- constructCrossSampleFrame(inFiles=rsControlFiles,cols2Ignore=c(3,4));
rsCanData <- constructCrossSampleFrame(inFiles=rsCancerFiles,cols2Ignore=c(3,4));
if(args$normFlagRPKM)
{
  rsRPKMData <- constructCrossSampleFrame(inFiles=c(rsControlFiles,rsCancerFiles),cols2Ignore=c(2,3));
}

remove(maControlFiles);
remove(maCancerFiles);
remove(rsControlFiles);
remove(rsCancerFiles);

#take intersection of genes between micro array and rnaseq set
print("Calculating intersection of genes.");
source("CoexpressionNetworkRProject/trim_TCGA_RNASeq_GeneNames.R");
#rsGenes <- trim_TCGA_RNASeq_GeneNames(rsGenes);
rsConData <- trim_TCGA_RNASeq_GeneNames(rsConData);
rsCanData <- trim_TCGA_RNASeq_GeneNames(rsCanData);
#rsRPKMData <- trim_TCGA_RNASeq_GeneNames(rsRPKMData);

sharedGenes <- intersect(row.names(maConData), row.names(rsConData));
maConData <- maConData[sharedGenes,];
rsConData <- rsConData[sharedGenes,];
maCanData <- maCanData[sharedGenes,];
rsCanData <- rsCanData[sharedGenes,];
#rsRPKMData <- rsRPKMData[sharedGenes,];
remove(sharedGenes);
conCount <- dim(maConData)[2];
canCount <- dim(maCanData)[2];

countData <- mapply(x=as.matrix(c(rsConData, rsCanData)), FUN=as.integer);
row.names(countData) <- row.names(rsConData);

#normalize
print("Begin normalization:")

#RMKM
if(args$normFlagRPKM)
{
  print("RPKM normalization:")
  source("http://bioconductor.org/biocLite.R")
  biocLite("easyRNASeq")
  library("easyRNASeq")
  
  biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
  biocLite("org.Hs.eg.db")
  hg19GeneLengths <- function(symbols)
  {
    require(TxDb.Hsapiens.UCSC.hg19.knownGene) 
    require(org.Hs.eg.db)
    exons.db = exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by='gene')    
    egs    = unlist(  mget(symbols[ symbols %in% keys(org.Hs.egSYMBOL2EG) ],org.Hs.egSYMBOL2EG) )
    sapply(egs,function(eg)
    {
      exons = exons.db[[eg]]
      if(is.null(exons))
      {return (0)}
      exons = reduce(exons)
      sum( width(exons) )
    })
    return as.numberic(egs);
  }
  rpkm <- function(expMatrix)
  {
    #RPKM = (10^9 * C)/(N * L), with
    #C = Number of reads mapped to a gene
    #N = Total mapped reads in the experiment
    #L = gene length in base-pairs for a gene
    lengths <- hg19GeneLengths(row.names(expMatrix));
    I<- intersect(names(lengths),row.names(expMatrix));
    expMatrix <- expMatrix[I,];
    lengths <- lengths[I];
    rpkmMat <- apply(X=expMatrix,MARGIN=2,FUN=function(x){return(x/sum(x))})
    rpkmMat <- 10^9 * rpkmMat 
    rpkmMat <- rpkmMat / as.numeric(lengths);
    return(rpkmMat);
  }
  
  rsData_RPKM <- rpkm(countData);
}

#DESeq variance stableizing transformation (VST) normalization
if(args$normFlagDESeq)
{
  print("DESeq normalization:");
  source("http://bioconductor.org/biocLite.R");
  biocLite("DESeq2");
  library("DESeq2");

  colData <- DataFrame(condition=c(rep(x="control",times=conCount), rep(x="cancer", times=canCount)), type=c(rep(x="single-read",times=conCount+canCount)));
  #row.names(colData) <- colnames(countData);
  dds <- DESeqDataSetFromMatrix(countData=countData, design = ~ condition,colData=colData);
  dds <- DESeq(dds);
  dds$condition <- factor(dds$condition, levels=c("cancer","control"));
  res <- results(dds);
  res$log2FoldChange<- -res$log2FoldChange;#reverse direction of test
  
  if(args$diffExprsFlaq)
  {
    #DESeq differential Expression
    resOrdered <- res[order(res$padj),];  
    head(resOrdered);
    n <- 150;
    write.csv(resOrdered,file=paste("DiffExpression DesSEQ.csv"),quote=FALSE,);
    rm(resOrdered);
  }
    normalizedCounts_DESeq2 <- t( t(counts(dds)) / sizeFactors(dds) );
    rsConData_norm_DESeq2 <- normalizedCounts_DESeq2[,1:conCount];
    rsCanData_norm_DESeq2 <- normalizedCounts_DESeq2[,conCount+1:canCount];
    row.names(rsConData_norm_DESeq2) <- row.names(normalizedCounts_DESeq2);
    row.names(rsCanData_norm_DESeq2) <- row.names(normalizedCounts_DESeq2);
    remove(normalizedCounts_DESeq2);
}

#EdgeR does not do cross sample normalization
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#library("edgeR")

#total ubiquetous normalization as per "Optimal scaling of Digital Transcriptomes" By Glusman et al.
if(args$normFlagUbi)
{
  print("Total ubiquitous normalization:");
  #setwd("CoexpressionNetworkProject/")
  source("CoexpressionNetworkRProject/ubiquitousNormalize.R");
  source("CoexpressionNetworkRProject/getUbiquitousGeneSet.R");
  source("CoexpressionNetworkRProject/getTrimmedSet_UbiquitousHelper.R");
  rsData_Ubi <- ubiquitousNormalize(countData,lowerPercentile=0.3,upperPercentile=0.85);
  #setwd("..")
}

#quantile
if(args$normFlagQuant)
{
  print("Quantile normalization:");
  source('http://bioconductor.org/biocLite.R');
  biocLite('preprocessCore');
  #load package
  library(preprocessCore);
  
  #quantile differential expression
  #...
  normalizedCounts_quant <- normalize.quantiles(x= countData, copy=TRUE);
  rsConData_norm_quant <- normalizedCounts_quant[,1:conCount];
  rsCanData_norm_quant <- normalizedCounts_quant[,conCount+1:canCount];
  remove(normalizedCounts_quant);
}

if(args$diffExprsFlaq)
{
  print("Gene Differential Expression:");
  
  #microarray differential expression
  #may not be working... produces rediculously small p-values
  source("http://bioconductor.org/biocLite.R");
  biocLite( "limma");
  library("limma");
  
  condition=c(rep(x="control",times=dim(x=rsConData)[2]), rep(x="cancer", times=canCount));
  combn <- factor(paste(pData(phenoData)[,1], pData(phenoData)[,2], sep = "_"));
  design <- model.matrix(~condition);# describe model to be fit
  
  fit <- lmFit(cbind(maConData, maCanData), design);# fit each probeset to model
  efit <- eBayes(fit);# empirical Bayes adjustment
  write.csv(topTable(efit, coef=2, number=length(efit$p.value)), file=paste("DiffExpression MicroArray.csv"),quote=FALSE);

  #significance comparison
  maPRank <- rank(efit$p.value[,2]);
  rsPRank <- rank(res$padj);
  names(rsPRank) <- row.names(res);
  maPRank <- maPRank[sort(names(maPRank))];
  rsPRank <- rsPRank[sort(names(rsPRank))];
  plot2Groups(GroupA=maPRank, GroupB=rsPRank, xlab="MicroArray Rank", ylab="RNASeq Rank", main="DEG Rank comparison by p-value", file="Comp_DEG_pVal_rank_across_tech.png");
  Rs_PRank <- cor(x=maPRank, y=rsPRank, method="spearman");
  rm(maPRank);
  rm(rsPRank);
  
  #log foldchange comparison
  maFC <- topTable(efit, coef=2, number=length(efit$p.value))[,1];
  names(maFC) <- row.names(topTable(efit, coef=2, number=length(efit$p.value)));
  maFCRank <- rank(maFC);
  rsFC <- res$log2FoldChange;
  names(rsFC) <- row.names(res);
  rsFCRank <- rank(rsFC);
  maFCRank <- maFCRank[sort(names(maFCRank))];
  rsFCRank <- rsFCRank[sort(names(rsFCRank))];
  plot2Groups(GroupA=maFCRank, GroupB=rsFCRank, xlab="MicroArray Rank", ylab="RNASeq Rank", main="Log2 Fold Change Rank comparison", file="Comp_FC_rank_across_tech.png");
  Rs_FCRank <- cor(x=maFCRank, y=rsFCRank, method="spearman");
  rm(maFCRank);
  rm(rsFCRank);
}

#subset by significance
if(args$diffExprsFlag)
{
  subBySig <- FALSE;
  if(subBySig)
  {
    print("Filtering genes by significance:");
    cutoff<- 0.0000001;
    filter <- efit.p.adj<cutoff & !is.na(efit.p.adj);
    maGenes<- as.matrix(cbind(efit.p.adj[filter], maConData[filter,],maCanData[filter,]));#attach adjusted p value now
    filter <- res$padj<cutoff & !is.na(res$padj);
    rsGenes<- as.matrix(cbind(res$padj[filter], rsConData_norm_DESeq2[filter,],rsCanData_norm_DESeq2[filter,]));#attach adjusted p value now
  }
  else #subset by rank top X most significant
  {
    print("Filtering genes by fold-change:")
    cutoff <- 12000;
    filter <- rank(efit.p.adj)<=cutoff & !is.na(efit.p.adj);
    maGenes<- as.matrix(cbind(efit.p.adj[filter], maConData[filter,],maCanData[filter,]));#attach adjusted p value now
    filter <- rank(res$padj)<=cutoff & !is.na(res$padj);
    rsGenes<- as.matrix(cbind(res$padj[filter], rsConData_norm_DESeq2[filter,],rsCanData_norm_DESeq2[filter,]));#attach adjusted p value now
  }
  
  #intersection of technology significant genes and output venn diagram gene lists
  filter <- intersect(row.names(maGenes), row.names(rsGenes));
  maIntGenes <- maGenes[filter,];
  rsIntGenes <- rsGenes[filter,];
  maUniGenes <- maGenes[setdiff(row.names(maGenes),row.names(maIntGenes)),];
  rsUniGenes <- rsGenes[setdiff(row.names(rsGenes),row.names(rsIntGenes)),];
  
  #top genes in each group
  head(sort(maIntGenes[,1]))
  head(sort(rsIntGenes[,1]))
  head(sort(maUniGenes[,1]))
  head(sort(rsUniGenes[,1]))
  #remove p.values from Gene lists
  maGenes <- maGenes[,-1];
  rsGenes <- rsGenes[,-1];
}
else
{
  print("No gene filtering performed.")
  maGenes <- cbind(maConData, maCanData);
  rsGenes <- cbind(rsConData_norm_DESeq2,rsCanData_norm_DESeq2);
  cutoff <- "allGenes";
}

if(args$QCFlag)
{
  print("Outputing quality control figures:");
  
  maRankData<-apply(maGenes,MARGIN=2,FUN=rank);
  rsRankData<-apply(rsGenes,MARGIN=2,FUN=rank);
  
  #plot means
  source("CoexpressionNetworkProject/plot2Groups.R");
  maGenes<- cbind(maConData, maCanData, maFC, maP);
  rsData_DESeq2 <- cbind(rsConData_norm_DESeq2, rsCanData_norm_DESeq2);
  rsData_quant <- cbind(rsConData_norm_quant, rsCanData_norm_quant);
  
  plot2Groups(rowMeans(maGenes, na.rm = TRUE), log(rowMeans(rsData_DESeq2, na.rm = TRUE)),main="Micro Array vs RNASeq DESeq gene means (91 paired patient samples)",xlab="Lowess Normalized MicroArray",ylab="Log RNASeq DEseq counts", file="Comp_gene_means_across_tech_DESeq.png", histA=TRUE, histB=TRUE);
  plot2Groups(apply(maGenes,1,median, na.rm = TRUE), apply(log(rsData_DESeq2),1,median, na.rm = TRUE),main="Micro Array vs RNASeq DESeq gene medians (91 paired patient samples)",xlab="Lowess Normalized MicroArray",ylab="Log RNASeq DESeq counts", file="Comp_gene_medians_across_tech_DESeq.png", histA=TRUE, histB=TRUE);
  plot2Groups(rowMeans(apply(maGenes,MARGIN=2,FUN=rank), na.rm = TRUE), rowMeans(apply(rsData_DESeq2,MARGIN=2,FUN=rank), na.rm = TRUE),main="Micro Array vs RNASeq DESeq gene mean ranks (91 paired patient samples)",xlab="Ranked Lowess Normalized MicroArray",ylab="Log RNASeq DESeq counts", file="Ranked Comp_gene_mean_rank_across_tech_DESeq.png");
  plot2Groups(apply(apply(maGenes,MARGIN=2,FUN=rank), 1,median, na.rm = TRUE), apply(apply(rsData_DESeq2,MARGIN=2,FUN=rank),1,median, na.rm = TRUE),main="Micro Array vs RNASeq DESeq gene median ranks (91 paired patient samples)",xlab="Ranked Lowess Normalized MicroArray",ylab="Ranked Log RNASeq DESeq counts", file="Comp_gene_median_rank_across_tech_DESeq.png");
  
  plot2Groups(rowMeans(maGenes, na.rm = TRUE), log(rowMeans(rsData_quant, na.rm = TRUE)),main="Micro Array vs RNASeq quant gene means (91 paired patient samples)",xlab="Lowess Normalized MicroArray",ylab="Log RNASeq quant counts", file="Comp_gene_means_across_tech_quant.png", histA=TRUE, histB=TRUE);
  
  GeneID=3245;
  GeneName= row.names(maGenes)[GeneID];
  plot2Groups(t(maGenes[GeneID,]), t(log(rsData_DESeq2[GeneID,])),main=paste("Micro Array vs RNASeq DESeq (", GeneName, ")"),xlab="Lowess Normalized MicroArray",ylab="Log RNASeq DESeq counts", file=paste("Comp_gene_across_tech_DESeq-",GeneName,".png"));
  #plot2Groups(t(maRankData[GeneID,]), t(rsRankData[GeneID,]),main=paste("Micro Array vs RNASeq ranks DESeq (", GeneName, ")"),xlab="Ranked Lowess Normalized MicroArray",ylab="Ranked Log RNASeq DESeq counts", file=paste("Comp_gene_rank_across_tech_DESeq-",GeneName,".png"));
  plot2Groups(t(maGenes[GeneID,]), t(log(rsData_quant[GeneID,])),main=paste("Micro Array vs RNASeq quant (", GeneName, ")"),xlab="Lowess Normalized MicroArray",ylab="Log RNASeq quant counts", file=paste("Comp_gene_across_tech_quant-",GeneName,".png"));
  
  #check normalizations with boxplots
  png(filename="MicroArray normalization check.png");
  boxplot(x=maGenes,names=seq(1,dim(maGenes)[2]), outcex=0.5, outpch=20, main="Patient box plots", xlab="Patient", ylab="expression value");
  dev.off();
  png(filename="RNASeq Count normalization check.png");
  boxplot(x=log(rsGenes),names=seq(1,dim(rsGenes)[2]), outcex=0.5, outpch=20, main="Patient box plots", xlab="Patient", ylab="log RNA seq counts");
  dev.off();
  png(filename="RNASeq RPKM normalization check.png");
  boxplot(x=log(rsRPKMData),names=seq(1,dim(rsRPKMData)[2]), outcex=0.5, outpch=20, main="Patient box plots", xlab="Patient", ylab="log RNA seq RPKM");
  dev.off();
  png(filename="RNASeq DESeq normalization check.png");
  boxplot(x=log(cbind(rsConData_norm_DESeq2, rsCanData_norm_DESeq2)),names=seq(1,dim(rsConData_norm_DESeq2)[2]+dim(rsCanData_norm_DESeq2)[2]), outcex=0.5, outpch=20, main="Patient box plots", xlab="Patient", ylab="log RNA seq DESeq normalized");
  dev.off();
  png(filename="RNASeq quantile normalization check.png");
  boxplot(x=log(cbind(rsConData_norm_quant, rsCanData_norm_quant)),names=seq(1,dim(rsConData_norm_quant)[2]+dim(rsCanData_norm_quant)[2]), outcex=0.5, outpch=20, main="Patient box plots", xlab="Patient", ylab="log RNA seq quan normalized");
  dev.off();
  png(filename="RNASeq total ubiquitous normalization check.png");
  boxplot(x=log(rsData_Ubi),names=seq(1,dim(rsConData_norm_quant)[2]+dim(rsCanData_norm_quant)[2]), outcex=0.5, outpch=20, main="Patient box plots", xlab="Patient", ylab="log RNA seq total ubiquitous normalized");
  dev.off();
}

#calculated correlation statistics
print("Constructing correlation matricies");
print("Calculating Microarray correlation matrix:");
maPearson <- cor(x=t(maGenes), method="pearson", use="complete.obs");
maPearson_hist <- hist(x=maPearson,breaks=100,plot=FALSE);
write.csv(x=maPearson,file=paste("Data/BRCA/maPearson_",cutoff,"_int.txt"));
if(!args$diffCoexFlag)
{
  rm(maPearson);
}

print("Calculating RNASeq DESeq correlation matrix:");
rsPearson_DESeq <- cor(x=t(rsGenes), method="pearson", use="complete.obs");
rsPearson_DESeq_hist <- hist(x=rsPearson_DESeq,breaks=100,plot=FALSE);
write.csv(x=rsPearson_DESeq,file=paste("Data/BRCA/rsPearson_",cutoff,"_int.txt"));
if(!args$diffCoexFlag)
{
  rm(rsPearson_DESeq);
}

print("Calculating RNASeq Quantile correlation matrix:");
rsPearson_quant <- cor(x=t(rsData_quant), method="pearson", use="complete.obs");
rsPearson_quant_hist <- hist(x=rsPearson_quant,breaks=100,plot=FALSE);
rm(rsPearson_quant);

print("Calculating RNASeq Total Ubiquitous correlation matrix:");
rsPearson_Ubi <- cor(x=t(rsData_Ubi), method="pearson", use="complete.obs");
rsPearson_Ubi_hist <- hist(x=rsPearson_Ubi,breaks=100,plot=FALSE);
write.csv(x=rsPearson_Ubi,file=paste("Data/BRCA/rsPearson_top12000_int.txt"));
rm(rsPearson_Ubi);

print("Calculating RNASeq Raw correlation matrix:");
rsPearson_raw <- cor(x=t(countData), method="pearson", use="complete.obs");
rsPearson_raw_hist <- hist(x=rsPearson_raw,breaks=100,plot=FALSE);
rm(rsPearson_raw);

print("Calculating RNASeq RPKM correlation matrix:");
rsPearson_RPKM <- cor(x=t(rsData_RPKM), method="pearson", use="complete.obs");
rsPearson_RPKM_hist <- hist(x=rsPearson_RPKM,breaks=100,plot=FALSE);
rm(rsPearson_RPKM);

#plot overlapping histogram of PCC
print("Outputting comparative PCC histogram:");
maxy <- 0.08
# Density plots
library("ggplot2")

PCC_density_data <- data.frame(
  PCC=c(
    maPearson_hist$mids, 
    rsPearson_DESeq_hist$mids, 
    rsPearson_quant_hist$mids, 
    rsPearson_Ubi_hist$mids, 
    rsPearson_raw_hist$mids, 
    rsPearson_RPKM_hist$mids), 
  density=c(
    maPearson_hist$counts/sum(maPearson_hist$counts),
    rsPearson_DESeq_hist$counts/sum(rsPearson_DESeq_hist$counts),
    rsPearson_quant_hist$counts/sum(rsPearson_quant_hist$counts),
    rsPearson_Ubi_hist$counts/sum(rsPearson_Ubi_hist$counts),
    rsPearson_raw_hist$counts/sum(rsPearson_raw_hist$counts),
    rsPearson_RPKM_hist$counts/sum(rsPearson_RPKM_hist$counts)),
  method=c(
    rep(x="microArray",times=length(maPearson_hist$counts)),
    rep(x="DESeq",times=length(rsPearson_DESeq_hist$counts)),
    rep(x="quantile",times=length(rsPearson_quant_hist$counts)),
    rep(x="Total Ubiquitous",times=length(rsPearson_Ubi_hist$counts)),
    rep(x="raw",times=length(rsPearson_raw_hist$counts)),
    rep(x="RPKM",times=length(rsPearson_RPKM_hist$counts))));

png(filename="RNASeq total ubiquitous normalization check.png");
ggplot(data=PCC_density_data, aes(x=PCC, y=density, group=method, colour=method)) + 
  geom_line(size=1, aes(linetype=method)) +
  ggtitle("PCC density comparison");
dev.off();

rm(PCC_density_data);

if(args$diffCoexFlaq)
{
  print("Calculating differential coexpression network.");
  diffPearson <- rsPearson - maPearson
  write.csv(x=diffPearson,file=paste("Data/BRCA/diffPearson.txt"));

  #output venn diagram gene-edge lists
}

#create iGraph and plot
library("igraph")

print("Outputting microArray coexpression matrix");
write.csv(x=maPearson,file=paste("Data/BRCA/maPearson_intersection_", cutoff, ".txt"));
range(maPearson);
png(filename=paste("Data/BRCA/maPearsonHist_intersection_", cutoff, ".png"));
hist(x=maPearson,breaks=20);
dev.off();

print("Constructing microArray iGraph.");
maGraph <-graph.adjacency(adjmatrix=maPearson*1000,mode="undirected", weighted=TRUE);
#add vertex attributes to graph
for(i in 1:length(V(maGraph)))
{
  name <- V(maGraph)[i]$name;
  V(maGraph)[i]$fc <- maFC[name];
  V(maGraph)[i]$p <- efit.p.adj[name];
}

print("Outputting microArray iGraph.");
write.graph(maGraph, file=paste("maGraph_", cutoff, ".graphml"), format="graphml" );
rm(maGraph);
rm(maPearson);

print("Outputting RNASeq coexpression matrix");
write.table(x=rsPearson,file=paste("Data/BRCA/rsPearson_intersection_", cutoff, ".txt"));
range(rsPearson);
png(filename=paste("Data/BRCA/rsPearsonHist_intersection_", cutoff, ".png"));
hist(x=rsPearson,breaks=20);
dev.off();

print("Constructing RNASeq iGraph.");
rsGraph <-graph.adjacency(adjmatrix=rsPearson*1000,mode="undirected", weighted=TRUE);
for(i in 1:length(V(rsGraph)))
{
  name <- V(rsGraph)[i]$name;
  V(rsGraph)[i]$fc <- rsFC[name];
  V(rsGraph)[i]$p <- res[name,]$padj;
}

print("Outputting RNASeq iGraph");
write.graph(rsGraph, file=paste("rsGraph_", cutoff, ".graphml"), format="graphml" );
rm(rsGraph);
rm(rsPearson);

if(args$diffCoexFlag)
{
  print("Outputting differential coexpression matrix");
  write.csv(x=diffPearson,file=paste("Data/BRCA/diffPearson_intersection_", cutoff, ".txt"));
  range(diffPearson);
  png(filename=paste("Data/BRCA/diffPearsonHist_intersection_", cutoff, ".png"));
  hist(x=maPearson,breaks=20);
  dev.off();
  
  print("Constructing differential coexpression iGraph");
  diffGraph <-graph.adjacency(adjmatrix=diffPearson*1000,mode="undirected", weighted=TRUE);
  for(i in 1:length(V(diffGraph)))
  {
    name <- V(diffGraph)[i]$name;
    V(diffGraph)[i]$RS_fc <- rsFC[name];
    V(diffGraph)[i]$RS_p <- res[name,]$padj;
    V(diffGraph)[i]$MA_fc <- maFC[name];
    V(diffGraph)[i]$MA_p <- efit.p.adj[name];
  }
  print("Outputting differential coexpression iGraph");
  write.graph(diffGraph, file=paste("diffGraph_", cutoff, ".graphml"), format="graphml" );
  rm(diffGraph);
  rm(diffPearson);
}
#calculate difference network and fold-change network
#difNet <- rsPearson - maPearson;
#fcNet <- rsPearson / maPearson;
#visualize network

##scale free networks via WGCNA
##Load WGCNA package
#library(WGCNA);
##Load additional necessary packages
#library(cluster);
#k=softConnectivity(datE=t(maGenes),power=6);
## Plot a histogram of k and a scale free topology plot
#sizeGrWindow(10,5);
#par(mfrow=c(1,2));
#png(filename="Data/BRCA/maPearson_WGCNA-power6_Hist.png");
#hist(k, main="Connectivity (MArray Pearson pow6)");
#dev.off();
#png(filename="Data/BRCA/maPearson_WGCNA-power6_ScaleFreePlot.png");
#scaleFreePlot(k, main="Check scale free topology (MArray Pearson pow6\n");
#dev.off();

#k=softConnectivity(datE=t(rsGenes),power=6);
## Plot a histogram of k and a scale free topology plot
#sizeGrWindow(10,5);
#par(mfrow=c(1,2));
#png(filename="Data/BRCA/rsPearson_WGCNA-power6_Hist.png");
#hist(k, main="Connectivity (RNASeq Pearson pow6)");
#dev.off();
#png(filename="Data/BRCA/rsPearson_WGCNA-power6_ScaleFreePlot.png");
#scaleFreePlot(k, main="Check scale free topology (RNASeq Pearson pow6)\n");
#dev.off();

quit();
