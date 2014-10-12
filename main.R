#Coexpression Network Project main script
#By: Brian Arand
#September 2014

#workingDirectories <- c("/home/barand/microArray_v_RNASeq/","C:/Users/Student/My Research/microArray v RNA Seq/");
#for(wd in workingDirectories)
#  setwd(wd);
print("Library Paths:");
.libPaths();

print(paste("Current working directory: ", getwd()));
print("Changing Current working directory.")

#setwd("C:/Users/Student/My Research/microArray v RNA Seq/");
setwd("/home/barand/microArray_v_RNASeq/")

print(paste("Current working directory: ", getwd()));

source("http://bioconductor.org/biocLite.R");

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
args$dataFromRDS <- initializeBooleanArg(arg=args$dataFromRDS, default=FALSE);
args$saveNormalizationRDS <- initializeBooleanArg(arg=args$saveNormalizationRDS, default=FALSE);

#import data
print("Importing data files.");
maDir <- "Data/BRCA/Batch 47/Expression-Genes/UNC__AgilentG4502A_07_3/Level_3/";
rsDir <- "Data/BRCA/Batch 47/RNASeqV2/UNC__IlluminaHiSeq_RNASeqV2/Level_3/";

#Metadata
metaData <- read.table(file="Data/BRCA/Batch 47/file_manifest.txt",header=TRUE,sep="\t");
metaData <- cbind(metaData, control=(substr(x=metaData[,5],start=14,stop=16)=='11'));#addControl bool
metaData[,6] <- gsub(pattern="-",replacement=".", x=metaData[,6]);#replace '-' with '.' to make mapping easier later

if(args$dataFromRDS)
{
  Data <- readRDS("normalization.RDS")
}
else
{
  maControlFiles <- paste(sep='', maDir ,as.character(metaData[metaData[,"control"] & metaData[,"Platform.Type"]=="Expression-Genes" ,"File.Name"]));
  maCancerFiles <- paste(sep='', maDir ,as.character(metaData[!metaData[,"control"] & metaData[,"Platform.Type"]=="Expression-Genes" ,"File.Name"]));
  rsControlFiles <- paste(sep='', rsDir ,as.character(metaData[metaData[,"control"] & metaData[,"Platform.Type"]=="RNASeqV2" & grepl(x=metaData[,"File.Name"], pattern="*.rsem.genes.results") ,"File.Name"]));
  rsCancerFiles <- paste(sep='', rsDir ,as.character(metaData[!metaData[,"control"] & metaData[,"Platform.Type"]=="RNASeqV2" & grepl(x=metaData[,"File.Name"], pattern="*.rsem.genes.results") ,"File.Name"]));
  
  source("CoexpressionNetworkRProject/constructCrossSampleFrame.R");
  Data <- list();
  Data$ma_con <- constructCrossSampleFrame(inFiles=maControlFiles,rows2Ignore=c(1));
  Data$ma_can <- constructCrossSampleFrame(inFiles=maCancerFiles,rows2Ignore=c(1));
  Data$rs_con <- constructCrossSampleFrame(inFiles=rsControlFiles,cols2Ignore=c(3,4));
  Data$rs_can <- constructCrossSampleFrame(inFiles=rsCancerFiles,cols2Ignore=c(3,4));
  
  if(is.null(Data$ma))
  {
    Data$ma <- cbind(Data$ma_con, Data$ma_can);
  }
  if(is.null(Data$rs_raw))
  {
    Data$rs_raw <- cbind(Data$rs_con, Data$rs_can);
    source("CoexpressionNetworkRProject/trim_TCGA_RNASeq_GeneNames.R");
    Data$rs_raw <- trim_TCGA_RNASeq_GeneNames(Data$rs_raw);
  }
  Data$conCount <- dim(Data$ma_con)[2];
  Data$canCount <- dim(Data$ma_can)[2];
  
  Data$ma_con <- NULL;
  Data$ma_can <- NULL;
  Data$rs_con <- NULL;
  Data$rs_can <- NULL;
  
  rm(maControlFiles, maCancerFiles, rsControlFiles, rsCancerFiles);
  
  #take intersection of genes between micro array and rnaseq set
  print("Calculating intersection of genes.");
  sharedGenes <- intersect(row.names(Data$ma), row.names(Data$rs_raw));
  Data$ma <- Data$ma[sharedGenes,];
  Data$rs_raw <- Data$rs_raw[sharedGenes,];
  remove(sharedGenes);
  
  Data$rs_raw <- matrix(data=mapply(x=as.matrix(Data$rs_raw), FUN=as.integer),nrow = dim(Data$rs_raw)[1],ncol=dim(Data$rs_raw)[2],dimnames = list(row.names(Data$rs_raw), colnames(Data$rs_raw)));
}
#normalize
print("Begin normalization:")

#RMKM
if(args$normFlagRPKM && is.null(Data$rs_RPKM))
{
  print("RPKM normalization:")
  #source("http://bioconductor.org/biocLite.R")
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
    return (as.numeric(egs));
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
  
  Data$rs_RPKM <- rpkm(Data$rs_raw);
}

#DESeq variance stableizing transformation (VST) normalization
if(args$normFlagDESeq && is.null(Data$rs_DESeq))
{
  print("DESeq normalization:");
  #source("http://bioconductor.org/biocLite.R");
  biocLite("DESeq2");
  library("DESeq2");

  colData <- DataFrame(condition=c(rep(x="control",times=Data$conCount), rep(x="cancer", times=Data$canCount)), type=c(rep(x="single-read",times=Data$conCount+Data$canCount)));
  #row.names(colData) <- colnames(countData);
  dds <- DESeqDataSetFromMatrix(countData=matrix(data=mapply(x=as.matrix(Data$rs_raw), FUN=as.integer),nrow = dim(Data$rs_raw)[1],ncol=dim(Data$rs_raw)[2],dimnames = list(row.names(Data$rs_raw), colnames(Data$rs_raw))), design = ~ condition,colData=colData);
  dds <- DESeq(dds);
  dds$condition <- factor(dds$condition, levels=c("cancer","control"));
  res <- results(dds);
  res$log2FoldChange<- -res$log2FoldChange;#reverse direction of test
  
  if(args$diffExprsFlag)
  {
    #DESeq differential Expression
    resOrdered <- res[order(res$padj),];  
    head(resOrdered);
    n <- 150;
    write.csv(resOrdered,file=paste("DiffExpression DesSEQ.csv"),quote=FALSE,);
    rm(resOrdered);
  }
  Data$rs_DESeq <- t( t(counts(dds)) / sizeFactors(dds) );
}

#EdgeR does not do cross sample normalization
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#library("edgeR")

#total ubiquetous normalization as per "Optimal scaling of Digital Transcriptomes" By Glusman et al.
if(args$normFlagUbi && is.null(Data$rs_Ubi))
{
  print("Total ubiquitous normalization:");
  #setwd("CoexpressionNetworkProject/")
  source("CoexpressionNetworkRProject/ubiquitousNormalize.R");
  source("CoexpressionNetworkRProject/getUbiquitousGeneSet.R");
  source("CoexpressionNetworkRProject/getTrimmedSet_UbiquitousHelper.R");
  Data$rs_Ubi <- ubiquitousNormalize(Data$rs_raw,lowerPercentile=0.3,upperPercentile=0.85);
  #setwd("..")
}

#quantile
if(args$normFlagQuant && Data$rs_quant == NULL)
{
  print("Quantile normalization:");
  source('http://bioconductor.org/biocLite.R');
  biocLite('preprocessCore');
  #load package
  library(preprocessCore);
  
  #quantile differential expression
  #...
  Data$rs_quant <- normalize.quantiles(x= Data$rs_raw, copy=TRUE);
}

if(args$saveNormalizationRDS)
{
  saveRDS(data, file="normalizationData.rds")
}

if(args$diffExprsFlag)
{
  print("Gene Differential Expression:");
  
  #microarray differential expression
  #may not be working... produces rediculously small p-values
  #source("http://bioconductor.org/biocLite.R");
  biocLite("limma");
  library("limma");
  
  condition=c(rep(x="control",times=Data$conCount), rep(x="cancer", times=Data$canCount));
  combn <- factor(paste(pData(phenoData)[,1], pData(phenoData)[,2], sep = "_"));
  design <- model.matrix(~condition);# describe model to be fit
  
  fit <- lmFit(Data$ma, design);# fit each probeset to model
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
    maGenes<- as.matrix(cbind(efit.p.adj[filter], Data$ma[filter,]));#attach adjusted p value now
    filter <- res$padj<cutoff & !is.na(res$padj);
    rsGenes<- as.matrix(cbind(res$padj[filter], Data$rs_DESeq[filter,]));#attach adjusted p value now
  }
  else #subset by rank top X most significant
  {
    print("Filtering genes by fold-change:")
    cutoff <- 12000;
    filter <- rank(efit.p.adj)<=cutoff & !is.na(efit.p.adj);
    maGenes<- as.matrix(cbind(efit.p.adj[filter], Data$ma[filter,]));#attach adjusted p value now
    filter <- rank(res$padj)<=cutoff & !is.na(res$padj);
    rsGenes<- as.matrix(cbind(res$padj[filter], Data$rs_DESeq[filter,]));#attach adjusted p value now
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
  maGenes <- Data$ma;
  rsGenes <- Data$rs_DESeq;
  cutoff <- "allGenes";
}

if(args$QCFlag)
{
  print("Outputing quality control figures:");
  
  maRankData<-apply(maGenes,MARGIN=2,FUN=rank);
  rsRankData<-apply(rsGenes,MARGIN=2,FUN=rank);
  
  #plot means
  source("CoexpressionNetworkProject/plot2Groups.R");
  
  plot2Groups(rowMeans(Data$ma, na.rm = TRUE), log(rowMeans(Data$rs_DESeq, na.rm = TRUE)),main="Micro Array vs RNASeq DESeq gene means (97 paired patient samples)",xlab="Lowess Normalized MicroArray",ylab="Log RNASeq DEseq counts", file="Comp_gene_means_across_tech_DESeq.png", histA=TRUE, histB=TRUE);
  plot2Groups(apply(Data$ma,1,median, na.rm = TRUE), apply(log(Data$rs_DESeq),1,median, na.rm = TRUE),main="Micro Array vs RNASeq DESeq gene medians (91 paired patient samples)",xlab="Lowess Normalized MicroArray",ylab="Log RNASeq DESeq counts", file="Comp_gene_medians_across_tech_DESeq.png", histA=TRUE, histB=TRUE);
  plot2Groups(rowMeans(apply(Data$ma,MARGIN=2,FUN=rank), na.rm = TRUE), rowMeans(apply(Data$rs_DESeq,MARGIN=2,FUN=rank), na.rm = TRUE),main="Micro Array vs RNASeq DESeq gene mean ranks (97 paired patient samples)",xlab="Ranked Lowess Normalized MicroArray",ylab="Log RNASeq DESeq counts", file="Ranked Comp_gene_mean_rank_across_tech_DESeq.png");
  plot2Groups(apply(apply(Data$ma,MARGIN=2,FUN=rank), 1,median, na.rm = TRUE), apply(apply(Data$rs_DESeq,MARGIN=2,FUN=rank),1,median, na.rm = TRUE),main="Micro Array vs RNASeq DESeq gene median ranks (97 paired patient samples)",xlab="Ranked Lowess Normalized MicroArray",ylab="Ranked Log RNASeq DESeq counts", file="Comp_gene_median_rank_across_tech_DESeq.png");
  
  plot2Groups(rowMeans(Data$ma, na.rm = TRUE), log(rowMeans(rsData_quant, na.rm = TRUE)),main="Micro Array vs RNASeq quant gene means (97 paired patient samples)",xlab="Lowess Normalized MicroArray",ylab="Log RNASeq quant counts", file="Comp_gene_means_across_tech_quant.png", histA=TRUE, histB=TRUE);
  
  #check normalizations with boxplots
  png(filename="MicroArray normalization check.png");
  boxplot(x=Data$ma,names=seq(1,dim(Data$ma)[2]), outcex=0.5, outpch=20, main="Patient box plots", xlab="Patient", ylab="expression value");
  dev.off();
  png(filename="RNASeq Count normalization check.png");
  boxplot(x=log(Data$rs_raw),names=seq(1,Data$conCount+Data$canCount), outcex=0.5, outpch=20, main="Patient box plots", xlab="Patient", ylab="log RNA seq counts");
  dev.off();
  png(filename="RNASeq RPKM normalization check.png");
  boxplot(x=log(Data$rs_RPKM),names=seq(1,Data$conCount+Data$canCount), outcex=0.5, outpch=20, main="Patient box plots", xlab="Patient", ylab="log RNA seq RPKM");
  dev.off();
  png(filename="RNASeq DESeq normalization check.png");
  boxplot(x=log(Data$rs_DESeq),names=seq(1,Data$conCount+Data$canCount), outcex=0.5, outpch=20, main="Patient box plots", xlab="Patient", ylab="log RNA seq DESeq normalized");
  dev.off();
  png(filename="RNASeq quantile normalization check.png");
  boxplot(x=log(Data$rs_quant),names=seq(1,Data$conCount+Data$canCount), outcex=0.5, outpch=20, main="Patient box plots", xlab="Patient", ylab="log RNA seq quan normalized");
  dev.off();
  png(filename="RNASeq total ubiquitous normalization check.png");
  boxplot(x=log(Data$rs_Ubi),names=seq(1,Data$conCount+Data$canCount), outcex=0.5, outpch=20, main="Patient box plots", xlab="Patient", ylab="log RNA seq total ubiquitous normalized");
  dev.off();
}

#calculated correlation statistics

correlationHistogram <- function(data, method, breaks=100, file)
{
  corrMat <- cor(x=data, method=method, use="complete.obs");
  hist <- hist(x=corrMat,breaks=breaks,plot=FALSE);
  write.csv(x=corrMat,file=file);
  return(list(corrMat=corrMat, hist=hist));
}

library("igraph");
library("ggplot2"); 

for(method in c("pearson","spearman"))
{
  print(paste("Computing ", method, " correlation for:"));

  print("Constructing correlation matricies");
  print("Calculating Microarray correlation matrix:");
  
  density_data <- data.frame();
  
  for(i in 1:length(Data))
  {
    normalization = Data[[i]];
    name <- names(Data)[i];
    profile <- correlationHistogram(data=Data[[i]], method=method, file=paste("Data/BRCA/", name,"_", method, "_",cutoff,"_int.txt"));
    density_data <- data.frame(cor=c(density_data$cor, profile$hist$mids),
    density=c(density_data$density, profile$hist$counts/sum(profile$hist$counts)),
    method=c(density_data$method, rep(x = name, times=length(profile$hist$counts))),stringsAsFactors=FALSE);
    profile$corrMat<-NULL;
  }
  
  #plot overlapping histogram of PCC
  print(paste("Outputting comparative ", method, " histogram:"));

  # Density plots

 # png(filename=);
  ggplot(data=density_data, aes(x=cor, y=density, group=method, colour=method)) + 
    geom_line(size=1, aes(linetype=method)) +
    ggtitle(paste(method, " density comparison"));
  ggsave(paste("Comparative density of network ", method, ".png"));
  #dev.off();
  
  rm(density_data);


  #coexpression networks direct comparison
  maCorrMat <- cor(x=Data$ma, method=method, use="complete.obs");
  rsCorrMat <- cor(x=Data$rs_DESeq, method=method, use="complete.obs");

  if(args$diffCoexFlag)
  {
    print("Calculating differential coexpression network.");
    diffCorrMat <- rsCorrMat - maCorrMat
    write.csv(x=diffCorrMat,file=paste("Data/BRCA/Differential Network ",method,".txt"));
    #calc histogram
    hist <- hist(x=diffCorrMat,breaks=100,plot=FALSE);
    density_data = data.frame("cor"=hist$mids,"density"=hist$counts/sum(hist$counts),"method"=rep(x="diff",times=length(hist$counts)));
    
    #output density distribution
    #png(filename=paste("Differential Network ", method, " density.png"));
    ggplot(data=density_data, aes(x=cor, y=density, group=method, colour=method)) + 
      geom_line(size=1, aes(linetype=method)) +
      ggtitle(paste("Density of ", method,"(rs) - ", method,"(ma)"));
    ggsave(paste("Differential Network ", method, " density.png"));
    #dev.off();
    
    rm(density_data, hist);
    #output venn diagram gene-edge lists
  }

  #create iGraph and plot

  print("Constructing microArray iGraph.");
  maGraph <-graph.adjacency(adjmatrix=maCorrMat*1000,mode="undirected", weighted=TRUE);
  rm(maCorrMat);
  #add vertex attributes to graph
  for(i in 1:length(V(maGraph)))
  {
    name <- V(maGraph)[i]$name;
    #V(maGraph)[i]$fc <- maFC[name];
    #V(maGraph)[i]$p <- efit.p.adj[name];
  }

  print("Outputting microArray iGraph.");
  write.graph(maGraph, file=paste(method, "_maGraph_", cutoff, ".graphml"), format="graphml" );
  rm(maGraph);

  print("Constructing RNASeq iGraph.");
  rsGraph <-graph.adjacency(adjmatrix=rsCorrMat*1000,mode="undirected", weighted=TRUE);
  rm(rsCorrMat);
  for(i in 1:length(V(rsGraph)))
  {
    name <- V(rsGraph)[i]$name;
    #V(rsGraph)[i]$fc <- rsFC[name];
    #V(rsGraph)[i]$p <- res[name,]$padj;
  }

  print("Outputting RNASeq iGraph");
  write.graph(rsGraph, file=paste(method, "_rsGraph_", cutoff, ".graphml"), format="graphml" );
  rm(rsGraph);

  if(args$diffCoexFlag)
  {
   
    print("Constructing differential coexpression iGraph");
    diffGraph <-graph.adjacency(adjmatrix=diffCorrMat*1000,mode="undirected", weighted=TRUE);
    rm(diffCorrMat);
    for(i in 1:length(V(diffGraph)))
    {
      name <- V(diffGraph)[i]$name;
      #V(diffGraph)[i]$RS_fc <- rsFC[name];
      #V(diffGraph)[i]$RS_p <- res[name,]$padj;
      #V(diffGraph)[i]$MA_fc <- maFC[name];
      #V(diffGraph)[i]$MA_p <- efit.p.adj[name];
    }
    print("Outputting differential coexpression iGraph");
    write.graph(diffGraph, file=paste(method, "_diffGraph_", cutoff, ".graphml"), format="graphml" );
    rm(diffGraph);
  }
}
#calculate difference network and fold-change network
#difNet <- rsCorrMat - maCorrMat;
#fcNet <- rsCorrMat / maCorrMat;
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
