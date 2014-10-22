setwd("/Users/Brian/Documents/Research/microArray v RNA Seq/BRCA/")
maClusts <- read.table(file="Data//BRCA//maPearson_top12000_int.clusters.csv",sep = ",",stringsAsFactors = FALSE)
rsClusts <- read.table(file="Data//BRCA//rsPearson_top12000_int.clusters.csv",sep = ",",stringsAsFactors = FALSE)
maClusts <- maClusts[,-1];
rsClusts <- rsClusts[,-1];
countsMat <- matrix(nrow = dim(maClusts)[1], ncol= dim(rsClusts)[1]);
countsMat[,]<-0;

for(i in 1:dim(maClusts)[1])
{
  for(j in 1:dim(rsClusts)[1])
  {
    print(paste0("Cross comparing microarray row ", i, " of ", dim(maClusts)[1]," and column ", j, " of ", dim(rsClusts)[1]));
    countsMat[i,j] <- length(intersect(maClusts[[i]],rsClusts[[j]]))
  }
}

maLengths <- mapply(maClusts, FUN=length)
rsLengths <- mapply(rsClusts, FUN=length)

library("ggplot2")

