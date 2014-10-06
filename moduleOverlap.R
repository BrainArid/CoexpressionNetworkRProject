setwd("/Users/Brian/Documents/Research/microArray v RNA Seq/BRCA/")
maClusts <- read.table(file="maPearson_top12000_int.clusters.csv",sep = ",",stringsAsFactors = FALSE)
rsClusts <- read.table(file="rsPearson_top12000_int.clusters.csv",sep = ",",stringsAsFactors = FALSE)
maClusts <- maClusts[,-1];
rsClusts <- rsClusts[,-1];
countsMat <- matrix(nrow = dim(maClusts)[1], ncol= dim(rsClusts)[1]);
countsMat[,]<-0;

for(i in 1:dim(maClusts)[1])
{
  print(paste("Cross comparing microarray row ", i, " of ", dim(maClusts)[1]));
  row = maClusts[i,maClusts[i,]!=''];
  for(r in 1:length(row))
  {
    value = row[r];
    for(j in 1:dim(rsClusts)[1])
    {
      rsClust <- rsClusts[j,rsClusts[j,]!='']
      if(value %in% rsClust)
      {
        countsMat[i, j] = countsMat[i, j] + 1;
      }
    }
  }
  print(paste("Found ", sum(countsMat[i,])));
}
