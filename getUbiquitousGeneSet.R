getUbiquitousGeneSet <-function(samples, lowerPercentile, upperPercentile)
{ 
  source("getTrimmedSet_UbiquitousHelper.R")
  i<-1;
  sample <- samples[,i];
  geneSet <- getTrimmedSet_UbiquitousHelper(sample, lowerPercentile, upperPercentile);
  ubiquitous <- names(geneSet);
  
  #for each sample get trimmed set
  for(i in 2:dim(samples)[2])
  {
    sample <- samples[,i];
    geneSet <- getTrimmedSet_UbiquitousHelper(sample, lowerPercentile, upperPercentile);
    ubiquitous <- intersect(x=ubiquitous,y=names(geneSet));
  }
  
  return(ubiquitous);
}