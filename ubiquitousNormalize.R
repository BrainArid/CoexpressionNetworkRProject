ubiquitousNormalize <- function(data, lowerPercentile, upperPercentile)
{
  #source("getUbiquitousGeneSet.R");
  
  ubiquitousGenes <- getUbiquitousGeneSet(data, lowerPercentile, upperPercentile);
  ubiData <- data[ubiquitousGenes,];
  
  #calculate average total counts
  counts <- vector(length=dim(data)[2]);
  for(i in 1:length(counts))
  {
    counts[i] <- sum(ubiData[,i])
  }
  meanCount <- mean(counts);
  for(i in 1:length(counts))
  {
    data[,i] <- data[,i] * (meanCount / counts[i]);
  }
  
  return(data);
}