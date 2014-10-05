constructCrossSampleFrame_privateHelper <- function(outFrame, inFiles, data, cols2Ignore, i)
{
  
  if(is.null(cols2Ignore))
  {
    colnames(outFrame)[i] <- paste(tail(strsplit(inFiles[i], split="/")[[1]],1), "_", colnames(data)[-1]);
    outFrame[,i] <- as.numeric(data[,-1]);
  }
  else
  {
    colnames(outFrame)[i] <- paste(tail(strsplit(inFiles[i], split="/")[[1]],1), "_", colnames(data)[-cols2Ignore][-1]);
    outFrame[,i] <- as.numeric(data[,-cols2Ignore][,-1]);
  }
  
  return(outFrame);
}

constructCrossSampleFrame <- function(inFiles, rows2Ignore=c(), cols2Ignore=c())
{
  data <- read.csv(file=inFiles[1], sep="\t", stringsAsFactors=FALSE);
  if(!is.null(rows2Ignore))
  {
    data <- data[-rows2Ignore,];
  }
  outFrame <- data.frame(matrix(ncol = length(inFiles)[1]*(dim(data)[2]-length(cols2Ignore)-1), nrow = dim(data)[1]));
  row.names(outFrame) <- data[,1];
  outFrame <- constructCrossSampleFrame_privateHelper(outFrame, inFiles, data, cols2Ignore, 1);
  i<-2
  for (file in inFiles[2:length(inFiles)])
  {
    data <- read.csv(file=file, sep="\t", stringsAsFactors=FALSE);
    if(!is.null(rows2Ignore))
    {
      data <- data[-rows2Ignore,];
    }
    outFrame <- constructCrossSampleFrame_privateHelper(outFrame, inFiles, data, cols2Ignore, i);
    i<-i+1
  }
    
  return(outFrame);
}