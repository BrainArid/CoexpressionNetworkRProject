trim_TCGA_RNASeq_GeneNames <- function(x)
{
  x <- x[(sapply(X=row.names(x), FUN=substr, start=1, stop=1)!="?"),];#remove Genes whose names are "?"
  potentialRowNames <- vector(length=dim(x)[1]);
  for(i in 1:length(row.names(x)))
  {
    potentialRowNames[i] <- head(strsplit(x=row.names(x)[i], split="|", fixed=TRUE)[[1]], 1);
  }
  #remove duplicates
  x <- x[!(duplicated(potentialRowNames)|duplicated(potentialRowNames, fromLast=TRUE)),];
  row.names(x) <- potentialRowNames[!(duplicated(potentialRowNames)|duplicated(potentialRowNames, fromLast=TRUE))];
  return(x);
}