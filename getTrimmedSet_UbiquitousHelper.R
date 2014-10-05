getTrimmedSet_UbiquitousHelper <-function(sample, lowerPercentile, upperPercentile)
{ 
  thresholds <- quantile(sample, c(lowerPercentile, upperPercentile));
  return(sample[thresholds[1] <= sample & sample <= thresholds[2]]);
}