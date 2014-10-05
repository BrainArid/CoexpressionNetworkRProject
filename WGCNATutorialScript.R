# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = "C:/Users/Student/My Research/microArray v RNA Seq/Data/WGCNA Tutorial Datasets/";
setwd(workingDir);
# Load the package
source("http://bioconductor.org/biocLite.R")
biocLite("impute")
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
femData = read.csv("FemaleLiver-Data/LiverFemale3600.csv");
# Read in the male liver data set
maleData = read.csv("MaleLiver-Data/LiverMale3600.csv");
# Take a quick look at what is in the data sets (caution, longish output):
dim(femData)
names(femData)
dim(maleData)
names(maleData)

# We work with two sets:
nSets = 2;
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Female liver", "Male liver")
shortLabels = c("Female", "Male")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = as.data.frame(t(femData[-c(1:8)])));
names(multiExpr[[1]]$data) = femData$substanceBXH;
rownames(multiExpr[[1]]$data) = names(femData)[-c(1:8)];
multiExpr[[2]] = list(data = as.data.frame(t(maleData[-c(1:8)])));
names(multiExpr[[2]]$data) = maleData$substanceBXH;
rownames(multiExpr[[2]]$data) = names(maleData)[-c(1:8)];
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)

# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Print information about the removed genes:
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                                              collapse = ", ")))
  for (set in 1:exprSize$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    # Remove the offending genes and samples
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  # Update exprSize
  exprSize = checkSets(multiExpr)
}

sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}

pdf(file = "Plots/SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
dev.off();

# Choose the "base" cut height for the female data set
baseHeight = 16
# Adjust the cut height for the male data set for the number of samples
cutHeights = c(16, 16*exprSize$nSamples[2]/exprSize$nSamples[1]);
# Re-plot the dendrograms including the cut lines
pdf(file = "Plots/SampleClustering.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
{
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);
  abline(h=cutHeights[set], col = "red");
}
dev.off();

#outlier removal
for (set in 1:nSets)
{
  # Find clusters cut by the line
  labels = cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set])
  # Keep the largest one (labeled by the number 1)
  keep = (labels==1)
  multiExpr[[set]]$data = multiExpr[[set]]$data[keep, ]
}
collectGarbage();
3# Check the size of the leftover data
exprSize = checkSets(multiExpr)
exprSize

#Load clinical train data
traitData = read.csv("FemaleLiver-Data/ClinicalTraits.csv");
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData[, -c(31, 16)];
allTraits = allTraits[, c(2, 11:36) ];
# See how big the traits are and what are the trait and sample names
dim(allTraits)
names(allTraits)
allTraits$Mice
# Form a multi-set structure that will hold the clinical traits.
Traits = vector(mode="list", length = nSets);
for (set in 1:nSets)
{
  setSamples = rownames(multiExpr[[set]]$data);
  traitRows = match(setSamples, allTraits$Mice);
  Traits[[set]] = list(data = allTraits[traitRows, -1]);
  rownames(Traits[[set]]$data) = allTraits[traitRows, 1];
}
collectGarbage();
# Define data set dimensions
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples;

save(multiExpr, Traits, nGenes, nSamples, setLabels, shortLabels,
     file = "Consensus-dataInput.RData");
#Tutorial PART 1 complete
