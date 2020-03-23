#This document presents an RNAseq differential expression workflow. We will start from the FASTQ files, align to the reference genome, prepare gene expression values as a count table by counting the sequenced fragments, perform differential gene expression analysis, and visually explore the results.

#We will use publicly available data from the article by Felix Haglund et al., J Clin Endocrin Metab 2012. The purpose of the experiment was to investigate the role of the estrogen receptor in parathyroid tumors. The investigators derived primary cultures of parathyroid adenoma cells from 4 patients. These primary cultures were treated with diarylpropionitrile (DPN), an estrogen receptor beta agonist, or with 4-hydroxytamoxifen (OHT). RNA was extracted at 24 hours and 48 hours from cultures under treatment and control.

#Aligning reads to a reference

#What we get from the sequencing machine is a set of FASTQ files that contain the nucleotide sequence of each read and a quality score at each position. These reads must first be aligned to a reference genome or transcriptome. It is important to know if the sequencing experiment was single-end or paired-end, as the alignment software will require the user to specify both FASTQ files for a paired-end experiment. The output of this alignment step is commonly stored in a file format called BAM.

#Here we use the TopHat2 spliced alignment software in combination with the Bowtie index available at the Illumina iGenomes.

#For example, the paired-end RNA-Seq reads for the parathyroidSE package were aligned using TopHat2 with 8 threads, with the call:

#We will use BAM files from parathyroidSE package to demonstrate how a count table can be constructed from BAM files.

library(parathyroidSE)
#Find extdata

extDataDir <- system.file("extdata", package = "parathyroidSE", mustWork = TRUE)
list.files( extDataDir )


#Typically, we have a table with experimental meta data for our samples. For these three files, it is as follows:

#Sample tables
sampleTable <- data.frame(
  sampleName = c( "Ctrl_24h_1", "Ctrl_48h_1", "DPN_24h_1" ),
  fileName = c( "SRR479052.bam", "SRR479053.bam",  "SRR479054.bam" ),
  treatment =  c( "Control", "Control", "DPN" ),
  time = c( "24h", "48h", "24h" ) )
#This is how the sample table should look like
sampleTable

#Construct the full paths to the files we want to perform the counting operation on:

bamFiles <- file.path( extDataDir, sampleTable$fileName )
bamFiles

#We can peek into one of the BAM files to see the naming style of the sequences (chromosomes). Here we use the BamFile function from the Rsamtools package.

library( "Rsamtools" )
seqinfo( BamFile( bamFiles[1] ) )

#Counting reads in genes

#To count how many read map to each gene, we need transcript annotation. Download the current GTF file with human gene annotation from Ensembl.

#From this file, the function makeTranscriptDbFromGFF from the GenomicFeatures package constructs a database of all annotated transcripts.

#For this lab you can use the truncated version of this file, called Homo_sapiens.GRCh37.75.subset.gtf.gz.

#Preparing gene models from GFF files

library( "GenomicFeatures" )


hse <- makeTxDbFromGFF("Homo_sapiens.GRCh37.75.subset.gtf", format="gtf")
#Extract exons for each gene
exonsByGene <- exonsBy( hse, by="gene" )
exonsByGene

#Reads counting

#The function summarizeOverlaps from the GenomicAlignments package will do this. 

library( "GenomicAlignments" )
se <- summarizeOverlaps( exonsByGene, BamFileList( bamFiles ), mode="Union",
                         singleEnd=FALSE, ignore.strand=TRUE, fragments=TRUE )

se

#We use the counting mode "Union", which indicates that those reads which overlap any portion of exactly one feature are counted. As this experiment produced paired-end reads, we specify singleEnd =FALSE. As protocol was not strand-specific, we specify ignore.strand = TRUE. fragments = TRUE indicates that we also want to count reads with unmapped pairs. This last argument is only for use with paired-end experiments.

# Details on how to read from the BAM files can be specified using the BamFileList function. For example, to control the memory, we could have specified that batches of 2 000 000 reads should be read at a time: 

BamFileList( bamFiles, yieldSize = 2000000 )

#SummarizedExperiment object : Output of counting

#We investigate the resulting SummarizedExperiment class by looking at the counts in the assay slot, the phenotypic data about the samples in colData slot (in this case an empty DataFrame), and the data about the genes in the rowData slot.

# Count data
head( assay(se) )

#Column sum
colSums( assay(se) )

#Phenotypic data
colData(se)

#Data about genes
rowData(se)

#The colData slot, so far empty, should contain all the meta data. We hence assign our sample table to it:

#Column data
colData(se) <- DataFrame( sampleTable )

#We can extract columns from the colData using the $ operator, and we can omit the colData to avoid extra keystrokes.

colData(se) <- DataFrame( sampleTable )
#We can extract columns from the colData using the $ operator, and we can omit the colData to avoid extra keystrokes.

#Column data
colData(se) <- DataFrame( sampleTable )

#We can extract columns from the colData using the $ operator, and we can omit the colData to avoid extra keystrokes.

colData(se)$treatment
se$treatment

#We can also use the sampleName table to name the columns of our data matrix:

colnames(se) <- sampleTable$sampleName
head( assay(se) )

#This SummarizedExperiment object se is then all we need to start our analysis.

#The DESeqDataSet, column metadata, and the design formula

#The data object class in DESeq2 is the DESeqDataSet, which is built on top of the SummarizedExperiment class. One main differences is that the assay slot is instead accessed using the count accessor, and the values in this matrix must be non-negative integers.

#A second difference is that the DESeqDataSet has an associated “design formula”. The design formula tells which variables in the column metadata table colData specify the experimental design and how these factors should be used in the analysis.

#The simplest design formula for differential expression would be ~ condition, where condition is a column in colData(dds) which specifies which of two (or more groups) the samples belong to.

#For the parathyroid experiment, we will specify ~ patient + treatment, which means that we want to test for the effect of treatment (the last factor), controlling for the effect of patient (the first factor).

#We now use R’s data command to load a prepared SummarizedExperiment that was generated from the publicly available sequencing data files associated with the Haglund et al. paper, described on page 1. The steps we used to produce this object were equivalent to those you worked through in the previous Section, except that we used the complete set of samples and all reads.

#Run data
data( "parathyroidGenesSE" )
se <- parathyroidGenesSE
colnames(se) <- se$run
# Check column data
colData(se)[1:5,1:4]

#Here we see that this object already contains an informative colData slot.

#When you work with your own data, you will have to add the pertinent sample / phenotypic information for the experiment at this stage. We highly recommend keeping this information in a comma-separated value (CSV) or tab-separated value (TSV) file, which can be exported from an Excel spreadsheet, and the assign this to the colData slot, as shown in the previous section.

str( metadata( rowData(se) ) )

#Once we have our fully annotated SummerizedExperiment object, we can construct a DESeqDataSet object from it, which will then form the staring point of the actual DESeq2 package

library( "DESeq2" )
ddsFull <- DESeqDataSet( se, design = ~ patient + treatment )


#Collapsing technical replicates

#There are a number of samples which were sequenced in multiple runs. For example, sample SRS308873 was sequenced twice.

#as.data.frame forces R to show us the full list
head(as.data.frame( colData( ddsFull )[ ,c("sample","patient","treatment","time") ] ), 12)

#A convenience function has been implemented to collapse, which can take an object, either SummarizedExperiment or DESeqDataSet, and a grouping factor, in this case the sample name, and return the object with the counts summed up for each unique sample. Optionally, we can provide a third argument, run, which can be used to paste together the names of the runs which were collapsed to create the new object.

ddsCollapsed <- collapseReplicates( ddsFull,
                                    groupby = ddsFull$sample,
                                    run = ddsFull$run )
head( as.data.frame( colData(ddsCollapsed)[ ,c("sample","runsCollapsed") ] ), 12 )

original <- rowSums( counts(ddsFull)[ , ddsFull$sample == "SRS308873" ] )
all( original == counts(ddsCollapsed)[ ,"SRS308873" ] )

#Running the DESeq2 pipeline

#Here we will analyze a subset of the samples, namely those taken after 48 hours, with either control, DPN or OHT treatment, taking into account the multifactor design.

#Preparing the data object for the analysis of interest

#First we subset the relevant columns from the full dataset:

dds <- ddsCollapsed[ , ddsCollapsed$time == "48h" ]

#Sometimes it is necessary to drop levels of the factors, in case that all the samples for one or more levels of a factor in the design have been removed. If time were included in the design formula, the following code could be used to take care of dropped levels in this column.


dds$time <- droplevels( dds$time )

#It will be convenient to make sure that Control is the first level in the treatment factor, so that the default log2 fold changes are calculated as treatment over control and not the other way around. The function relevel achieves this:

dds$treatment <- relevel( dds$treatment, "Control" )

#A quick check whether we now have the right samples:

as.data.frame( colData(dds) )

#In order to speed up some annotation steps below, it makes sense to remove genes which have zero counts for all samples. This can be done by simply indexing the dds object:

dds <- dds[ rowSums( counts(dds) ) > 0 , ]


#Running the pipeline

#Let’s recall what design we have specified:

#design
design(dds)

#Run DESeq : Modeling counts with patient and treatment effects
dds <- DESeq(dds)

#This function will print out a message for the various steps it performs. Briefly these are: the estimation of size factors (which control for differences in the library size of the sequencing experiments), the estimation of dispersion for each gene, and fitting a generalized linear model.

#A DESeqDataSet is returned which contains all the fitted information within it, and the following section describes how to extract out results tables of interest from this object.


#Inspecting the results table

#Calling results without any arguments will extract the estimated log2 fold changes and p values for the last variable in the design formula. If there are more than 2 levels for this variable – as is the case in this analysis – results will extract the results table for a comparison of the last level over the first level. The following section describes how to extract other comparisons.


res <- results( dds )
res

#As res is a DataFrame object, it carries metadata with information on the meaning of the columns:

mcols(res, use.names=TRUE)


#The first column, baseMean, is a just the average of the normalized count values, dividing by size factors, taken over all samples. The remaining four columns refer to a specific contrast, namely the comparison of the levels DPN versus Control of the factor variable treatment.

#See the help page for results (by typing ?results) for information on how to obtain other contrasts.

#The column log2FoldChange is the effect size estimate. It tells us how much the gene’s expression seems to have changed due to treatment with DPN in comparison to control. This value is reported on a logarithmic scale to base 2: for example, a log2 fold change of 1.5 means that the gene’s expression is increased by a multiplicative factor of 21.5 ≈ 2.82.

#Of course, this estimate has an uncertainty associated with it, which is available in the column lfcSE, the standard error estimate for the log2 fold change estimate. The column p value indicates wether the observed difference between treatment and control is significantly different.

#We note that a subset of the p values in res are NA (“notavailable”). This is DESeq’s way of reporting that all counts for this gene were zero, and hence not test was applied. In addition, p values can be assigned NA if the gene was excluded from analysis because it contained an extreme count outlier. For more information, see the outlier detection section of the advanced vignette.

#We can examine the counts and normalized counts for the gene with the smallest p value:

# SmallestPvalue
idx <- which.min(res$pvalue)
counts(dds)[idx, ]


#The results for a comparison of any two levels of a variable can be extracted using the contrast argument to results. The user should specify three values: The name of the variable, the name of the level in the numerator, and the name of the level in the denominator.

#Here we extract results for the log2 of the fold change of DPN/Control:

#Save result table
resOHT <- res
#Other comparisons
res <- results( dds, contrast = c("treatment", "DPN", "Control") )
res


#Adding gene names

#Our result table only uses Ensembl gene IDs, but gene names may be more informative. Bioconductor’s annotation packages help with mapping various ID schemes to each other.

#We load the annotation package org.Hs.eg.db:

library( "org.Hs.eg.db" )

#This is the organism annotation package (“org”“) for Homo sapiens (”Hs“”), organized as an AnnotationDbi package (“db”“), using Entrez Gene IDs (”eg“) as primary key.

#To get a list of all available key types, use

columns(org.Hs.eg.db)

#Converting IDs with the native functions from the AnnotationDbi package is currently a bit cumbersome, so we provide the following convenience function (without explaining how exactly it works):

#ids = list of IDS
#fromKey = key type; toKey = key type we want to convert to
#db = the AnnotationDb object to use.
#ifMultiple = the argument specifies what to do if one source ID maps to several target IDs:
#should the function return an NA or simply the first of the multiple IDs?
convertIDs <- function( ids, fromKey, toKey, db, ifMultiple=c( "putNA", "useFirst" ) ) {
  stopifnot( inherits( db, "AnnotationDb" ) )
  ifMultiple <- match.arg( ifMultiple )
  suppressWarnings( selRes <- AnnotationDbi::select( 
    db, keys=ids, keytype=fromKey, columns=c(fromKey,toKey) ) )
  if( ifMultiple == "putNA" ) {
    duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]   
    selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ] }
  return( selRes[ match( ids, selRes[,1] ), 2 ] )
}


res$hgnc_symbol <- convertIDs( row.names(res), "ENSEMBL", "SYMBOL", org.Hs.eg.db )
res$entrezid <- convertIDs( row.names(res), "ENSEMBL", "ENTREZID", org.Hs.eg.db )
head(res, 4)


#Further points
#Multiple testing

#DESeq2 uses the so-called Benjamini-Hochberg (BH) adjustment for multiple testing problem; in brief, this method calculates for each gene an adjusted p value which answers the following question: if one called significant all genes with a p value less than or equal to this gene’s p value threshold, what would be the fraction of false positives (the false discovery rate, FDR) among them (in the sense of the calculation outlined above)? These values, called the BH-adjusted p values, are given in the column padj of the results object.

#Hence, if we consider a fraction of 10% false positives acceptable, we can consider all genes with an adjusted p value below 10%=0.1 as significant. How many such genes are there?

sum( res$padj < 0.1, na.rm=TRUE )

#We subset the results table to these genes and then sort it by the log2 fold change estimate to get the significant genes with the strongest down-regulation:

resSig <- subset(res, res$padj < 0.1 )
head( resSig[ order( resSig$log2FoldChange ), ], 4)


#and with the strongest upregulation :

head( resSig[ order( -resSig$log2FoldChange ), ], 4)
#Diagnostic plot

#MA plot

#A so-called MA plot provides a useful overview for an experiment with a two-group comparison:

plotMA( res, ylim = c(-3, 3) )

#The MA-plot represents each gene with a dot. The x axis is the average expression over all samples, the y axis the log2 fold change of normalized counts (i.e the average of counts normalized by size factor) between treatment and control. Genes with an adjusted p value below a threshold (here 0.1, the default) are shown in red.

#This plot demonstrates that only genes with a large average normalized count contain sufficient information to yield a significant call.

#Dispersion plot

#Also note DESeq2 shrinkage estimation of log fold changes (LFCs): When count values are too low to allow an accurate estimate of the LFC, the value is “shrunken”" towards zero to avoid that these values, which otherwise would frequently be unrealistically large, dominate the top-ranked log fold change. Whether a gene is called significant depends not only on its LFC but also on its within-group variability, which DESeq2 quantifies as the dispersion. For strongly expressed genes, the dispersion can be understood as a squared coefficient of variation: a dispersion value of 0.01 means that the gene’s expression tends to differ by typically $\sqrt{0.01}=10\%$ between samples of the same treatment group. For weak genes, the Poisson noise is an additional source of noise, which is added to the dispersion. The function plotDispEsts visualizes DESeq2’s dispersion estimates:

plotDispEsts( dds, ylim = c(1e-6, 1e1) )

#The black points are the dispersion estimates for each gene as obtained by considering the information from each gene separately. Unless one has many samples, these values fluctuate strongly around their true values. Therefore, we fit the red trend line, which shows the dispersions’ dependence on the mean, and then shrink each gene’s estimate towards the red line to obtain the final estimates (blue points) that are then used in the hypothesis test. The blue circles above the main “cloud”" of points are genes which have high gene-wise dispersion estimates which are labelled as dispersion outliers. These estimates are therefore not shrunk toward the fitted trend line.

#Histogram of the p-value

hist( res$pvalue, breaks=20, col="grey" )
#Independent filtering

#The MA plot highlights an important property of RNA-Seq data. For weakly expressed genes, we have no chance of seeing differential expression, because the low read counts suffer from so high Poisson noise that any biological effect is drowned in the uncertainties from the read counting. We can also show this by examining the ratio of small p values (say, less than, 0.01) for genes binned by mean normalized count:

# create bins using the quantile function
qs <- c( 0, quantile( res$baseMean[res$baseMean > 0], 0:7/7 ) )
# "cut" the genes into the bins
bins <- cut( res$baseMean, qs )
# rename the levels of the bins using the middle point
levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
# calculate the ratio of p values less than .01 for each bin
ratios <- tapply( res$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) )
# plot these ratios
barplot(ratios, xlab="mean normalized count", ylab="ratio of small p values")


# At first sight, there may seem to be little benefit in filtering out these genes. After all, the test found them to be non-significant anyway. However, these genes have an influence on the multiple testing adjustment, whose performance improves if such genes are removed. By removing the weakly-expressed genes from the input to the FDR procedure, we can find more genes to be significant among those which we keep, and so improved the power of our test. This approach is known as independent filtering. 


#The DESeq software automatically performs independent filtering which maximizes the number of genes which will have adjusted p value less than a critical value (by default, alpha is set to 0.1). This automatic independent filtering is performed by, and can be controlled by, the results function. We can observe how the number of rejections changes for various cutoffs based on mean normalized count. The following optimal threshold and table of possible values is stored as an attribute of the results object.

attr(res,"filterThreshold")

plot(attr(res,"filterNumRej"),type="b",
     xlab="quantiles of 'baseMean'",
     ylab="number of rejections")


# The term independent highlights an important caveat. Such filtering is permissible only if the filter criterion is independent of the actual test statistic. Otherwise, the filtering would invalidate the test and consequently the assumptions of the BH procedure. This is why we filtered on the average over all samples: this filter is blind to the assignment of samples to the treatment and control group and hence independent. 


#Exporting results

#You can easily save the results table in a CSV file, which you can then load with a spreadsheet program such as Excel:

res[1:2,]
write.csv( as.data.frame(res), file="results.csv" )



#Gene-set enrichment analysis

#Do the genes with a strong up- or down-regulation have something in common? We perform next a gene-set enrichment analysis (GSEA) to examine this question.

#We here present a relatively simplistic approach, to demonstrate the basic ideas, but note that a more careful treatment will be needed for more definitive results.

#We use the gene sets in the Reactome database:


#source("http://bioconductor.org/biocLite.R")
#biocLite("reactome.db")
library( "reactome.db" )


#This database works with Entrez IDs, so we will need the entrezid column that we added earlier to the res object.

#First, we subset the results table, res, to only those genes for which the Reactome database has data (i.e, whose Entrez ID we find in the respective key column of reactome.db and for which the DESeq2 test gave an adjusted p value that was not NA.

res2 <- res[ res$entrezid %in% keys( reactome.db, "ENTREZID" ) & !is.na( res$padj ) , ]
head(res2)

#Using select, a function from AnnotationDbi for querying database objects, we get a table with the mapping from Entrez IDs to Reactome Path IDs :

reactomeTable <- AnnotationDbi::select( reactome.db, 
                                        keys=as.character(res2$entrezid), keytype="ENTREZID", 
                                        columns=c("ENTREZID","REACTOMEID") )
head(reactomeTable)

#The next code chunk transforms this table into an incidence matrix. This is a Boolean matrix with one row for each Reactome Path and one column for each unique gene in res2, which tells us which genes are members of which Reactome Paths.

incm <- do.call( rbind, with(reactomeTable, tapply( 
  ENTREZID, factor(REACTOMEID), function(x) res2$entrezid %in% x ) ))
colnames(incm) <- res2$entrez
str(incm)

#We remove all rows corresponding to Reactome Paths with less than 20 or more than 80 assigned genes.

within <- function(x, lower, upper) (x>=lower & x<=upper)
incm <- incm[ within(rowSums(incm), lower=20, upper=80), ]

#To test whether the genes in a Reactome Path behave in a special way in our experiment, we calculate a number of statistics, including a t-statistic to see whether the average of the genes’ log2 fold change values in the gene set is different from zero. To facilitate the computations, we define a little helper function:


testCategory <- function( reactomeID ) {
  isMember <- incm[ reactomeID, ]
  data.frame( 
    reactomeID  = reactomeID,
    numGenes    = sum( isMember ),
    avgLFC      = mean( res2$log2FoldChange[isMember] ),
    sdLFC       = sd( res2$log2FoldChange[isMember] ),
    zValue      = mean( res2$log2FoldChange[isMember] ) /sd( res2$log2FoldChange[isMember] ),
    strength    = sum( res2$log2FoldChange[isMember] ) / sqrt(sum(isMember)),
    pvalue      = t.test( res2$log2FoldChange[ isMember ] )$p.value,
    reactomeName = reactomePATHID2NAME[[reactomeID]],
    stringsAsFactors = FALSE ) }

#The function can be called with a Reactome Path ID:

testCategory("109606")
#Working with rlog-transformed data


#The rlog transform

#Many common statistical methods for exploratory analysis of multidimensional data, especially methods for clustering (e.g., principal-component analysis and the like), work best for (at least approximately) homoskedastic data; this means that the variance of an observable quantity (i.e., here, the expression strength of a gene) does not depend on the mean. In RNA-Seq data, however, variance grows with the mean. For example, if one performs PCA directly on a matrix of normalized read counts, the result typically depends only on the few most strongly expressed genes because they show the largest absolute differences between samples. A simple and often used strategy to avoid this is to take the logarithm of the normalized count values plus a small pseudocount; however, now the genes with low counts tend to dominate the results because, due to the strong Poisson noise inherent to small count values, they show the strongest relative differences between samples.

#As a solution, DESeq2 offers the regularized-logarithm transformation, or rlog for short. For genes with high counts, the rlog transformation differs not much from an ordinary log2 transformation. For genes with lower counts, however, the values are shrunken towards the genes’ averages across all samples. Using an empirical Bayesian prior in the form of a ridge penalty, this is done such that the rlog-transformed data are approximately homoskedastic.

#Note that the rlog transformation is provided for applications other than differential testing. For differential testing we recommend the DESeq function applied to raw counts, as described earlier in this vignette, which also takes into account the dependence of the variance of counts on the mean value during the dispersion estimation step.

#The function rlog returns a SummarizedExperiment object which contains the rlog-transformed values in its assay slot:

rld <- rlog( dds )
assay(rld)[ 1:3, 1:3]


#To show the effect of the transformation, we plot the first sample against the second, first simply using the log2 function (after adding 1, to avoid taking the log of zero), and then using the rlog-transformed values.


plot( log2( 1+counts(dds, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3 )
plot( assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3 )


#Sample distances

#A useful first step in an RNA-Seq analysis is often to assess overall similarity between samples.

#We use the R function dist to calculate the Euclidean distance between samples. To avoid that the distance measure is dominated by a few highly variable genes, and have a roughly equal contribution from all genes, we use it on the rlog-transformed data:

sampleDists <- dist( t( assay(rld) ) )
as.matrix( sampleDists )[ 1:3, 1:3 ]


# Note the use of the function t to transpose the data matrix. We need this because dist calculates distances between data rows and our samples constitute the columns. 

#We visualize the distances in a heatmap, using the function heatmap.2 from the gplots package.

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$treatment, 
                                     rld$patient, sep="-" )
colnames(sampleDistMatrix) <- NULL   
library( "gplots" )
library( "RColorBrewer" )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)


#Another way to visualize sample-to-sample distances is a principal-components analysis (PCA). In this ordination method, the data points (i.e., here, the samples) are projected onto the 2D plane such that they spread out optimally.


colours <- c(rgb(1:3/4,0,0),rgb(0,1:3/4,0),rgb(0,0,1:3/4),rgb(1:3/4,0,1:3/4))
plotPCA( rld, intgroup = c("patient","treatment"), col=colours )

#From both visualizations, we see that the differences between patients is much larger than the difference between treatment and control samples of the same patient. This shows why it was important to account for this paired design (``paired’’, because each treated sample is paired with one control sample from the same patient).

#We did so by using the design formula ~ patient + treatment when setting up the data object in the beginning. Had we used an un-paired analysis, by specifying only , we would not have found many hits, because then, the patient-to-patient differences would have drowned out any treatment effects.

#Gene clustering

#In the above heatmap, the dendrogram at the side shows us a hierarchical clustering of the samples. Such a clustering can also be performed for the genes.

#Since the clustering is only relevant for genes that actually carry signal, one usually carries it out only for a subset of most highly variable genes. Here, for demonstration, let us select the 35 genes with the highest variance across samples:

library( "genefilter" )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )


#The heatmap becomes more interesting if we do not look at absolute expression strength but rather at the amount by which each gene deviates in a specific sample from the gene’s average across all samples. Hence, we center and scale each genes’ values across samples, and plot a heatmap.

heatmap.2( assay(rld)[ topVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
           ColSideColors = c( Control="gray", DPN="darkgreen", OHT="orange" )[
             colData(rld)$treatment ] )

#We can now see on the heatmap, blocks of genes which covary across patients. Often, such a heatmap is insightful, even though here, seeing these variations across patients is of limited value because we are rather interested in the effects between the treatments from each patient. 

#Analyze more datasets: use the function defined in the following code chunk to download a processed count matrix from the ReCount website.

#The following function takes a name of the dataset from the ReCount website, e.g. “hammer”, and returns a SummarizedExperiment object.

recount2SE <- function(name) {
  filename <- paste0(name,"_eset.RData")
  if (!file.exists(filename)) download.file(paste0(
    "http://bowtie-bio.sourceforge.net/recount/ExpressionSets/",
    filename),filename)
  load(filename)
  e <- get(paste0(name,".eset"))
  se <- SummarizedExperiment(SimpleList(counts=exprs(e)),
                             colData=DataFrame(pData(e)))
  se                   
}

#For a more in-depth explanation of the advanced details, we advise you to proceed to the vignette of the DESeq2 package package, Differential analysis of count data. For a treatment of exon-level differential expression, we refer to the vignette of the DEXSeq package, Analyzing RN-seq data for differential exon usage with the DEXSeq package.

#DEXSeq for differential exon usage. See the accompanying vignette, Analyzing RNA-seq data for differential exon usage with the DEXSeq package, which is similar to the style of this tutorial.

#Other Bioconductor packages for RNA-Seq differential expression:

#edgeR, limma, DSS, BitSeq (transcript level), EBSeq, cummeRbund (for importing and visualizing Cufflinks results), monocle (single-cell analysis). More at http://bioconductor.org/packages/release/BiocViews.html#___RNASeq

#Methods for gene set testing: romer and roast in limma, permutation based: safe 
#Packages for normalizing for covariates (e.g., GC content): cqn, EDASeq 
#Packages for detecting batches: sva, RUVSeq 
#Generating HTML results tables with links to outside resources (gene descriptions): ReportingTools

#Session Info

#As last part of this document, we call the function , which reports the version numbers of R and all the packages used in this session. It is good practice to always keep such a record as it will help to trace down what has happened in case that an R script ceases to work because a package has been changed in a newer version.

  