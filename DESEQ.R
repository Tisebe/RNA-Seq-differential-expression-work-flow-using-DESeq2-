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
