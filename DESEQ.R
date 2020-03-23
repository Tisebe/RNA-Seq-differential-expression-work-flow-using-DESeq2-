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

