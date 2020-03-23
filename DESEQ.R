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

