


RNASeq Practical
========================================================
author: CSC Bioinformatics Core
date: 17 December 2015
width: 1640
height: 1200
autosize: true
font-import: <link href='http://fonts.googleapis.com/css?family=Slabo+27px' rel='stylesheet' type='text/css'>
font-family: 'Slabo 27px', serif;
css:style.css


Contents
====================================
- Building index
- Splice aware alignment
- Quality assessment
- Assigning reads to genomic features
- Differential gene expression analysis
- Exploring and saving results
- Differential splicing analysis
- Gene ontology and pathway enrichment analysis


External libraries
====================================


```r
### load libraries

library(Rsubread)
library(limma)
library(edgeR)
library(DESeq2)
library(goseq)
library("RColorBrewer")
library(ggplot2)
library(KEGG.db)
library("org.Mm.eg.db")
library(gplots)
```

Analysis Parameters
========================================================


```r
#  Define analysis parameters

genome<-"mm9" 
strandspecific<-0 
factor1<-"Group"
dir="/Users/skhadaya/Downloads/sratoolkit.2.5.5-mac64/bin/"
isPairedEnd<-TRUE
nthreads=8
bamdir="/Users/skhadaya/Downloads/rnaseqpractical/"
```


Reading sample data
========================================================

```r
# Read in target file using a funcion from limma package

targets <- readTargets("targets.txt")
```

========================================================

```r
targets
```

```
  Sample Group Batch        InputFile       InputFile2 OutputFile
1   Viv1   Viv     a SRR2001243.fastq SRR2001244.fastq   Viv1.bam
2   Viv2   Viv     a SRR2001245.fastq SRR2001246.fastq   Viv2.bam
3   Viv3   Viv     b SRR2001247.fastq SRR2001248.fastq   Viv3.bam
4   Hfd1   Hfd     a SRR2001249.fastq SRR2001250.fastq   Hfd1.bam
5   Hfd2   Hfd     b SRR2001251.fastq SRR2001252.fastq   Hfd2.bam
6   Hfd3   Hfd     a SRR2001253.fastq SRR2001254.fastq   Hfd3.bam
```

========================================================


```r
# Get full path names of raw data files


read1=paste0(dir,targets$InputFile)

read2=paste0(dir,targets$InputFile2)
```

========================================================


```r
read1[1]
```

```
[1] "/Users/skhadaya/Downloads/sratoolkit.2.5.5-mac64/bin/SRR2001243.fastq"
```

```r
read2[1]
```

```
[1] "/Users/skhadaya/Downloads/sratoolkit.2.5.5-mac64/bin/SRR2001244.fastq"
```

========================================================


```r
# get the directory with rsubread genome index

index="index"
```

Building index
========================================================


```r
We will use rsubread package to perform read alignment and counting.

#  Step 1: Index building

ref <- system.file("extdata","reference.fa",package="Rsubread")

buildindex(basename="reference_index",reference=ref)
```

Read Alignment
========================================================


```r
# align reads using subjunc function of rsubreads package

if (isPairedEnd == TRUE)
{
  subjunc(index=index,
  readfile1=read1,
  readfile2=read2,
  input_format="gzFASTQ",
  output_format="BAM",
  output_file=paste0(bamdir,targets$OutputFile),
  nthreads=nthreads,
  unique=TRUE,
  indels=5)
}
```

========================================================


```r
else {
  subjunc(index=index,
  readfile1=read1,
  input_format="gzFASTQ",
  output_format="BAM",
  output_file=paste0(bamdir,targets$OutputFile),
  useMetaFeatures=TRUE,
  nthreads=nthreads,
  unique=TRUE,
  indels=5)
}
```

Quality Assessment
========================================================


```r
# Quality scores give the probabilities of read bases being incorrectly called, which is useful for examining the quality of sequencing data. This function extract quality strings and convert them to Phred scores.

# Base quality scores
qscore <- qualityScores(filename=read1[1],offset=33,nreads=1000)
```

========================================================


```r
plot(colMeans(qscore),ylim=c(0,40))
```

![plot of chunk unnamed-chunk-12](practical-figure/unnamed-chunk-12-1.png) 

========================================================


```r
# GC content

# The basewise calculation is useful for examining the GC bias towards the base position in the read.

atgc<-atgcContent(filename=read1[1],basewise = TRUE)

plot(atgc[1,],col="red",ylim=c(0,0.4),ylab="content")
points(atgc[2,],col="green")
points(atgc[3,],col="yellow")
points(atgc[4,],col="brown")
legend("topright", legend=c("A", "T","G","C"), col=c("red","green","yellow","brown"), pch=15)
```

![plot of chunk unnamed-chunk-13](practical-figure/unnamed-chunk-13-1.png) 

========================================================



```r
# produce alignment statistics

propmapped(paste0(bamdir,targets$OutputFile[1]))
```

```
The input file is opened as a BAM file.
The fragments in the input file are being counted.
Finished. All records: 49694842; all fragments: 49694842; mapped fragments: 32819270; the mappability is 66.04%
```

```
                                             Samples NumTotal NumMapped
1 /Users/skhadaya/Downloads/rnaseqpractical/Viv1.bam 49694842  32819270
  PropMapped
1   0.660416
```

Counting Reads
========================================================


```r
# count numbers of reads mapped to genes

anno_for_featurecount<-paste0(bamdir,"genes.gtf")

fc <-featureCounts(
files=paste0(bamdir,targets$OutputFile)
,
annot.ext=anno_for_featurecount,
isGTFAnnotationFile=TRUE,
useMetaFeatures=TRUE,
GTF.featureType="exon",
GTF.attrType="gene_id",
nthreads=nthreads,
strandSpecific=strandspecific,
isPairedEnd=isPairedEnd)
```

========================================================


```r
head(fc$counts)
```

```
                   X.Users.skhadaya.Downloads.rnaseqpractical.Viv1.bam
ENSMUSG00000090025                                                   0
ENSMUSG00000064842                                                   0
ENSMUSG00000051951                                                   0
ENSMUSG00000089699                                                   0
ENSMUSG00000088390                                                   0
ENSMUSG00000089420                                                   0
                   X.Users.skhadaya.Downloads.rnaseqpractical.Viv2.bam
ENSMUSG00000090025                                                   0
ENSMUSG00000064842                                                   0
ENSMUSG00000051951                                                   1
ENSMUSG00000089699                                                   0
ENSMUSG00000088390                                                   0
ENSMUSG00000089420                                                   0
                   X.Users.skhadaya.Downloads.rnaseqpractical.Viv3.bam
ENSMUSG00000090025                                                   0
ENSMUSG00000064842                                                   0
ENSMUSG00000051951                                                   1
ENSMUSG00000089699                                                   0
ENSMUSG00000088390                                                   0
ENSMUSG00000089420                                                   0
                   X.Users.skhadaya.Downloads.rnaseqpractical.Hfd1.bam
ENSMUSG00000090025                                                   0
ENSMUSG00000064842                                                   0
ENSMUSG00000051951                                                   3
ENSMUSG00000089699                                                   0
ENSMUSG00000088390                                                   0
ENSMUSG00000089420                                                   0
                   X.Users.skhadaya.Downloads.rnaseqpractical.Hfd2.bam
ENSMUSG00000090025                                                   0
ENSMUSG00000064842                                                   0
ENSMUSG00000051951                                                   0
ENSMUSG00000089699                                                   0
ENSMUSG00000088390                                                   0
ENSMUSG00000089420                                                   0
                   X.Users.skhadaya.Downloads.rnaseqpractical.Hfd3.bam
ENSMUSG00000090025                                                   0
ENSMUSG00000064842                                                   0
ENSMUSG00000051951                                                   0
ENSMUSG00000089699                                                   0
ENSMUSG00000088390                                                   0
ENSMUSG00000089420                                                   0
```


========================================================


```r
head(fc$annotation)
```

```
              GeneID     Chr                           Start
1 ENSMUSG00000090025       1                         3044314
2 ENSMUSG00000064842       1                         3092097
3 ENSMUSG00000051951 1;1;1;1 3195982;3203520;3411783;3660633
4 ENSMUSG00000089699     1;1                 3456668;3503486
5 ENSMUSG00000088390       1                         3668961
6 ENSMUSG00000089420       1                         3773818
                              End  Strand Length
1                         3044814       +    501
2                         3092206       +    110
3 3197398;3207049;3411982;3661579 -;-;-;-   6094
4                 3456768;3503634     +;+    250
5                         3669024       -     64
6                         3773879       -     62
```


========================================================


```r
  # collect sample information
                              
 colData<-cbind(targets$OutputFile,targets$Group,targets$Batch)
 
rownames(colData)<-colData[,1]
colnames(colData)<-c("name","Group","Batch")
                             
colData<-data.frame(colData)
```

========================================================


```r
 colData
```

```
             name Group Batch
Viv1.bam Viv1.bam   Viv     a
Viv2.bam Viv2.bam   Viv     a
Viv3.bam Viv3.bam   Viv     b
Hfd1.bam Hfd1.bam   Hfd     a
Hfd2.bam Hfd2.bam   Hfd     b
Hfd3.bam Hfd3.bam   Hfd     a
```

Prepare deseqdataset object
========================================================


```r
 # The class used by the DESeq2 package to store the read counts is DESeqDataSet
 
# construct deseqdataset object
 
  dds<-DESeqDataSetFromMatrix(
 countData= fc$counts,colData=targets,
 design=~Group)
```

Differential expression analysis
========================================================


```r
# The standard differential expression analysis steps are wrapped into a single function, DESeq
# Perform normalization, fitting to the model
 
    dds<-DESeq(dds)
```



========================================================


```r
# The function DESeq runs the following functions in order:
   
# 1. estimation of size factors: 
#  estimateSizeFactors()

# The sizeFactors vector assigns to each column of the count matrix a value, the size factor, such that  count  values  in  the  columns  can  be  brought  to  a  common  scale  by  dividing  by  the  corresponding size factor.

sizeFactors(dds)
```

```
        1         2         3         4         5         6 
1.2430187 0.7755226 1.0501449 0.9457439 1.0124687 1.0515602 
```

========================================================


```r
# 2. estimation of dispersion: 
#  estimateDispersions()

# This function obtains dispersion estimates for Negative Binomial distributed data.

head(dispersions(dds))
```

```
[1]       NA       NA 2.381315       NA       NA       NA
```

```r
plotDispEsts(dds)
```

![plot of chunk unnamed-chunk-23](practical-figure/unnamed-chunk-23-1.png) 

========================================================


```r
# 3. Negative Binomial GLM fitting and Wald statistics: 
#  nbinomWaldTest()

# This function tests for significance of coefficients in a Negative Binomial GLM, using previously calculated sizeFactors (or normalizationFactors ) and dispersion estimates.
```

Transformation of count data
========================================================


```r
# The function rlog , stands for regularized log , transforming the original count data to the log2 scale
# Aim of this transformation,  the rlog and  the  VST,  is  to  remove  the  dependence  of  the  variance
# on  the  mean,  so that data is suitable for visualization. 
# particularly  the  high  variance  of  the  logarithm  of  count  data  when  the  mean  is  low. 
# 

      rld<-rlog(dds)  
```

Principal component plot of the samples
========================================================


```r
      plotPCA(rld, intgroup="Group")
```

![plot of chunk unnamed-chunk-26](practical-figure/unnamed-chunk-26-1.png) 


Heatmap of sample to sample distances
========================================================


```r
      # The assay function is used to extract the matrix of normalized values
rlogcount <- assay(rld)
      
rlogcount <- rlogcount[!rowSums(rlogcount) == 0,]

colnames(rlogcount) <-  paste0(colData(dds)$sample)                          
sampleDists <- as.matrix(dist(t(rlogcount)))
 showcols <- brewer.pal(8, "Dark2")[1:length(unique(colData(dds)$Group))]
```

========================================================


```r
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
col=colorpanel(100, "black", "white"),
ColSideColors=showcols[colData(dds)$Group], 
RowSideColors=showcols[colData(dds)$Group],
margin=c(10, 10), main="Sample Distance Matrix")
```

========================================================

![plot of chunk unnamed-chunk-29](practical-figure/unnamed-chunk-29-1.png) 


Exporting results
========================================================


```r
 # Results tables are generated using the function results , which extracts a results table with log2 fold changes, p  values and adjusted p values
 
    res<-results(dds) 
 # Order results by adjusted p value 
 
        resOrdered<-res[order(res$padj),]
    
      write.table(resOrdered,file="DEgenes.txt",sep="\t")
```

========================================================


```r
  summary(res)
```

```

out of 22605 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)     : 990, 4.4% 
LFC < 0 (down)   : 1481, 6.6% 
outliers [1]     : 19, 0.084% 
low counts [2]   : 5493, 24% 
(mean count < 3)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

```r
  # How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)
```

```
[1] 2471
```

Exploring results
========================================================


```r
      # MA plot
      
  # The  function plotMA shows  the  log2  fold  changes  attributable  to  a  given  variable  over  the  mean of normalized counts.  Points will be colored red if the adjusted p value is less than 0.1.  Points which fall out of the window are plotted as open triangles pointing either up or down.
      
    plotMA(res, main="DESeq2", ylim=c(-4,4))
```


========================================================

![plot of chunk unnamed-chunk-33](practical-figure/unnamed-chunk-33-1.png) 


========================================================


```r
      # Plot counts
   # Plot of normalized counts for a single gene on log scale

    plotCounts(dds,gene=which.min(res$padj),intgroup="Group")
```

![plot of chunk unnamed-chunk-34](practical-figure/unnamed-chunk-34-1.png) 


Multi factor designs
========================================================


```r
#   Experiments  with  more  than  one  factor  influencing  the  counts  can  be  analyzed  using  design  formula  that  include  the  additional  variables. 

targets
```

```
  Sample Group Batch        InputFile       InputFile2 OutputFile
1   Viv1   Viv     a SRR2001243.fastq SRR2001244.fastq   Viv1.bam
2   Viv2   Viv     a SRR2001245.fastq SRR2001246.fastq   Viv2.bam
3   Viv3   Viv     b SRR2001247.fastq SRR2001248.fastq   Viv3.bam
4   Hfd1   Hfd     a SRR2001249.fastq SRR2001250.fastq   Hfd1.bam
5   Hfd2   Hfd     b SRR2001251.fastq SRR2001252.fastq   Hfd2.bam
6   Hfd3   Hfd     a SRR2001253.fastq SRR2001254.fastq   Hfd3.bam
```



========================================================


```r
ddsMF<-DESeqDataSetFromMatrix(countData= fc$counts,colData= colData,design=~ Batch + Group)
      
ddsMF <- DESeq(ddsMF)

resMF <- results(ddsMF)
```



========================================================


```r
resMForder<-resMF[order(resMF$padj),]
head(resMForder)
```

```
log2 fold change (MAP): Group Viv vs Hfd 
Wald test p-value: Group Viv vs Hfd 
DataFrame with 6 rows and 6 columns
                     baseMean log2FoldChange     lfcSE      stat
                    <numeric>      <numeric> <numeric> <numeric>
ENSMUSG00000026475  3141.2571      -4.683048 0.1975390 -23.70695
ENSMUSG00000024526   465.8790      -5.842030 0.2505497 -23.31685
ENSMUSG00000032080 14291.8178      -4.562778 0.2262454 -20.16738
ENSMUSG00000028051   339.3991      -3.352374 0.1984522 -16.89260
ENSMUSG00000069170   343.9283      -3.660025 0.2184593 -16.75381
ENSMUSG00000038576   842.8800       2.144259 0.1348818  15.89733
                          pvalue          padj
                       <numeric>     <numeric>
ENSMUSG00000026475 3.056829e-124 5.230846e-120
ENSMUSG00000024526 2.991369e-120 2.559415e-116
ENSMUSG00000032080  1.894070e-90  1.080377e-86
ENSMUSG00000028051  5.100105e-64  2.181825e-60
ENSMUSG00000069170  5.311611e-63  1.817846e-59
ENSMUSG00000038576  6.612782e-57  1.885965e-53
```

========================================================


```r
# It  is  also  possible  to  retrieve  the  log2  fold  changes, p values  and  adjusted p values  of  the Batch variable

resMFbatch <- results(ddsMF, contrast=c("Batch","a","b"))
resMFbatchOrder<-resMFbatch[order(resMFbatch$padj),]
head(resMFbatchOrder)
```

```
log2 fold change (MAP): Batch a vs b 
Wald test p-value: Batch a vs b 
DataFrame with 6 rows and 6 columns
                    baseMean log2FoldChange     lfcSE      stat
                   <numeric>      <numeric> <numeric> <numeric>
ENSMUSG00000067219  172.4854      1.7566066 0.1795349  9.784204
ENSMUSG00000020122 5025.7129     -0.8415695 0.1209184 -6.959815
ENSMUSG00000041567 1349.3249     -0.9372828 0.1351759 -6.933804
ENSMUSG00000019577  529.3305      1.0465173 0.1557041  6.721191
ENSMUSG00000049580 1487.8130     -0.8759096 0.1331018 -6.580750
ENSMUSG00000020102 1494.3816      0.7013548 0.1148805  6.105079
                         pvalue         padj
                      <numeric>    <numeric>
ENSMUSG00000067219 1.316264e-22 1.696138e-18
ENSMUSG00000020122 3.407200e-12 1.759684e-08
ENSMUSG00000041567 4.096735e-12 1.759684e-08
ENSMUSG00000019577 1.802449e-11 5.806588e-08
ENSMUSG00000049580 4.680805e-11 1.206337e-07
ENSMUSG00000020102 1.027500e-09 2.206727e-06
```

========================================================


```r
summary(resMFbatch)
```

```

out of 22605 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)     : 43, 0.19% 
LFC < 0 (down)   : 39, 0.17% 
outliers [1]     : 0, 0% 
low counts [2]   : 9719, 43% 
(mean count < 28)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```


Likelyhood Ratio Test
========================================================


```r
DESeq2 offers  two  kinds  of  hypothesis  tests: 
  
  1. The  Wald  test,  where  we  use  the  estimated  standard  error  of a  log2  fold  change  to  test  if  it  is  equal  to  zero.

  2. The  likelihood  ratio  test  (LRT).  
  The  LRT  examines  two models  for  the  counts,  a full model  with  a  certain  number  of  terms  and  a
reduced model,  in  which  some  of the terms of the full model are removed.  
The test determines if the increased likelihood of the data using the extra terms in the full model is more than expected if those extra terms are truly zero.
```

========================================================


```r
The LRT is therefore useful for testing multiple terms at once, for example testing 3 or more levels of a factor
at once, or all interactions between two variables.

The likelihood ratio test can be specified using the test argument to DESeq , which substitutes
nbinomWaldTest with nbinomLRT.  
In  this  case,  the  user  needs  to  provide  a  reduced  formula,  e.g.   one  in  which  a  number  of terms  from design(dds) are  removed.  
```

DESeq using LRT
========================================================


```r
ddsLRT<-DESeqDataSetFromMatrix(countData= fc$counts,colData= colData,design=~ Batch + Group)
      
ddsLRT <- DESeq(ddsLRT, test="LRT", full=~Batch+ Group, reduced=~Group)

resLRT<-results(ddsLRT)

resLRTorder<-resLRT[order(resLRT$padj),]
```

========================================================


```r
head(resLRTorder)
```

```
log2 fold change (MLE): Group Viv vs Hfd 
LRT p-value: '~ Batch + Group' vs '~ Group' 
DataFrame with 6 rows and 6 columns
                    baseMean log2FoldChange     lfcSE      stat
                   <numeric>      <numeric> <numeric> <numeric>
ENSMUSG00000067219  172.4854     -1.9971724 0.2802905  87.45952
ENSMUSG00000020122 5025.7129      1.3153495 0.1303957  51.01771
ENSMUSG00000041567 1349.3249     -1.4974541 0.1532013  50.92142
ENSMUSG00000049580 1487.8130     -0.6392233 0.1493971  45.58019
ENSMUSG00000019577  529.3305     -1.0547483 0.1903876  41.49217
ENSMUSG00000022181 3384.0488      1.1178976 0.1687633  38.04588
                         pvalue         padj
                      <numeric>    <numeric>
ENSMUSG00000067219 8.602240e-21 1.362939e-16
ENSMUSG00000020122 9.153612e-13 5.077393e-09
ENSMUSG00000041567 9.613846e-13 5.077393e-09
ENSMUSG00000049580 1.465163e-11 5.803509e-08
ENSMUSG00000019577 1.183441e-10 3.750087e-07
ENSMUSG00000022181 6.910046e-10 1.596878e-06
```

========================================================


```r
summary(resLRT)
```

```

out of 22605 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)     : 32, 0.14% 
LFC < 0 (down)   : 40, 0.18% 
outliers [1]     : 0, 0% 
low counts [2]   : 6761, 30% 
(mean count < 5)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```


Differential splicing analysis
========================================================


```r
 # Perform exon level read counting 

fcexo<-featureCounts(files=paste0(bamdir,targets$OutputFile),
 annot.ext=anno_for_featurecount,
 strandSpecific=strandspecific,
isGTFAnnotationFile=TRUE,
GTF.featureType="exon",
GTF.attrType="gene_id",
nthreads=8,
isPairedEnd=isPairedEnd,
useMetaFeatures=FALSE, 
allowMultiOverlap=TRUE)
```

========================================================


```r
  head(fcexo$counts)
```

```
                   X.Users.skhadaya.Downloads.rnaseqpractical.Viv1.bam
ENSMUSG00000090025                                                   0
ENSMUSG00000064842                                                   0
ENSMUSG00000051951                                                   0
ENSMUSG00000051951                                                   0
ENSMUSG00000051951                                                   0
ENSMUSG00000051951                                                   0
                   X.Users.skhadaya.Downloads.rnaseqpractical.Viv2.bam
ENSMUSG00000090025                                                   0
ENSMUSG00000064842                                                   0
ENSMUSG00000051951                                                   0
ENSMUSG00000051951                                                   0
ENSMUSG00000051951                                                   1
ENSMUSG00000051951                                                   1
                   X.Users.skhadaya.Downloads.rnaseqpractical.Viv3.bam
ENSMUSG00000090025                                                   0
ENSMUSG00000064842                                                   0
ENSMUSG00000051951                                                   0
ENSMUSG00000051951                                                   0
ENSMUSG00000051951                                                   0
ENSMUSG00000051951                                                   0
                   X.Users.skhadaya.Downloads.rnaseqpractical.Hfd1.bam
ENSMUSG00000090025                                                   0
ENSMUSG00000064842                                                   0
ENSMUSG00000051951                                                   1
ENSMUSG00000051951                                                   0
ENSMUSG00000051951                                                   1
ENSMUSG00000051951                                                   2
                   X.Users.skhadaya.Downloads.rnaseqpractical.Hfd2.bam
ENSMUSG00000090025                                                   0
ENSMUSG00000064842                                                   0
ENSMUSG00000051951                                                   0
ENSMUSG00000051951                                                   0
ENSMUSG00000051951                                                   0
ENSMUSG00000051951                                                   0
                   X.Users.skhadaya.Downloads.rnaseqpractical.Hfd3.bam
ENSMUSG00000090025                                                   0
ENSMUSG00000064842                                                   0
ENSMUSG00000051951                                                   0
ENSMUSG00000051951                                                   0
ENSMUSG00000051951                                                   0
ENSMUSG00000051951                                                   0
```


========================================================


```r
   # assemble all the counts into an edgeR DGEList object:
       
   dge <-DGEList(counts=fcexo$counts, genes=fcexo$annotation)

     dim(dge)
```

```
[1] 689492      6
```

========================================================


```r
  head(dge$genes) 
```

```
              GeneID Chr   Start     End Strand Length
1 ENSMUSG00000090025   1 3044314 3044814      +    501
2 ENSMUSG00000064842   1 3092097 3092206      +    110
3 ENSMUSG00000051951   1 3195982 3197398      -   1417
4 ENSMUSG00000051951   1 3196604 3197398      -    795
5 ENSMUSG00000051951   1 3203520 3205713      -   2194
6 ENSMUSG00000051951   1 3203690 3206425      -   2736
```

========================================================


```r
  # Filtering step  
  # Keep exons that have more than 10 cpm in total
  
        sumall <- rowSums(dge$counts)
        dge<-dge[sumall>10,,keep.lib.sizes=FALSE]
```

========================================================


```r
  # Apply TMM normalization
        
        dge <- calcNormFactors(dge)
                          
       dge$sample
```

```
                                                    group lib.size
X.Users.skhadaya.Downloads.rnaseqpractical.Viv1.bam     1 69920297
X.Users.skhadaya.Downloads.rnaseqpractical.Viv2.bam     1 46995421
X.Users.skhadaya.Downloads.rnaseqpractical.Viv3.bam     1 66420700
X.Users.skhadaya.Downloads.rnaseqpractical.Hfd1.bam     1 49782723
X.Users.skhadaya.Downloads.rnaseqpractical.Hfd2.bam     1 55146770
X.Users.skhadaya.Downloads.rnaseqpractical.Hfd3.bam     1 56267940
                                                    norm.factors
X.Users.skhadaya.Downloads.rnaseqpractical.Viv1.bam    1.0474196
X.Users.skhadaya.Downloads.rnaseqpractical.Viv2.bam    0.8981800
X.Users.skhadaya.Downloads.rnaseqpractical.Viv3.bam    0.9324318
X.Users.skhadaya.Downloads.rnaseqpractical.Hfd1.bam    1.0866839
X.Users.skhadaya.Downloads.rnaseqpractical.Hfd2.bam    0.9825045
X.Users.skhadaya.Downloads.rnaseqpractical.Hfd3.bam    1.0677291
```

========================================================


```r
#  A multi-dimensional scaling plot can be used to find the presence
# of outlier samples in our experiment.
       
  plotMDS(dge, labels=1:6, col=as.numeric(colData$Group), main="MDS plot")
  legend("topright", legend=c("Viv", "Hfd"), col=1:2, pch=15)
```

![plot of chunk unnamed-chunk-51](practical-figure/unnamed-chunk-51-1.png) 

========================================================


```r
  # Create design matrix
         design <- model.matrix(~targets$Group)
        design
         
  # Apply voom to convert the read counts to log2-cpm with associated weights:
         
        v <- voom(dge,design,plot=TRUE)
```


========================================================


```
  (Intercept) targets$GroupViv
1           1                1
2           1                1
3           1                1
4           1                0
5           1                0
6           1                0
attr(,"assign")
[1] 0 1
attr(,"contrasts")
attr(,"contrasts")$`targets$Group`
[1] "contr.treatment"
```

![plot of chunk unnamed-chunk-53](practical-figure/unnamed-chunk-53-1.png) 


========================================================


```r
  # Linear modelling
        
fit <- lmFit(v,design)
                          
  # Test for differences in exon retention
ex <- diffSplice(fit, geneid="GeneID")
```

```
Total number of exons:  419426 
Total number of genes:  18245 
Number of genes with 1 exon:  4317 
Mean number of exons in a gene:  23 
Max number of exons in a gene:  996 
```

========================================================

![plot of chunk unnamed-chunk-55](practical-figure/unnamed-chunk-55-1.png) 


========================================================

![plot of chunk unnamed-chunk-56](practical-figure/unnamed-chunk-56-1.png) 

========================================================


```r
   # Rank genes by evidence of differential splicing                       
    spliced<-topSplice(ex,coef=2,test="F")

   spliced
```

```
                   GeneID Chr Strand NExons         F      P.Value
626840 ENSMUSG00000032187   9      +    157  5.804282 2.260673e-57
199013 ENSMUSG00000022364  15      +    112  3.959148 2.515890e-25
533504 ENSMUSG00000038521   6      -     74  5.252492 4.075225e-25
547612 ENSMUSG00000002992   7      -     21 15.658949 1.610242e-21
121467 ENSMUSG00000020604  11      +     33  9.261986 1.658529e-21
291098 ENSMUSG00000075044  19      -     64  5.134996 2.032149e-21
103422 ENSMUSG00000017428  11      +    126  3.042660 1.442538e-18
105769 ENSMUSG00000059439  11      +    203  2.380457 1.277043e-17
258845 ENSMUSG00000079507  17      +     18 13.458592 6.083390e-17
191549 ENSMUSG00000022119  14      -    114  2.853708 3.600086e-15
                FDR
626840 3.148665e-53
199013 1.752066e-21
533504 1.891991e-21
547612 4.619998e-18
121467 4.619998e-18
291098 4.717295e-18
103422 2.870238e-15
105769 2.223331e-14
258845 9.414385e-14
191549 5.014199e-12
```

========================================================


```r
    #   Show individual exons that are enriched or depleted relative to other exons in the same gene.
    
  altUsed<-topSplice(ex,coef=2,test="t")
```

========================================================


```r
altUsed
```

```
                   GeneID Chr     Start       End Strand Length     logFC
547610 ENSMUSG00000002992   7  20266654  20266759      -    106  1.195129
547612 ENSMUSG00000002992   7  20266654  20266772      -    119  1.189237
322122 ENSMUSG00000035875   2  34956375  34956528      -    154  2.228447
533407 ENSMUSG00000038521   6 124483914 124483992      -     79  1.936305
533408 ENSMUSG00000038521   6 124483914 124483992      -     79  1.936305
533409 ENSMUSG00000038521   6 124483914 124483992      -     79  1.936305
533410 ENSMUSG00000038521   6 124483914 124483992      -     79  1.936305
533411 ENSMUSG00000038521   6 124483914 124483992      -     79  1.936305
322123 ENSMUSG00000035875   2  34956848  34957063      -    216  2.359795
139027 ENSMUSG00000021069  12  71331958  71332475      -    518 -2.456829
               t      P.Value          FDR
547610 10.339976 4.430733e-17 1.160100e-11
547612 10.291838 5.589377e-17 1.160100e-11
322122  8.546973 9.534507e-15 1.319287e-09
533407  7.735572 1.526915e-13 7.922954e-09
533408  7.735572 1.526915e-13 7.922954e-09
533409  7.735572 1.526915e-13 7.922954e-09
533410  7.735572 1.526915e-13 7.922954e-09
533411  7.735572 1.526915e-13 7.922954e-09
322123  7.540859 3.291195e-12 1.518005e-07
139027 -7.251057 7.278650e-12 3.021433e-07
```

========================================================


```r
   #  Save results
  
write.table(spliced,file="AltSplicedgenes.txt",sep="\t")
  
write.table(altUsed,file="DiffUsedExons.txt",sep="\t")
```

Gene Ontology and Pathway Enrichment Analysis
========================================================


```r
We will perform GO analysis using goseq package.

 In order to perform a GO analysis of your RNA-seq data,goseq only requires a simple named vector, which contains two pieces of information.

 1. Measured genes
 all genes for which RNA-seq data was gathered for your experiment.  Each element of your vector should be named by a unique gene identifier.

 2.Differentially expressed genes
 each element of your vector should be either a 1 or 0, where 1 indicates that the gene is differentially expressed and 0 that it is not.
```

========================================================


```r
 resdat<- resOrdered[complete.cases(resOrdered$padj),]
     degenes<-as.integer(resdat$padj<0.05)
    names(degenes)<-rownames(resdat)
    
    # remove duplicate gene names
  degenes<-degenes[match(unique(names(degenes)),names(degenes))]
  table(degenes)
```

```
degenes
    0     1 
15100  1993 
```

========================================================

```r
   # fitting the probability weighting function (PWF)
# We first need to obtain a weighting for each gene, depending on its length, given by the PWF
  
   pwf=nullp(degenes,genome,'ensGene',plot.fit = FALSE)
```

========================================================


```r
  head(pwf)
```

```
                   DEgenes bias.data       pwf
ENSMUSG00000024526       1      1118 0.1161615
ENSMUSG00000032080       1      1451 0.1257491
ENSMUSG00000026475       1      2341 0.1308108
ENSMUSG00000069170       1      7047 0.1308108
ENSMUSG00000042041       1      1167 0.1179737
ENSMUSG00000034634       1       682 0.0950137
```

========================================================


```r
   plotPWF(pwf)
```

![plot of chunk unnamed-chunk-65](practical-figure/unnamed-chunk-65-1.png) 

========================================================

```r
    # change the Keggpath id to name in the goseq output
                              
    xx <- as.list(KEGGPATHID2NAME)
    temp <- cbind(names(xx),unlist(xx))
    
    addKeggTogoseq <- function(JX,temp){
      for(l in 1:nrow(JX)){
     if(JX[l,1] %in% temp[,1]){
     JX[l,"term"] <- temp[temp[,1] %in% JX[l,1],2]
      JX[l,"ontology"] <- "KEGG"
                         }
                 }
      return(JX)
    }
```

========================================================

```r
#  Calculate  the  over  and  under  expressed  GO
# categories among DE genes

functional_analysis=goseq(pwf,genome,
'ensGene',test.cats=c("GO:BP","GO:MF","KEGG"))
```

========================================================

```r
head(functional_analysis)
```

```
        category over_represented_pvalue under_represented_pvalue
10400 GO:0044281            9.013811e-27                        1
2623  GO:0006082            2.234892e-26                        1
10082 GO:0043436            5.867345e-25                        1
6187  GO:0019752            4.825825e-24                        1
7781  GO:0032787            3.800220e-23                        1
3290  GO:0007155            2.030335e-20                        1
      numDEInCat numInCat                                  term ontology
10400        321     1461      small molecule metabolic process       BP
2623         209      801        organic acid metabolic process       BP
10082        203      787             oxoacid metabolic process       BP
6187         191      734     carboxylic acid metabolic process       BP
7781         140      475 monocarboxylic acid metabolic process       BP
3290         220      950                         cell adhesion       BP
```

========================================================

```r
restemp<-addKeggTogoseq(functional_analysis,temp)   

    head(restemp)
```

```
        category over_represented_pvalue under_represented_pvalue
10400 GO:0044281            9.013811e-27                        1
2623  GO:0006082            2.234892e-26                        1
10082 GO:0043436            5.867345e-25                        1
6187  GO:0019752            4.825825e-24                        1
7781  GO:0032787            3.800220e-23                        1
3290  GO:0007155            2.030335e-20                        1
      numDEInCat numInCat                                  term ontology
10400        321     1461      small molecule metabolic process       BP
2623         209      801        organic acid metabolic process       BP
10082        203      787             oxoacid metabolic process       BP
6187         191      734     carboxylic acid metabolic process       BP
7781         140      475 monocarboxylic acid metabolic process       BP
3290         220      950                         cell adhesion       BP
```


========================================================

```r
write.table(restemp,file="GO_Kegg_Wallenius.txt", row.names=F,sep="\t")
```

Session Information
========================================================

```r
    sessionInfo()
```

