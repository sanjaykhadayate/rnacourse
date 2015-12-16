# load libraries

library(Rsubread)
library(limma)
library(edgeR)
library(DESeq2)
library(goseq)
library(RColorBrewer)
library(ggplot2)
library(KEGG.db)
library(org.Mm.eg.db)
library(gplots)

#  Some useful variables
genome<-"mm9" 
strandspecific<-0 
factor1<-"Group"
dir="/home/ubuntu/referenceData"
isPairedEnd<-TRUE
nthreads=4
bamdir="/home/ubuntu/referenceData"

# read in target file using a funcion from limma package
targets <- readTargets("targets.txt")
targets

# Get full path names of raw data files
read1=paste0(dir,targets$InputFile[6])
read2=paste0(dir,targets$InputFile2[6])

# get the name of rsubread genome index

index="mm9_index"

#Step 1: Index building
ref <- system.file("extdata","reference.fa",package="Rsubread")
buildindex(basename="reference_index",reference=ref)

# align reads using rsubreads package
if (isPairedEnd == TRUE)
{
  subjunc(index=index,readfile1=read1,readfile2=read2,input_format="gzFASTQ",output_format="BAM",output_file=paste0(bamdir,targets$OutputFile),nthreads=nthreads,unique=TRUE,indels=5)
  
}else {
  subjunc(index=index,readfile1=read1,input_format="gzFASTQ",output_format="BAM",output_file=paste0(bamdir,targets$OutputFile),nthreads=nthreads,unique=TRUE,indels=5)
}


# Base quality scores
qscore <- qualityScores(filename=read1[1],offset=33,nreads=1000)
plot(colMeans(qscore),ylim=c(0,40))

# GC content

atgc<-atgcContent(filename=read1[1],basewise = TRUE)

plot(atgc[1,],col="red",ylim=c(0,0.4),ylab="content")
points(atgc[2,],col="green")
points(atgc[3,],col="yellow")
points(atgc[4,],col="brown")
legend("topright", legend=c("A", "T","G","C"), col=c("red","green","yellow","brown"), pch=15)


# produce alignment statistics
alignstats<-propmapped(targets$OutputFile[1])
write.table(alignstats,file="AlignmentSummary.txt",sep="\t")

# count numbers of reads mapped to iGenome genes
anno_for_featurecount<-paste0(bamdir,"genes.gtf")

  fc <-featureCounts(files=paste0(bamdir,targets$OutputFile),annot.ext=anno_for_featurecount,isGTFAnnotationFile=TRUE,
   GTF.featureType="exon",useMetaFeatures=TRUE, GTF.attrType="gene_id",nthreads=nthreads,strandSpecific=strandspecific,isPairedEnd=isPairedEnd)
                              
  head(fc$counts)
  
  head(fc$annotation)
  
  write.table(fc$counts,file="genecounts.txt",sep="\t")
  
  # collect sample information
                              
 colData<-cbind(targets$OutputFile,targets$Group,targets$Batch)
 
 rownames(colData)<-colData[,1]
colnames(colData)<-c("name","Group","Batch")
                             
 colData<-data.frame(colData)
 colData
 
# construct deseqdataset object
  dds<-DESeqDataSetFromMatrix(countData= fc$counts,colData= targets,design=~Group)
                            
# Perform normalization, fitting to the model
    dds<-DESeq(dds)
    
    sizeFactors(dds)
    
    head(dispersions(dds))
    
    plotDispEsts(dds)
    
    # Some quality control  #####
      # Do a PCA plot
    rld<-rlog(dds)  
    plotPCA(rld, intgroup="Group")
    
      
  #  Sample clustering plot
         rlogcount <- assay(rld)
         rlogcount <- rlogcount[!rowSums(rlogcount) == 0,]
         colnames(rlogcount) <-  paste0(colData(dds)$Sample)                          
     sampleDists <- as.matrix(dist(t(rlogcount)))
     showcol <- brewer.pal(8, "Dark2")[1:length(unique(colData(dds)$Group))]
   
    heatmap.2(as.matrix(sampleDists), key=F, trace="none",
  col=colorpanel(100, "black", "white"),
  ColSideColors=showcol[colData(dds)$Group], RowSideColors=showcol[colData(dds)$Group],
                     margin=c(10, 10), main="Sample Distance Matrix")
        
                                                          
    ## Getting and saving results
                              
    res<-results(dds)
    
     summary(res)
    
      # How many adjusted p-values were less than 0.1?
     
        sum(res$padj < 0.1, na.rm=TRUE)   
        
        resOrdered<-res[order(res$padj),]
    
      write.table(resOrdered,file="DEgenes.txt",sep="\t")
    
  # MA plot
      plotMA(res, main="DESeq2", ylim=c(-4,4))
      
  # plot counts
      plotCounts(dds, gene=which.min(res$padj), intgroup="Group")
      
      
  # Multi Factor Designs
      
      ddsMF<-DESeqDataSetFromMatrix(countData= fc$counts,colData= colData,design=~ Batch + Group)
      
      
      ddsMF <- DESeq(ddsMF)
      
      resMF <- results(ddsMF)
      
      resMForder<-resMF[order(resMF$padj),]
      head(resMForder)
    
      
      resMFbatch <- results(ddsMF, contrast=c("Batch","a","b"))
      
      resMFbatchOrder<-resMFbatch[order(resMFbatch$padj),]
      
      head(resMFbatchOrder)
      
      summary(resMFbatch)
      
      
      
      ddsLRT<-DESeqDataSetFromMatrix(countData= fc$counts,colData= colData,design=~ Batch + Group)
      
      ddsLRT <- DESeq(ddsLRT, test="LRT", full=~Batch+ Group, reduced=~Group)
      
      resLRT<-results(ddsLRT)
      
      resLRTorder<-resLRT[order(resLRT$padj),]
      
      head(resLRTorder)
      
      summary(resLRT)
      
      
      
  # Perform splicing analysis

 # Perform exon level read counting 
   fcexo<-featureCounts(files=paste0(bamdir,targets$OutputFile),annot.ext=anno_for_featurecount,strandSpecific=strandspecific,
     isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_id",nthreads=8,
      isPairedEnd=isPairedEnd,useMetaFeatures=FALSE, allowMultiOverlap=TRUE)
       
   head(fcexo$counts)
   
   write.table(fcexo$counts,file="exoncounts.txt",sep="\t")
   
   # Get DEGList object
       dge <- DGEList(counts=fcexo$counts, genes=fcexo$annotation)
       dim(dge)
       
       head(dge$genes)
       
  # Filtering step     
        sumall <- rowSums(dge$counts)
        dge <- dge[sumall>10,,keep.lib.sizes=FALSE]
        dim(dge)
        
  # Calculate normalization factors and perform normalization
        dge <- calcNormFactors(dge)
               
        dge$sample
        
        plotMDS(dge, labels=1:6, col=as.numeric(colData$Group), main="MDS plot")
        legend("topright", legend=c("Viv", "Hfd"), col=1:2, pch=15)
        
        
         design <- model.matrix(~targets$Group)
        voomnorm <- voom(dge,design,plot=TRUE)
        
  # Fit linear model       
         fit <- lmFit(voomnorm,design)
                          
  # Test for differences in exon retention
        exondiff <- diffSplice(fit, geneid="GeneID")
  
        
        plotSplice(exondiff, geneid="ENSMUSG00000002992", genecol="GeneID")
        
       # Rank genes by evidence of differential splicing                       
    spliced<-topSplice(exondiff,coef=2,test="F")
             
    #   Show individual exons that are enriched or depleted relative to other exons in the 
    # same gene.
    
  altUsed<-topSplice(exondiff,coef=2,test="t")
  
  
  
   #  Save results
  
    write.table(spliced,file="AltSplicedgenes.txt",sep="\t")
    write.table(altUsed,file="DiffUsedExons.txt",sep="\t")
                              
                             
     # Gene Ontology and pathway enrichment analysis
                             
                              
    resOrdered<- resOrdered[complete.cases(resOrdered$padj),]
     degenes<-as.integer(resOrdered$padj<0.05)
    names(degenes)<-rownames(resOrdered)
    
    # remove duplicate gene names
  degenes<-degenes[match(unique(names(degenes)),names(degenes))]
  table(degenes)
  
                         
   # fitting the probability weighting function (PWF)
   pwf=nullp(degenes,genome,'ensGene')
   
   head(pwf)
   
   plotPWF(pwf)
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
    
    functional_analysis=goseq(pwf,genome,'ensGene',test.cats=c("GO:BP","GO:MF","KEGG"))
   
    head(functional_analysis) 
    
    restemp<-addKeggTogoseq(functional_analysis,temp)    # switch Keggpathid to name
    
    head(restemp)
    
    write.table(restemp,file="GO_Kegg_Wallenius.txt",row.names=F,sep="\t")
                              
                              
    sessionInfo()
                              
