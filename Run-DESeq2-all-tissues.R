### Run DESeq2 - separate groups

library(tidyverse)
library(DESeq2)
library(edgeR)
library(sva)

setwd("~/projects/Kiskinis/Target ALS/scripts/")


all.md<-read.table("../output/DESeq2/2023/2023-metadata-v2-fixed.txt", sep="\t", header=T, stringsAsFactors = F, quote="")
all.md$ExternalSampleId<-str_replace_all(all.md$ExternalSampleId, "-", ".")
### Fix "ALS/FTD" group
#all.md$Subject.Group2<-ifelse(all.md$Subject.Group2=="ALS/FTD", "ALS.FTD", all.md$Subject.Group2)
### Set Healthy as reference
all.md$Subject.Group2<-as.factor(all.md$Subject.Group2)
all.md$Subject.Group2<-relevel(all.md$Subject.Group2, ref = "Control")


### Total patients
#nodup<-all.md[,c("ExternalSubjectId", "mutation3")]
#nodup<-nodup[!duplicated(nodup),]
#table(nodup$mutation3)

all.md<-all.md[which(all.md$mutation!="NEK1"),]
all.md<-all.md[which(all.md$mutation!="FUS"),]
### Remove age==unknown
all.md<-all.md[which(all.md$Age.at.Death!="Unknown"),]
### Remove age==Not Applicable (is this person still alive?!)
all.md<-all.md[which(all.md$Age.at.Death!="Not Applicable"),]
all.md<-all.md[which(all.md$RIN!="na"),]
all.md<-all.md[which(all.md$RIN!="NA "),]


all.md$RIN<-as.numeric(all.md$RIN)

write.table(all.md, "../output/DESeq2/2023/2023-metadata-v2-fixed-filtered.txt", sep="\t", col.names = T, row.names = F, quote=F)

### Read in counts

allcounts<-read.table("../data/counts/withNYGCsamples/2023-counts-CNS-tissues-v2.txt", sep="\t", header=T, stringsAsFactors = F, quote="")
allcounts$external_gene_name<-ifelse(allcounts$external_gene_name=="", allcounts$ensembl_gene_id, allcounts$external_gene_name)
rownames(allcounts)<-allcounts$ensembl_gene_id

geneinfo<-allcounts[,1:3]

### ALS/ALS-FTD/FTD as broad groups
### ALS/ALS-FTD/FTD breakdown by tissue
table(all.md$Sample.group3, all.md$Subject.Group2)
### FTD = don't run for Occipital (1 sample), Cortex.Motor (2 samples), Hippocampus (0 samples). Keep in Spinal (3 samples)
### ALS/FTD - don't run for Temporal (1 sample)

### Samples for each
# Cerebellum - all
# Cortex_Frontal - all
# Cortex_Occipital - ALS, ALS/FTD, Control
# Cortex_Temporal - ALS, Control, FTD
# Cortex.Motor - ALS, ALS/FTD, Control
# Hippocampus - ALS, ALS/FTD, Control
# Spinal - all

### Run DE in loop
tissuesofinterest<-unique(all.md$Sample.group3)

svList<-list()

for(i in tissuesofinterest){
  tissue=i
  exclude.groups<-NULL
  
  cat(paste0("Running tissue: ", tissue, "\n"))
  
  ### Get metadata for this tissue
  tiss.md<-all.md[which(all.md$Sample.group3==tissue),]
  
  ### Get non-control groups
  subjectgroups<-unique(tiss.md$Subject.Group2)
  subjectgroups<-subjectgroups[subjectgroups!="Control"]
  
  ### Take only samples from this tissue
  #tiss.md<-all.md[which(all.md$Sample.group3==tissue),]
   
  ### If a subject group has less than 3 samples, exclude it 
  exclude.groups=names(table(tiss.md$Subject.Group2))[which(table(tiss.md$Subject.Group2)<3)]
  
  if(length(exclude.groups)>0){
    cat(paste0("Excluding ", exclude.groups, " from ", tissue, " analysis. Fewer than 3 samples.\n"))
  }
  
  ### If excluding groups with too few samples, remove from md
  tiss.md<-tiss.md[!(tiss.md$Subject.Group2 %in% exclude.groups),]
 
  ### Remove age==unknown
  tiss.md<-tiss.md[which(tiss.md$Age.at.Death!="Unknown"),]
  ### Remove age==Not Applicable (is this person still alive?!)
  tiss.md<-tiss.md[which(tiss.md$Age.at.Death!="Not Applicable"),]
  
  ### Set sex as a factor
  tiss.md$Sex<-as.factor(tiss.md$Sex)
  #tiss.md$Prep<-as.factor(tiss.md$Prep)
  
  ### Set age at death as a numeric (also change "90 or Older" to 90)
  tiss.md$Age.at.Death<-as.numeric(ifelse(tiss.md$Age.at.Death=="90 or Older", 90, tiss.md$Age.at.Death))
  
  ### Set RIN as numeric
  tiss.md$RIN<-as.numeric(tiss.md$RIN)
  ### Drop empty levels
  tiss.md$Subject.Group2<-as.factor(as.character(tiss.md$Subject.Group2))
  tiss.md$Subject.Group2<-relevel(tiss.md$Subject.Group2, ref = "Control")
  
  
  ### Take only samples from this tissue
  tiss.counts<-allcounts[,colnames(allcounts) %in% tiss.md$ExternalSampleId]
  
  nSamp<-nrow(tiss.md)
  minSamp<-table(tiss.md$Subject.Group2)[which(table(tiss.md$Subject.Group2)==min(table(tiss.md$Subject.Group2)))]
  
  ### Make sure order of samples in metadata matches order of samples in counts
  all(tiss.md$ExternalSampleId==colnames(tiss.counts))
  tiss.counts<-tiss.counts[,match(tiss.md$ExternalSampleId, colnames(tiss.counts))]
  all(tiss.md$ExternalSampleId==colnames(tiss.counts))
  
  ### Remove low expression genes (although not necessary for DESeq2)
  keep<-rowSums(tiss.counts>5)>=ceiling(nSamp/2)
  #keep<-rowSums(tiss.counts>5)>=minSamp
  
  #table(keep)
  tiss.filt<-tiss.counts[keep,]
  
  
  
  ### normalize and run vst for svaseq
  dds <- DESeqDataSetFromMatrix(countData = tiss.filt, colData = tiss.md, design = ~1)
  dds <- DESeq(dds)
  
  vsd<-vst(dds)
  vsd.out<-assay(vsd)
  
  vsd.out2<-vsd.out %>% as.data.frame %>% rownames_to_column("Gene")
  vsd.out2<-merge(geneinfo, vsd.out2, by=1)
  
  outname=paste0("../output/DESeq2/2023/sva-vsd/noPrep/disease-groups/vsd/", tissue, "-AllGroups.vs.Healthy.txt")
  write.table(vsd.out2, file = outname, sep = "\t", col.names = T, row.names = F, quote=F)
  
  ### Get CPM
  y<-DGEList(counts = tiss.filt)
  y<-calcNormFactors(y)
  cpmOut<-as.data.frame(cpm(y))
  cpmOut<-merge(geneinfo, cpmOut, by.x=1, by.y="row.names")
  outname=paste0("../output/DESeq2/2023/sva-vsd/noPrep/disease-groups/TMMCPM/", tissue, "-AllGroups.vs.Healthy.txt")
  write.table(cpmOut, file = outname, sep = "\t", col.names = T, row.names = F, quote=F)
  
  ### Run svaseq
  ### Don't include Prep in design
  mod<-model.matrix(~tiss.md$Subject.Group2 + tiss.md$Sex + tiss.md$RIN + tiss.md$Age.at.Death)
  mod0<-model.matrix(~tiss.md$Sex + tiss.md$RIN + tiss.md$Age.at.Death)
  
  set.seed(42)
  svseq <- svaseq(vsd.out, mod = mod, mod0 = mod0)
  
  ### Save SVs for future reference
  sv.out<-as.data.frame(svseq$sv)
  colnames(sv.out)<-paste0("SV", 1:ncol(sv.out))
  rownames(sv.out)<-tiss.md$ExternalSampleId
  sv.out %>% rownames_to_column("Sample") -> sv.out
  outname<-paste0("../output/DESeq2/2023/sva-vsd/noPrep/disease-groups/SVs/", tissue, "-AllGroups-SVs.txt")
  write.table(sv.out, file = outname, sep = "\t", col.names = T, row.names = F, quote=F)
  
  svtmp<-data.frame(tissue=tissue, numSVs=ncol(svseq$sv))
  svList[[i]]<-svtmp
  
  tiss.md2<-cbind(tiss.md, svseq$sv)
  mdname=paste0("../output/DESeq2/2023/sva-vsd/noPrep/disease-groups/metadata/", tissue, "-AllGroups-metadata.txt")
  write.table(tiss.md2, file = mdname, sep="\t", col.names = T, row.names = F, quote=F)
  
  ### Create design matrix with mutation, sex, prep, rin, age, and SVs
  design<-model.matrix(~tiss.md$Subject.Group2 + tiss.md$Sex + tiss.md$RIN + tiss.md$Age.at.Death + svseq$sv)
  #colnames(design)
  
  ### Create DESeq object with design matrix and run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = tiss.filt, colData = tiss.md, design = design)
  dds <- DESeq(dds)
  
  all(colnames(tiss.filt)==tiss.md2$ExternalSampleId)
  #resultsNames(dds)
  
  ### Save ALS 
  if(!("ALS" %in% exclude.groups)){
    sva.out<-as.data.frame(results(dds, name = "tiss.md.Subject.Group2ALS"))
    sva.out %>% rownames_to_column("GeneID") -> sva.out
    sva.out<-merge(geneinfo, sva.out, by=1)
    
    outname=paste0("../output/DESeq2/2023/sva-vsd/noPrep/disease-groups/", tissue, "-ALS.vs.Healthy.txt")
    write.table(sva.out, file = outname, sep = "\t", col.names = T, row.names = F, quote=F)
    
  }
  
  ### Save ALSFTD
  if(!("ALSFTD" %in% exclude.groups)){
    sva.out<-as.data.frame(results(dds, name = "tiss.md.Subject.Group2ALSFTD"))
    sva.out %>% rownames_to_column("GeneID") -> sva.out
    sva.out<-merge(geneinfo, sva.out, by=1)
    
    outname=paste0("../output/DESeq2/2023/sva-vsd/noPrep/disease-groups/", tissue, "-ALSFTD.vs.Healthy.txt")
    write.table(sva.out, file = outname, sep = "\t", col.names = T, row.names = F, quote=F)
    
  }
  
  ### Save ALSFTD
  if(!("FTD" %in% exclude.groups)){
    sva.out<-as.data.frame(results(dds, name = "tiss.md.Subject.Group2FTD"))
    sva.out %>% rownames_to_column("GeneID") -> sva.out
    sva.out<-merge(geneinfo, sva.out, by=1)
    
    outname=paste0("../output/DESeq2/2023/sva-vsd/noPrep/disease-groups/", tissue, "-FTD.vs.Healthy.txt")
    write.table(sva.out, file = outname, sep = "\t", col.names = T, row.names = F, quote=F)
    
  }
  
  cat(paste0("Done with ", tissue, "\n"))
  
}

### Individual spinal sections
### Run DE in loop
tissuesofinterest<-c("Spinal_Cord_Cervical", "Spinal_Cord_Thoracic", "Spinal_Cord_Lumbar")
svList<-list()

for(i in tissuesofinterest){
  
  tissue=i
  
  cat(paste0("Running tissue: ", tissue, "\n"))
  
  tiss.md<-all.md[which(all.md$Sample.group2==tissue),]

    ### Take only samples from this tissue
    tiss.md<-all.md[which(all.md$Sample.group2==tissue),]
    
    exclude.groups=names(table(tiss.md$Subject.Group2))[which(table(tiss.md$Subject.Group2)<3)]
    
    if(length(exclude.groups)>0){
      cat(paste0("Excluding ", exclude.groups, " from ", tissue, " analysis. Fewer than 3 samples.\n"))
    }
    
    tiss.md<-tiss.md[!(tiss.md$Subject.Group2 %in% exclude.groups),]
    
    ### Remove age==unknown
    tiss.md<-tiss.md[which(tiss.md$Age.at.Death!="Unknown"),]
    ### Remove age==Not Applicable (is this person still alive?!)
    tiss.md<-tiss.md[which(tiss.md$Age.at.Death!="Not Applicable"),]
    
    tiss.md$RIN<-as.numeric(tiss.md$RIN)
    tiss.md$Sex<-as.factor(tiss.md$Sex)
    tiss.md$Prep<-as.factor(tiss.md$Prep)
    tiss.md$Age.at.Death<-as.numeric(ifelse(tiss.md$Age.at.Death=="90 or Older", 90, tiss.md$Age.at.Death))
    ### Drop empty levels
    tiss.md$Subject.Group2<-as.factor(as.character(tiss.md$Subject.Group2))
    tiss.md$Subject.Group2<-relevel(tiss.md$Subject.Group2, ref = "Control")
    
    ### Take only samples from this tissue
    tiss.counts<-allcounts[,colnames(allcounts) %in% tiss.md$ExternalSampleId]
    
    nSamp<-nrow(tiss.md)
    
    ### Make sure order of samples in metadata matches order of samples in counts
    all(tiss.md$ExternalSampleId==colnames(tiss.counts))
    tiss.counts<-tiss.counts[,match(tiss.md$ExternalSampleId, colnames(tiss.counts))]
    all(tiss.md$ExternalSampleId==colnames(tiss.counts))
    
    ### Remove low expression genes (although not necessary for DESeq2)
    keep<-rowSums(tiss.counts>5)>=ceiling(nSamp/2)
    #table(keep)
    tiss.filt<-tiss.counts[keep,]
    
    
    
    ### normalize and run vst for svaseq
    dds <- DESeqDataSetFromMatrix(countData = tiss.filt, colData = tiss.md, design = ~1)
    dds <- DESeq(dds)
    
    vsd<-vst(dds)
    vsd.out<-assay(vsd)
    
    vsd.out2<-vsd.out %>% as.data.frame %>% rownames_to_column("Gene")
    vsd.out2<-merge(geneinfo, vsd.out2, by=1)
    
    outname=paste0("../output/DESeq2/2023/sva-vsd/noPrep/disease-groups/vsd/", tissue, "-AllGroups.vs.Healthy.txt")
    write.table(vsd.out2, file = outname, sep = "\t", col.names = T, row.names = F, quote=F)
    
    ### Get CPM
    y<-DGEList(counts = tiss.filt)
    y<-calcNormFactors(y)
    cpmOut<-as.data.frame(cpm(y))
    cpmOut<-merge(geneinfo, cpmOut, by.x=1, by.y="row.names")
    outname=paste0("../output/DESeq2/2023/sva-vsd/noPrep/disease-groups/TMMCPM/", tissue, "-AllGroups.vs.Healthy.txt")
    write.table(cpmOut, file = outname, sep = "\t", col.names = T, row.names = F, quote=F)
    
    ### Run svaseq
    ### Don't include Prep in design
    mod<-model.matrix(~tiss.md$Subject.Group2 + tiss.md$Sex + tiss.md$RIN + tiss.md$Age.at.Death)
    mod0<-model.matrix(~tiss.md$Sex + tiss.md$RIN + tiss.md$Age.at.Death)
    
    set.seed(42)
    svseq <- svaseq(vsd.out, mod = mod, mod0 = mod0)
    
    ### Save SVs for future reference
    sv.out<-as.data.frame(svseq$sv)
    colnames(sv.out)<-paste0("SV", 1:ncol(sv.out))
    rownames(sv.out)<-tiss.md$ExternalSampleId
    sv.out %>% rownames_to_column("Sample") -> sv.out
    outname<-paste0("../output/DESeq2/2023/sva-vsd/noPrep/disease-groups/SVs/", tissue, "-AllGroups-SVs.txt")
    write.table(sv.out, file = outname, sep = "\t", col.names = T, row.names = F, quote=F)
    
    svtmp<-data.frame(tissue=tissue, numSVs=ncol(svseq$sv))
    svList[[i]]<-svtmp
    
    tiss.md2<-cbind(tiss.md, svseq$sv)
    mdname=paste0("../output/DESeq2/2023/sva-vsd/noPrep/disease-groups/metadata/", tissue, "-AllGroups-metadata.txt")
    write.table(tiss.md2, file = mdname, sep="\t", col.names = T, row.names = F, quote=F)
    
    
    ### Create design matrix with mutation, sex, prep, rin, age, and SVs
    design<-model.matrix(~tiss.md$Subject.Group2 + tiss.md$Sex + tiss.md$RIN + tiss.md$Age.at.Death + svseq$sv)
    #colnames(design)
    
    ### Create DESeq object with design matrix and run DESeq2
    dds <- DESeqDataSetFromMatrix(countData = tiss.filt, colData = tiss.md, design = design)
    dds <- DESeq(dds)
    #resultsNames(dds)
    
    ### Save ALS 
    if(!("ALS" %in% exclude.groups)){
      sva.out<-as.data.frame(results(dds, name = "tiss.md.Subject.Group2ALS"))
      sva.out %>% rownames_to_column("GeneID") -> sva.out
      sva.out<-merge(geneinfo, sva.out, by=1)
      
      outname=paste0("../output/DESeq2/2023/sva-vsd/noPrep/disease-groups/", tissue, "-ALS.vs.Healthy.txt")
      write.table(sva.out, file = outname, sep = "\t", col.names = T, row.names = F, quote=F)
      
    }
    
    ### Save ALSFTD
    if(!("ALSFTD" %in% exclude.groups)){
      sva.out<-as.data.frame(results(dds, name = "tiss.md.Subject.Group2ALSFTD"))
      sva.out %>% rownames_to_column("GeneID") -> sva.out
      sva.out<-merge(geneinfo, sva.out, by=1)
      
      outname=paste0("../output/DESeq2/2023/sva-vsd/noPrep/disease-groups/", tissue, "-ALSFTD.vs.Healthy.txt")
      write.table(sva.out, file = outname, sep = "\t", col.names = T, row.names = F, quote=F)
      
    }
    
    ### Save ALSFTD
    if(!("FTD" %in% exclude.groups)){
      sva.out<-as.data.frame(results(dds, name = "tiss.md.Subject.Group2FTD"))
      sva.out %>% rownames_to_column("GeneID") -> sva.out
      sva.out<-merge(geneinfo, sva.out, by=1)
      
      outname=paste0("../output/DESeq2/2023/sva-vsd/noPrep/disease-groups/", tissue, "-FTD.vs.Healthy.txt")
      write.table(sva.out, file = outname, sep = "\t", col.names = T, row.names = F, quote=F)
      
    }
  cat(paste0("Done with ", tissue, "\n"))
  
}
