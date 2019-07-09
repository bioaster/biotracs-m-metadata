### PERFORM BATCHING FOR MOSAIC MOUSE STUDY
### VARIABLES TO SHUFFLE : treatment AND TIME POINT
###
### BE CAREFUL FOR OTHER STUDIES: NO LONGITUDINAL FOLLOW-UP FOR THE MOUSES

rm(list=ls())
path="U:/OMICS/Lyon/MetaboProteo/MOSAIC_BA020/BA020-WP99-Omics-01/Data/Metadata/Mouse/"

### 1) LOAD LIBRARY
library(xlsx)
library(gridExtra)
library(grid)
library(lattice)
library(ggplot2)
library(MASS)
library(reshape)
library(parallel)

### 2) LOAD DATA
# covdesc = read.xlsx(paste0(path,"MetaData_v4.csv"),sheetIndex=5,stringsAsFactors=F)
covdesc = read.csv(paste0(path,"MetaData_v4.csv"),stringsAsFactors=F)
colnames(covdesc) =c("subject_index","	noe_id_sample","partner_subject_id","noe_id_mouse","omics_radomization","technical_id","technical_randmization","MouseNumber","batch_library_transcripto","batch_sequencing_miRNA","batch_sequencing_mRNA","sample_type"	,"day","day_alt2","day_alt3","treatment","treatment_alt2","treatment_alt3","treatment_alt4","taxonomy_class","taxonomy_species","partner_room","partner_intra_study_id","partner_variable_alt1","partner_variable_alt2","Comment","SamplingDateTime","PlasmaVolume")

chisq.test( table(covdesc$day,covdesc$treatment )) # pas de difference mais c'est normal, on a toutes les donnees

# on garde les echantillons communs transcriptos et metabo
# idx_MT=which((covdesc$validtranscripto=="oui")*(covdesc$validmetabo=="oui")==1)
# table_transcripto_metabo=table(covdesc$day[idx_MT],covdesc$treatment[idx_MT])

# chisq.test(table_transcripto_metabo) # pas de difference
# print(table_transcripto_metabo)
# on rajoute un echantillon traitement B à D+7 et un echantillon B à D+1
# set.seed(1234)
# idx1=sample(which( ((covdesc$day=="D+7")*(covdesc$treatment=="B")*(covdesc$validtranscripto=="oui")*(covdesc$validmetabo=="non"))==1),size=1)
# idx2=which( ((covdesc$day=="D+1")*(covdesc$treatment=="B")*(covdesc$validtranscripto=="oui")*(covdesc$validmetabo=="non"))==1)

# covdesc=covdesc[c(idx_MT,idx1,idx2),]

### 3) COMPUTE OBSERVED FREQUENCIES
obsFreq = lapply(c("noe_id_mouse","day","treatment"), function(cov) table(covdesc[,cov])/nrow(covdesc))
names(obsFreq) = c("noe_id_mouse","day","treatment")

### 4) VAR DEFINITION
#timePoints = unique(covdesc$day)
nbBatches = 3
batchSize = 66 # libraries
num_sequencing=c(rep(1:12,each=15),rep(13,12)) # sequencing
num_sequencing_miRNA=c(rep(1:2,each=39),rep(3:5,each=38))

pvThres = 0.9 # on considere que c'est bien reparti si toutes les pval sont superieures

### 5) PERFORM BATCHING
equalobsFreq = F
while(!equalobsFreq){
  
  ibatches = NULL    # list with index per batch
  i2batch = rep(0,nrow(covdesc)) # vector with batch number for each sample
  for(i in 1:nbBatches){
    i2sampleFrom = setdiff(1:nrow(covdesc),unlist(ibatches))
    if(length(i2sampleFrom)>=batchSize){
      ibatches[[i]] = sample(i2sampleFrom,batchSize)
    }#else{ibatches[[i]] = i2sampleFrom}
    i2batch[ibatches[[i]]] = i
  }
  
  tableday=sapply(ibatches,function(x) table(c(names(obsFreq[["day"]]),covdesc[x,"day"]))-1)
  tabletreatment=sapply(ibatches,function(x) table(c(names(obsFreq[["treatment"]]),covdesc[x,"treatment"]))-1)
  # 
  # ibatches_sequencing=lapply(unique(num_sequencing),function(x){unlist(ibatches)[num_sequencing==x]})
  # tableday_seq = sapply(ibatches_sequencing,function(x) table(c(names(obsFreq[["day"]]),covdesc[x,"day"]))-1)
  # tabletreatment_seq = sapply(ibatches_sequencing,function(x) table(c(names(obsFreq[["treatment"]]),covdesc[x,"treatment"]))-1)
  # # 
  # ibatches_sequencing_miRNA=lapply(unique(num_sequencing_miRNA),function(x){unlist(ibatches)[num_sequencing_miRNA==x]})
  # tableday_seq_miRNA = sapply(ibatches_sequencing_miRNA,function(x) table(c(names(obsFreq[["day"]]),covdesc[x,"day"]))-1)
  # tabletreatment_seq_miRNA = sapply(ibatches_sequencing_miRNA,function(x) table(c(names(obsFreq[["treatment"]]),covdesc[x,"treatment"]))-1)
  
  
  #       if(all(tableday_seq>0) )      # chaque batch contient tous les day
  #               if(all(tabletreatment_seq>0) )  # chaque batch contient tous les traitements
  #               {
  pvalues = c(
    #fisher.test(sapply(ibatches,function(x) table(c(names(obsFreq[["mouseID"]]),covdesc[x,"mouseID"]))-1),simulate.p.value=TRUE)$p.value,
    fisher.test(tableday,simulate.p.value=TRUE)$p.value
    ,fisher.test(tabletreatment,simulate.p.value=TRUE)$p.value
    #,fisher.test(sapply(ibatches,function(x) table(c(names(obsFreq[["Sexe"]]),covdesc[x,"Sexe"]))-1))$p.value
    #,kruskal.test(covdesc$Age, as.factor(i2batch))$p.value)
  )
  print(pvalues)
  
  if(all(pvalues>pvThres)){  # if library batches well balanced ==> sequencing
    pvalues_sequencing = c(
      #fisher.test(sapply(ibatches_sequencing,function(x) table(c(names(obsFreq[["mouseID"]]),covdesc[x,"mouseID"]))-1),simulate.p.value=TRUE)$p.value,
      fisher.test(tableday_seq,simulate.p.value=TRUE)$p.value
      ,fisher.test(tabletreatment_seq,simulate.p.value=TRUE)$p.value)
    equalobsFreq=T
    print(pvalues_sequencing)
    if(all(pvalues_sequencing>pvThres)){
      pvalues_sequencing_miRNA = c(
        #fisher.test(sapply(ibatches_sequencing_miRNA,function(x) table(c(names(obsFreq[["mouseID"]]),covdesc[x,"mouseID"]))-1),simulate.p.value=TRUE)$p.value,
        fisher.test(tableday_seq_miRNA ,simulate.p.value=TRUE)$p.value
        ,fisher.test(tabletreatment_seq_miRNA,simulate.p.value=TRUE)$p.value)
      print(pvalues_sequencing_miRNA)
      if(all(pvalues_sequencing_miRNA>pvThres)) {
        equalobsFreq=T}
    }
  }
  # }
}

### 6) PLOT COVARIATES DISTRIBUTION IN BATCHES (library)
df = lapply(c("noe_id_mouse","day","treatment"), function(cov){
data.frame(rbind(t(sapply(1:nbBatches, function(i) (table(c(names(obsFreq[[cov]]),covdesc[ibatches[[i]],cov]))-1)/table(covdesc[,cov])))
), c(paste("Batch",1:nbBatches,sep="")))
})

for(i in 1:length(df)){colnames(df[[i]]) = c(names(obsFreq[[i]]),"Batch")}
df.final = lapply(df,melt)
for(i in 1:length(df.final)){colnames(df.final[[i]]) = c("Batch",c("noe_id_mouse","day","treatment")[i],"Frequency")}

#p1 <- ggplot(data=df.final[[1]], aes(fill=mouseID, y=Frequency, x=Batch)) + geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p2 <- ggplot(data=df.final[[2]], aes(fill=Batch, y=Frequency, x=day)) + geom_bar(stat="identity") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p3 <- ggplot(data=df.final[[3]], aes(fill=Batch, y=Frequency, x=treatment)) + geom_bar(stat="identity") + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
#grid.arrange(p1,p2,p3,  ncol=3)

#pdf(paste0(path,"mouse_batching_library.pdf"))
grid.arrange(p2,p3,  ncol=2)
#dev.off()


### 7) PLOT COVARIATES DISTRIBUTION IN BATCHES (Sequencing)
df_seq = lapply(c("mouseID","day","treatment"), function(cov){
data.frame(rbind(t(sapply(1:max(num_sequencing), function(i) (table(c(names(obsFreq[[cov]]),covdesc[ibatches_sequencing[[i]],cov]))-1)/table(covdesc[,cov])))), paste("Batch_seq",1:max(num_sequencing),sep=""))
})
for(i in 1:length(df_seq)){colnames(df_seq[[i]]) = c(names(obsFreq[[i]]),"Batch_seq")}
df_seq.final = lapply(df_seq,melt)
for(i in 1:length(df_seq.final)){colnames(df_seq.final[[i]]) = c("Batch_seq",c("mouseID","day","treatment")[i],"Frequency")}

p1_seq <- ggplot(data=df_seq.final[[1]], aes(fill=mouseID, y=Frequency, x=Batch_seq)) + geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p2_seq <- ggplot(data=df_seq.final[[2]], aes(fill=Batch_seq, y=Frequency, x=day)) + geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p3_seq <- ggplot(data=df_seq.final[[3]], aes(fill=Batch_seq, y=Frequency, x=treatment)) + geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

#pdf(paste0(path,"mouse_batching_sequencing.pdf"))
grid.arrange(p2_seq,p3_seq,  ncol=2)
#dev.off()



### 8) PLOT COVARIATES DISTRIBUTION IN BATCHES (miRNA Sequencing)
df_seq_miRNA = lapply(c("mouseID","day","treatment"), function(cov){
data.frame(rbind(t(sapply(1:max(num_sequencing_miRNA), function(i) (table(c(names(obsFreq[[cov]]),covdesc[ibatches_sequencing_miRNA[[i]],cov]))-1)
/table(covdesc[,cov])))), paste("Batch_seq",1:max(num_sequencing_miRNA),sep=""))
})
for(i in 1:length(df_seq)){colnames(df_seq[[i]]) = c(names(obsFreq[[i]]),"Batch_seq")}
df_seq_miRNA.final = lapply(df_seq_miRNA,melt)
for(i in 1:length(df_seq_miRNA.final)){colnames(df_seq_miRNA.final[[i]]) = c("Batch_seq",c("mouseID","day","treatment")[i],"Frequency")}

p1_seq_miRNA <- ggplot(data=df_seq_miRNA.final[[1]], aes(fill=mouseID, y=Frequency, x=Batch_seq)) + geom_bar(stat="identity", position=position_dodge()) +
theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p2_seq_miRNA <- ggplot(data=df_seq_miRNA.final[[2]], aes(fill=day, y=Frequency, x=Batch_seq)) + geom_bar(stat="identity", position=position_dodge()) + 
theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p3_seq_miRNA <- ggplot(data=df_seq_miRNA.final[[3]], aes(fill=treatment, y=Frequency, x=Batch_seq)) + geom_bar(stat="identity", position=position_dodge())
+ theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

#pdf(paste0(path,"mouse_batching_sequencing_miRNA.pdf"))
grid.arrange(p2_seq_miRNA,p3_seq_miRNA,  ncol=2)
#dev.off()


### 9) WRITE BATCHING IN A CSV

df.final = cbind(data.frame(covdesc[unlist(ibatches),]),batch_library=rep(1:nbBatches,sapply(ibatches,length)))
df.final$batch_sequencing_mRNA=num_sequencing
df.final$batch_sequencing_miRNA=num_sequencing_miRNA



write.table(df.final,sep=";",row.names=FALSE,file=paste0(path,"mouse_batching_meatboproteo_v1.csv"))
