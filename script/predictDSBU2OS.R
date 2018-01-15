# Raphael Mourad
# University Paul Sabatier, France
# raphael.mourad@ibcg.biotoul.fr
# 15/01/2018

# This R script predicts double-strand DNA breaks (DSBs) in U2OS cells using U2OS epigenomics and chromatin data with an NHEK-trained model.
# First, NHEK DSBCapture DSB data and NHEK ChIP-seq and NHEK DNase-seq were used to train a classifier (random forest or lasso logistic regression).
# Second, U2OS epigenomics and chromatin data is used to predict U2OS DSBs with an NHEK-trained model. 
# Third, U2OS DSB predictions are compared to U2OS DSBCapture DSBs.



setwd("DSBpred_Rcode")


library(pROC)
library(glmnet)
library(ranger)
library(Matrix)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
library(JASPAR2016)
library(TFBSTools)
source("script/miscFunctions.R")


Genome=BSgenome.Hsapiens.UCSC.hg19
model="human"
SeqinfoGenome=seqinfo(BSgenome.Hsapiens.UCSC.hg19)
Chr.V=c(paste0("chr",1:22),"chrX")
SeqInfo=SeqinfoGenome[Chr.V]


# DATA IMPORT ------------------------------------------------------------
# U2OS epigenome and chromatin data
fileAnnot=list.files("data/Epigenome_U2OS/")
AnnotNames=as.vector(sapply(fileAnnot,function(x){strsplit(x,'_')[[1]][1]}))

GenomicFeatureList.GR=list()
for(i in 1:length(AnnotNames)){
 GenomicFeatureList.GR[[i]] <- sort(unique(readGFBed(paste0("data/Epigenome_U2OS/",fileAnnot[i]),SeqInfo)))
 print(paste0(AnnotNames[i]," : ",length(GenomicFeatureList.GR[[i]])))
}
names(GenomicFeatureList.GR)=AnnotNames

# U2OS DSBCapture peaks
DSB_endo_pos.GR=sort(readGFBed("data/DSB_U2OS/breakome_DSBcap_hg19_U2OS_endo.bed",SeqInfo))
DSB_endo_neg.GR=sort(readGFBed("data/DSB_U2OS/breakome_DSBcap_hg19_U2OS_endo_neg.bed",SeqInfo))
start(DSB_endo_pos.GR)=start(DSB_endo_pos.GR)-100
end(DSB_endo_pos.GR)=end(DSB_endo_pos.GR)+100
start(DSB_endo_neg.GR)=start(DSB_endo_neg.GR)-100
end(DSB_endo_neg.GR)=end(DSB_endo_neg.GR)+100
DSB_endo_pos.GR=trim(DSB_endo_pos.GR)
DSB_endo_neg.GR=trim(DSB_endo_neg.GR)
DSB_endo_neg.GR=sort(DSB_endo_neg.GR[sample(1:length(DSB_endo_neg.GR),length(DSB_endo_pos.GR))])


# NHEK-TRAINED MODEL ------------------------------------------------------------
# This random forest was previously computed using NHEK data with the script predictDSB.R and the mode="EpigenomeForU2OS".
load(file="results/predEpigenomeForU2OS/RF_EpigenomeForU2OS_10vars.RData")


# MAPPING FEATURES -------------------------------------------------------------- 
bin.Mat=c(rep(1,length(DSB_endo_pos.GR)),rep(0,length(DSB_endo_neg.GR)))
for(i in 1:length(GenomicFeatureList.GR)){
 GRi=GenomicFeatureList.GR[[i]]
 annotPosi=annotateLoci(DSB_endo_pos.GR,GRi)
 annotNegi=annotateLoci(DSB_endo_neg.GR,GRi)
 annoti=c(annotPosi,annotNegi)
 annoti[annoti>1]=1
 bin.Mat=cbind(bin.Mat,annoti)
 rm(annoti)
 print(paste0(AnnotNames[i]," annotated"))
}
colnames(bin.Mat)=c("Breaks",AnnotNames)


# PREDICTIONS IN U2OS USING NHEK-TRAINED MODEL ----------------------------------------
dir.create("results/predEpigenomeU2OS")

pred=predict(RFall,bin.Mat)$predictions
rocRF=roc(as.factor(bin.Mat[,1]),pred,ci=T)
aucRF=pROC::auc(rocRF)

fileRF="results/predEpigenomeU2OS/rocRF_Epigenome_U2OS_trainedNHEK.pdf"
pdf(fileRF,4,4)
plot(rocRF,main=paste0("AUC: ",round(aucRF,3)))
dev.off()


# PREDICTIONS IN U2OS USING U2OS-TRAINED MODEL ----------------------------------------
RF_U2OS=ranger("Breaks~.",data=as.data.frame(bin.Mat),importance="permutation")
rocRF_U2OS=roc(as.factor(bin.Mat[,1]),predict(RF_U2OS,bin.Mat)$predictions,ci=T)
aucRF_U2OS=pROC::auc(rocRF_U2OS)

fileRF_U2OS="results/predEpigenomeU2OS/rocRF_Epigenome_U2OS_trainedU2OS.pdf"
pdf(fileRF_U2OS,4,4)
plot(rocRF_U2OS,main=paste0("AUC: ",round(aucRF_U2OS,3)))
dev.off()

fileVI_U2OS="results/predEpigenomeU2OS/barplot_VI_predModelU2OS.pdf"
pdf(fileVI_U2OS,6,6)
barplot(sort(RF_U2OS$variable.importance,decreasing=F),horiz=T,las=2)
dev.off()














