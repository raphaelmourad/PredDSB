# Raphael Mourad
# University Paul Sabatier, France
# raphael.mourad@ibcg.biotoul.fr
# 18/10/2017

# This R script predicts double-strand DNA breaks (DSBs) in NHEK cells using NHEK epigenomics and chromatin data.
# NHEK DSBCapture DSB data and NHEK ChIP-seq and NHEK DNase-seq are used to train a classifier (random forest or lasso logistic regression).
# Then a testing DSB dataset is used to validate DSB predictions in NHEK.
# Note that this script can be used to predict NHEK DSBs using DNA motif and DNA shape.
# Note also that this script can be used to train a model using NHEK data. The NHEK-trained model can then be used for U2OS prediction (see script predictDSBU2OS.R). 


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


# OPTION ------------------------------------------------------------
mode="Epigenome" # To predict using Epigenomic and Chromatin data (ChIP-seq and DNase-seq data)
#mode="EpigenomeForU2OS" # To train a model using Epigenomic and Chromatin data that are available in both NHEK and U2OS cells, i.e.: DNA-seq, CTCF, H3K4me1/3, H3K9me3, H3K27ac, H3K27me3, H3K36me3 and POL2B. 
#mode="Motif" # To predict using DNA motif data only. 
#mode="Motif+Shape" # To predict using DNA motif data and DNA shape. 


# FILES ------------------------------------------------------------

fileAnnot=list.files("data/Epigenome/")
AnnotNames=as.vector(sapply(fileAnnot,function(x){strsplit(x,'_')[[1]][1]}))


# DATA IMPORT --------------------------------------------------------------------

# Import breaks
fileBedBreaksPos="data/DSB/breakome_DSBcap_hg19_20kseq.bed"
dataBreaksPos.GR=sort(readGFBed(fileBedBreaksPos,SeqInfo))
fileBedBreaksNeg="data/DSB/breakome_DSBcap_hg19_20kseq_neg.bed"
dataBreaksNeg.GR=sort(readGFBed(fileBedBreaksNeg,SeqInfo))
dataBreaks.GR=c(dataBreaksPos.GR,dataBreaksNeg.GR)

# Import other data
if(mode=="Epigenome"){
 GenomicFeatureList.GR=list()
 for(i in 1:length(AnnotNames)){
  GenomicFeatureList.GR[[i]] <- sort(unique(readGFBed(paste0("data/Epigenome/",fileAnnot[i]),SeqInfo)))
  print(paste0(AnnotNames[i]," : ",length(GenomicFeatureList.GR[[i]])))
 }
 names(GenomicFeatureList.GR)=AnnotNames
}else if(mode=="EpigenomeForU2OS"){
 GenomicFeatureList.GR=list()
 for(i in 1:length(AnnotNames)){
  GenomicFeatureList.GR[[i]] <- sort(unique(readGFBed(paste0("data/Epigenome/",fileAnnot[i]),SeqInfo)))
  print(paste0(AnnotNames[i]," : ",length(GenomicFeatureList.GR[[i]])))
 }
 names(GenomicFeatureList.GR)=AnnotNames
 GenomicFeatureList.GR=GenomicFeatureList.GR[names(GenomicFeatureList.GR)%in%c("CTCF","DNase","H2az","H3K27ac","H3K27me3","H3K36me3","H3K4me1","H3K4me3","H3K9me3","POL2B")]
 AnnotNames=names(GenomicFeatureList.GR)
}else if(mode=="Motif" | mode=="Motif+Shape"){
 # Load PFMs
 opts <- list()
 opts[["species"]] <- 9606 # human
 opts[["all_versions"]] <- TRUE
 PFMatrixList <- getMatrixSet(JASPAR2016, opts)
}


# MAPPING FEATURES --------------------------------------------------------------------

# Map features to DSB and non-DSB regions 
bin.Mat=c(rep(1,length(dataBreaksPos.GR)),rep(0,length(dataBreaksNeg.GR)))
if(mode=="Epigenome" | mode=="EpigenomeForU2OS"){
 for(i in 1:length(GenomicFeatureList.GR)){
  GRi=GenomicFeatureList.GR[[i]]
  annotPosi=annotateLoci(dataBreaksPos.GR,GRi)
  annotNegi=annotateLoci(dataBreaksNeg.GR,GRi)
  annoti=c(annotPosi,annotNegi)
  annoti[annoti>1]=1
  bin.Mat=cbind(bin.Mat,annoti)
  rm(annoti)
  print(paste0(AnnotNames[i]," annotated"))
 }
 colnames(bin.Mat)=c("Breaks",AnnotNames)
}else if(mode=="Motif" | mode=="Motif+Shape"){
 if(F){# It can takes 20h to annotate DNA motifs to DSB and non-DSB regions!
 matMotifs=findMotifs(dataBreaks.GR,PFMatrixList,BSgenome.Hsapiens.UCSC.hg19.masked)
 save(matMotifs,file="data/Motif/matMotifs.RData")
 }# Preprocessed data are available.
 load(file="data/Motif/matMotifs.RData")
 bin.Mat=cbind(bin.Mat,matMotifs)
 rownames(bin.Mat)=1:nrow(bin.Mat)
 colnames(bin.Mat)[1]="Breaks"

 if(mode=="Motif+Shape"){
  load("data/Shape/matShape.RData")
  bin.Mat=cbind(bin.Mat,matShape)
 }
}


# PREDICTIONS --------------------------------------------------------------------

dir.create(paste0("results/pred",mode))

# Training and testing sets
dataDSB=data.frame(bin.Mat)
rownames(dataDSB)=1:nrow(dataDSB)
idxs=sample(1:nrow(dataDSB),3e4)
dataDSBlearn=dataDSB[sort(idxs),]
dataDSBtest=dataDSB[-idxs,]

# Random Forests
RFall=ranger("Breaks~.",data=dataDSBlearn,importance="permutation")
rocRFall=roc(as.factor(dataDSBtest[,1]),predict(RFall,dataDSBtest)$predictions,ci=T)
aucRF=pROC::auc(rocRFall)

varimp=data.frame(Feature=names(RFall$variable.importance),VariableImportance=RFall$variable.importance)
varimp=varimp[order(varimp[,2],decreasing=T),]
file_varimp=paste0("results/pred",mode,"/varimpRF_",mode,".csv")
write.table(varimp,file=file_varimp,row.names=F,sep='\t',quote=F)

file_rocRF=paste0("results/pred",mode,"/rocRF_",mode,".pdf")
pdf(file_rocRF,4,4)
plot(rocRFall,main=paste0("AUC: ",round(aucRF,4)))
dev.off()

if(mode=="EpigenomeForU2OS"){
 save(RFall,file=paste0("results/pred",mode,"/RF_",mode,"_10vars.RData"))
}


# Lasso logistic regression
CVLasso=cv.glmnet(as(dataDSBlearn[,-1],"Matrix"),dataDSBlearn[,1],family="binomial",parallel=F)
lambda=CVLasso$lambda.min # CVLasso$lambda.min or CVLasso$lambda.1se
CVLassoError=CVLasso$cvm[which(CVLasso$lambda==lambda)]
devLasso=deviance.glmnet(CVLasso$glmnet.fit)[which(CVLasso$lambda==lambda)]
coefLasso=CVLasso$glmnet.fit$beta[,which(CVLasso$lambda==lambda)]
coefLassoMat=data.frame(Variable=names(coefLasso),Coefficient=round(coefLasso,5))
write.table(coefLassoMat,file=paste0("results/pred",mode,"/coefLassoMat_",mode,".csv"),row.names=F,sep='\t',quote=F)

rocLasso=roc(as.factor(dataDSBtest[,1]),predict(CVLasso,as(dataDSBtest[,-1],"Matrix")),ci=T)
aucLasso=pROC::auc(rocLasso)

file_rocLasso=paste0("results/pred",mode,"/rocLasso_",mode,".pdf")
pdf(file_rocLasso,4,4)
plot(rocLasso,main=paste0("AUC: ",round(aucLasso,4)))
dev.off()

if(mode=="EpigenomeForU2OS"){
 save(CVLasso,file=paste0("results/pred",mode,"/Lasso_",mode,"_10vars.RData"))
}




