# Raphael Mourad
# University Paul Sabatier, France
# raphael.mourad@ibcg.biotoul.fr
# 18/10/2017


# Script to compute DNA shape features of DSB and non-DSB regions.
# NB: shape data have already been precomputed and stored in the folder data/Shape.


setwd("DSBpred_Rcode")


library(BSgenome.Hsapiens.UCSC.hg19)
library(DNAshapeR)
source("script/miscFunctions.R")



SeqinfoGenome=seqinfo(BSgenome.Hsapiens.UCSC.hg19)
Chr.V=c(paste0("chr",1:22),"chrX")
SeqInfo=SeqinfoGenome[Chr.V]


# Import breaks
fileBedBreaksPos="data/DSB/breakome_DSBcap_hg19_20kseq.bed"
dataBreaksPos.GR=sort(readGFBed(fileBedBreaksPos,SeqInfo))
fileBedBreaksNeg="data/DSB/breakome_DSBcap_hg19_20kseq_neg.bed"
dataBreaksNeg.GR=sort(readGFBed(fileBedBreaksNeg,SeqInfo))

dataBreaksPos.GR=resize(dataBreaksPos.GR,1,fix="center")
dataBreaksPos.GR=resize(dataBreaksPos.GR,100,fix="center")
dataBreaksNeg.GR=resize(dataBreaksNeg.GR,1,fix="center")
dataBreaksNeg.GR=resize(dataBreaksNeg.GR,100,fix="center")


# Compute DNA shape
dir.create("data/Shape")
getFasta(dataBreaksPos.GR,BSgenome=BSgenome.Hsapiens.UCSC.hg19,width=100,filename="data/Shape/breakome_DSBcap_hg19_shape.fa")
getFasta(dataBreaksNeg.GR,BSgenome=BSgenome.Hsapiens.UCSC.hg19,width=100,filename="data/Shape/breakome_DSBcap_hg19_neg_shape.fa")

getShape(filename="data/Shape/breakome_DSBcap_hg19_shape.fa")
getShape(filename="data/Shape/breakome_DSBcap_hg19_neg_shape.fa")

shapeName=c("HelT","MGW","ProT","Roll")
dnaShapePos=NULL
dnaShapeNeg=NULL
for(i in 1:4){
 dnaShapePosi=readShape(filename=paste0("data/Shape/breakome_DSBcap_hg19_shape.fa.",shapeName[i]))
 dnaShapeNegi=readShape(filename=paste0("data/Shape/breakome_DSBcap_hg19_neg_shape.fa.",shapeName[i]))

 meanPosi=apply(dnaShapePosi,1,median,na.rm=T)
 meanNegi=apply(dnaShapeNegi,1,median,na.rm=T)

 qPosi=t(apply(dnaShapePosi,1,quantile,probs=seq(0,1,0.1),na.rm=T))
 qNegi=t(apply(dnaShapeNegi,1,quantile,probs=seq(0,1,0.1),na.rm=T))

 varPosi=apply(dnaShapePosi,1,var,na.rm=T)
 varNegi=apply(dnaShapeNegi,1,var,na.rm=T)

 dnaShapePosi=cbind(qPosi,varPosi)
 dnaShapeNegi=cbind(qNegi,varNegi)
 colnames(dnaShapePosi)=c(paste0(shapeName[i],seq(0,1,0.1)),paste0(shapeName[i],"var"))
 colnames(dnaShapeNegi)=c(paste0(shapeName[i],seq(0,1,0.1)),paste0(shapeName[i],"var"))

 dnaShapePos=cbind(dnaShapePos,dnaShapePosi)
 dnaShapeNeg=cbind(dnaShapeNeg,dnaShapeNegi)
}

matShape=as.data.frame(rbind(dnaShapePos,dnaShapeNeg))
save(matShape,file="data/Shape/matShape.RData")

file.remove(c("data/Shape/breakome_DSBcap_hg19_neg_shape.fa","data/Shape/breakome_DSBcap_hg19_neg_shape.fa.HelT",
              "data/Shape/breakome_DSBcap_hg19_neg_shape.fa.MGW","data/Shape/breakome_DSBcap_hg19_neg_shape.fa.ProT",
              "data/Shape/breakome_DSBcap_hg19_neg_shape.fa.Roll","data/Shape/breakome_DSBcap_hg19_shape.fa",
              "data/Shape/breakome_DSBcap_hg19_shape.fa.HelT","data/Shape/breakome_DSBcap_hg19_shape.fa.MGW",
              "data/Shape/breakome_DSBcap_hg19_shape.fa.ProT","data/Shape/breakome_DSBcap_hg19_shape.fa.Roll"))


















