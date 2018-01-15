# Raphael Mourad
# University Paul Sabatier, France
# raphael.mourad@ibcg.biotoul.fr
# 18/10/2017


# Script to randomly select genomic sequences that matched sizes, GC and repeat contents of DSB sites using R package gkmSVM. These random regions are used as negative observations (regions without double-strand breaks) during model learning. 
# NB: These regions have been already computed and stored in the file "breakome_DSBcap_hg19_20kseq_neg.bed".

setwd("/DSBpred_Rcode")


library(gkmSVM)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg19.masked)
library(Biostrings)




# Draw the regions.
Genome=BSgenome.Hsapiens.UCSC.hg19
fileBedBreaks="data/DSB/breakome_DSBcap_hg19_20kseq.bed"
fileFastaPos="data/DSB/breakome_DSBcap_hg19_20kseq_pos.fa"
fileBedNeg="data/DSB/breakome_DSBcap_hg19_20kseq_neg.bed"
fileFastaNeg="data/DSB/breakome_DSBcap_hg19_20kseq_neg.fa"

genNullSeqs(inputBedFN=fileBedBreaks,nMaxTrials=10,xfold=2,genomeVersion="hg19",
	outputPosFastaFN=fileFastaPos,outputBedFN=fileBedNeg,outputNegFastaFN=fileFastaNeg)


file.remove(c("data/DSB/breakome_DSBcap_hg19_20kseq_pos.fa","data/DSB/breakome_DSBcap_hg19_20kseq_neg.fa"))


