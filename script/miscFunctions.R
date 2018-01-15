# Raphael Mourad
# University Paul Sabatier, France
# raphael.mourad@ibcg.biotoul.fr
# 18/10/2017


# Miscellaneous functions. 


# Read bed file containing ChIP-seq peak data
readGFBed<- function(GFBedFile,seqInfoChr){

 if(class(seqInfoChr)!="Seqinfo"){print("seqInfoChr is not a seqinfo object!");return(0)}
 if(length(grep(pattern='.bed',x=GFBedFile))==0){print("GFBedFile is not a bed file!");return(0)}

 dataGF=read.table(GFBedFile,sep='\t',header=F)[,c(1:3)]
 dataGF2=dataGF[dataGF[,1]%in%seqnames(seqInfoChr),]
 chr=as.character(dataGF2[,1])
 posLeft=as.numeric(dataGF2[,2])
 posRight=as.numeric(dataGF2[,3])

 GenomicFeature.GR=NULL
 for(i in 1:length(seqnames(seqInfoChr))){
  chri=seqnames(seqInfoChr)[i]
  if(sum(chr==chri)){
   GenomicFeature.GRi <- GRanges(seqnames=chri,IRanges(start=posLeft[chr==chri],end=posRight[chr==chri]))
   if(is.null(GenomicFeature.GR)){
    GenomicFeature.GR=GenomicFeature.GRi
   }else{
    GenomicFeature.GR=c(GenomicFeature.GR,GenomicFeature.GRi)
   }
  }
 }

 seqlevels(GenomicFeature.GR) = seqlevels(seqInfoChr)
 seqinfo(GenomicFeature.GR)=seqInfoChr

 return(GenomicFeature.GR)
}


# Map a genomic feature to loci. 
annotateLoci <- function(loci.GR, GenomicFeature.GR)
{
 if(!is(loci.GR, "GenomicRanges")){stop("'loci.GR' must be a GenomicRanges object")}
 if(!is(GenomicFeature.GR, "GenomicRanges")){stop("'GenomicFeature.GR' must be a GenomicRanges object")}
 if(any(is.na(seqlengths(GenomicFeature.GR)))){stop("'seqlengths(GenomicFeature.GR)' contains NAs")}

 if(is.null(GenomicFeature.GR$score)){
  cvg <- coverage(GenomicFeature.GR)
 }else{
  cvg <- coverage(GenomicFeature.GR, weight="score")
 }

 chr=seqnames(seqinfo(loci.GR))
 loci.list=list()
 for(i in 1:length(chr)){
  loci.list[[chr[i]]]=loci.GR[seqnames(loci.GR)==chr[i]]
 }
 loci.GRL=GenomicRangesList(loci.list)

 if(length(chr)==1){
  averageBin=viewMeans(Views(cvg[[chr]],ranges(loci.GRL[[chr]])))
 }else{
  views_list <- RleViewsList(lapply(names(cvg),function(seqname)Views(cvg[[seqname]], ranges(loci.GRL[[seqname]]))))
  averageBin=NULL
  for(i in 1:length(views_list)){
   averageBin=c(averageBin,viewMeans(views_list)[[i]])
  }
 }

 return(averageBin)
}


# Function to find protein binding sites predicted from DNA motifs in a sequence
findMotifs<-function(seq.GR,PFML,genome){

 seq=as.list(getSeq(genome, seqnames(seq.GR), start=start(seq.GR), end=end(seq.GR), as.character=TRUE))
 matMotifCount=NULL
 name_id=NULL

 for(i in 1:length(PFML)){ 

  PFMi=PFML[[i]]
  PWMi=as.matrix(toPWM(PFMi))
  PWMci=reverseComplement(PWMi)
  idi=ID(PFMi)
  namei=name(PFMi)
  name_id=c(name_id,paste0(namei,"_",idi))
  
  posMotifPos=unlist(lapply(lapply(seq,function(x){matchPWM(PWMi,x)}),length))
  posMotifNeg=unlist(lapply(lapply(seq,function(x){matchPWM(PWMci,x)}),length))
  posMotifAll=c(posMotifPos,posMotifNeg)
  
  matMotifCount=cbind(matMotifCount,posMotifAll)
  
  print(paste0(namei," annotated"))
 }
 colnames(matMotifCount)=name_id

 return(matMotifCount)
}

