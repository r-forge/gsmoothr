processAffyAnnotation <- function(csvFile, skip=19, getRefseq=FALSE, ..., verbose=TRUE) {

  if (verbose)
    cat("Reading file:", csvFile,"\n")
  anno <- read.csv(csvFile,sep=",",skip=skip,header=TRUE,comment.char="",stringsAsFactors=FALSE)


  if (verbose)
    cat("Gathering source of annotation.\n")
  mr<-as.list(anno$mrna_assignment)
  names(mr)<-anno$probeset_id
  
  ms<-sapply(mr,FUN=function(u) {
    s<-strsplit(u,split=" // ")[[1]];
    if(length(s)>1) {
      if(s[2]=="ENSEMBL") {
        gtype <- strsplit(s[3]," ")[[1]][1]
        return( paste(s[2],gtype,sep="-") )
      } else {
        return(s[2])
      }
    } else {
      return(s[1])
    }
  })
  
  if (verbose)
    cat("Parsing gene symbol.\n")
  g<-as.list(anno$gene_assignment)
  names(g)<-anno$probeset_id
  h<-sapply(g,FUN=function(u) { s<-strsplit(u,split=" // ", fixed=TRUE)[[1]]; ifelse(length(s)>1,s[2],s[1])})
  
  anno <- anno[,c("probeset_id","seqname","strand","start","stop","total_probes","category")]
  if( getRefseq ) {
    r<-sapply(g,FUN=function(u) strsplit(u,split=" // ", fixed=TRUE)[[1]][1])
		data.frame(anno,id=r,symbol=h,type=ms)
  } else {
		data.frame(anno,symbol=h,type=ms)
  }
}