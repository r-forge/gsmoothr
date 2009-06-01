getDetails<-function(plm,probesets,id="7952953",mart=NULL,verbose=FALSE,biomart="ensembl",verticalBars=FALSE,
                     dataset="hsapiens_gene_ensembl",colours="blue",lwd=1,o=NULL,psId="Transcript.Cluster.ID") {
  require(GenomeGraphs)
  require(biomaRt)
  if(is.null(mart))
    mart=useMart(biomart=biomart, dataset=dataset)
  #cs <- extract(getDataSet(plm),1:3)
  #rs <- extract(getResidualSet(plm),1:3)
  cs <- getDataSet(plm)
  rs <- getResidualSet(plm)
  cdf<-getCdf(cs)
  if (is.integer(id)) {
    unit <- id
	id <- getUnitNames(cdf)[unit]
  } else {
    unit <- which( getUnitNames(cdf)== id)
  }
  ind<-unlist(getCellIndices(getCdf(cs),units=unit),use.names=FALSE)
  d<-log2(extractMatrix(cs,cells=ind,verbose=verbose))
  r<-log2(extractMatrix(rs,cells=ind,verbose=verbose))
  if(is.null(o))
    o <- seq_len(ncol(d))
  if(length(colours)==1)
    colours <- rep(colours,ncol(d))
  if(length(lwd)==1)
    lwd <- rep(lwd,ncol(d))
  lwd <- lwd[o]
  colours <- colours[o]
  d <- d[,o]
  r <- r[,o]
  #nm<-getUnitGroupNamesFromUgcMap(cdf,getUnitGroupCellMap(cdf,units=unit))
  #ugnm<-unique(nm$groupName)
  m<-which(probesets[[psId]]==id)  # assume this matches with the order in the CDF!!!
  #nProbe<-probesets[match(ugnm,probesets$probeset_id),"probe_count"]
  #names(nProbe)<-ugnm  # can re-order these based on info in 'probesets'!!!
  pm<-probesets[m,]
  if( any( colnames(pm) %in% "probe_count" ))
    nProbe <- pm[,"probe_count"]
  else
    nProbe <- rep(1,length(m))
  ea1<-new("ExonArray", intensity = d, probeStart = as.numeric(pm[,"start"]), probeEnd=as.numeric(pm[,"stop"]),
                   probeId = rep("",length(m)), nProbes = nProbe, dp = DisplayPars(color = colours, lwd=lwd,
                   mapColor = "dodgerblue2",plotMap=FALSE, probeSetLwd=as.numeric(verticalBars)), displayProbesets = TRUE)
  ea2<-new("ExonArray", intensity = r, probeStart = as.numeric(pm[,"start"]), probeEnd=as.numeric(pm[,"stop"]),
                   probeId = rep("",length(m)), nProbes = nProbe, dp = DisplayPars(color = colours, lwd=lwd,
                   mapColor = "dodgerblue2", probeSetLwd=as.numeric(verticalBars)), displayProbesets = FALSE)
  ga <- new("GenomeAxis", add53 = TRUE)
  #print(pm)
  print(pm[1,])
  ch <- gsub("chr","",pm[1,"seqname"])
  if(nchar(ch) > 2) {
    ti<-new("Title", title = paste(getUnitNames(cdf)[unit]), dp = DisplayPars(color = "darkred"))
	return(list(ti,ea1,ga,ea2))
  }
  gr<-new("GeneRegion", chromosome = ch,
                   start = as.numeric(min(pm[1,"start"])), end = as.numeric(pm[1,"stop"]), strand = as.character(pm[1,"strand"]), biomart = mart)
				   print(gr)
  #gid<-gr@ens[1,"gene_stable_id"]
  gid<-names(sort(table(gr@ens[,"ensembl_gene_id"]),decreasing=TRUE)[1])
  print(gid)
  b<-getBM("hgnc_symbol",filters="ensembl_gene_id",values=gid,mart=mart)
  ti<-new("Title", title = paste(getUnitNames(cdf)[unit],gid,b[1,1],sep=" -- "), dp = DisplayPars(color = "darkred"))
  if( is.null(gid) )
    return(list(ti,gr))
  tr<-new("Transcript",id=gid,biomart=mart,dp=DisplayPars(plotId=TRUE))
  list(ti,ea1,ga,ea2,gr,tr)
}
