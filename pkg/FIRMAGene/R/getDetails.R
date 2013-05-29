setMethodS3("getDetails", "ProbeLevelModel",
   function(plm,probesets,id="7952953",mart=NULL,verticalBars=FALSE,geneSymbolId="external_gene_id",
            transcriptClusterId="Transcript.Cluster.ID",colours="blue",lwd=1,o=NULL, verbose=FALSE) {
  library("GenomeGraphs")
  library("biomaRt")
  cs <- getDataSet(plm)
  rs <- getResidualSet(plm)
  cdf <- getCdf(cs)
  if (is.integer(id)) {
    unit <- id
	  id <- getUnitNames(cdf)[unit]
  } else {
    unit <- which( getUnitNames(cdf)== id)
  }
  ind <- unlist(getCellIndices(getCdf(cs),units=unit),use.names=FALSE)
  d <- log2(extractMatrix(cs,cells=ind,verbose=verbose))
  r <- log2(extractMatrix(rs,cells=ind,verbose=verbose))
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
  m <- which(probesets[[transcriptClusterId]]==id)  # assume this matches with the order in the CDF!!!
  pm <- probesets[m,]
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
  if (verbose)
    print(pm[1,])
  ch <- gsub("chr","",pm[1,"seqname"])
  gr<-new("GeneRegion", chromosome = ch,
                   start = as.numeric(min(pm[1,"start"])), end = as.numeric(pm[1,"stop"]), strand = as.character(pm[1,"strand"]), biomart = mart)
  gid <- names(sort(table(gr@ens[,"ensembl_gene_id"]),decreasing=TRUE)[1])
  if (verbose)
    print(gid)
  if( nrow(gr@ens)==0 )
    gr <- NULL
	if (verbose)
    print(gr)
 if(is.null(gid))
      b <- ""
 else
      b <- getBM(geneSymbolId, filters = "ensembl_gene_id", values = gid, mart = mart)
  ti <- new("Title", title = paste(getUnitNames(cdf)[unit],gid,b[1],sep=" -- "), dp = DisplayPars(color = "darkred"))
  if( !is.null(gid) )
    tr <- new("Transcript",id=gid,biomart=mart,dp=DisplayPars(plotId=TRUE))
  else
    tr <- NULL
  returnList <- list(ti,ea1,ga,ea2,gr,tr)
  keep <- !sapply(returnList,is.null)
  return(returnList[keep])  
}) # getDetails
