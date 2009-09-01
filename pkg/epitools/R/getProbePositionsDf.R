getProbePositionsDf <- function(cdf, ..., verbose=-20) {

    require(aroma.affymetrix)

	chrs <- paste( "chr", c(1:22,"X","Y","M"), sep="" )
    names(chrs) <- 1:25
	
    ind <- getCellIndices(cdf,...,useNames=FALSE,unlist=TRUE,verbose=verbose)
	acp <- AromaCellPositionFile$byChipType(getChipType(cdf))
    ch <- acp[ind,1,drop=TRUE]
    sp <- acp[ind,2,drop=TRUE]

    data.frame(chr=chrs[as.character(ch)],position=sp,index=ind,stringsAsFactors=FALSE)
}

ndfFilePath <- dir(designDir, "090410_HG18_JS_ChIP.ndf", full = TRUE)
ndfAll <- processNDF(ndfFilePath)
ndfAll <- ndfAll[grepl("^chr|^SPIKE|^RANDOM", ndfAll$chr), ]
columnsOnChip <- max(RGTobyAzaCopy$genes$X)
matchIndices <- match(ndfAll$index, RGTobyAzaCopy$genes$Y * columnsOnChip + RGTobyAzaCopy$genes$X)
RGTobyAzaCopy <- RGTobyAzaCopy[matchIndices, ]
colnames(RGTobyAzaCopy$G) <- RGTobyAzaCopy$targets$Cy3
rownames(RGTobyAzaCopy$G) <- RGTobyAzaCopy$genes$PROBE_ID
RGTobyAzaCopy$G <- normalizeQuantiles(RGTobyAzaCopy$G)


data.frame(chr=ch,position=sp,stringsAsFactors=FALSE)

arrayColourSpotNames <- RG$targets$Cy3
designMatrix <- matrix(0, nrow = length(arrayColourSpotNames), ncol = 1,
dimnames = list(arrayColourSpotNames, "LNCAP Aza - PREC Aza"))
LNCapIP <- grep("LNCaP", arrayColourSpotNames)
designMatrix[LNCapIP, 1] <- 1/length(LNCapIP)
PrECIP <- grep("PrEC", arrayColourSpotNames)
designMatrix[PrECIP, 1] <- -1/length(PrECIP)