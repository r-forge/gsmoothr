annotationBlocksLookup <- function(probes, annotation, probeIndex=NULL) {
#probes = dataframe of $chr, $position and $strand ("+" or "-")
#annotation = dataframe of $chr, $start, $end, $strand ("+" or "-") and $name or rownames = annotation name

	processChunk <- function(probePositions, annotation) {
		#initialise return variables
		numAnnot = nrow(annotation)
		annotProbes = list(indexes=vector(mode='list', length=numAnnot), offsets=vector(mode='list', length=numAnnot))
		if (is.null(probePositions)||(length(probePositions)==0)) return(annotProbes) #if no probes in this chunk return empty annotations	
		chromosomeSize = max(probePositions)

		#allocate vector of chromosome size
		chromosomeLookup <- rep.int(NA, chromosomeSize) 

		#set positions of probes
		chromosomeLookup[probePositions] = names(probePositions)

		for (i in 1:numAnnot) {
			tempProbes = chromosomeLookup[annotation$start[i]:annotation$end[i]]
			annotProbes$indexes[[i]] = as.integer(na.omit(tempProbes))
			#adjust offsets if at start/end of a chromosome
			annotProbes$offsets[[i]] = which(!is.na(tempProbes))
		}
		return(annotProbes)
	}

	processChromosome <- function(probePositions, annotation) {
		numAnnot = nrow(annotation)
		annotProbes = list(indexes=vector(mode='list', length=numAnnot), offsets=vector(mode='list', length=numAnnot))
		if (length(probePositions)==0) return(annotProbes) #if no probes on this chromosome return empty annotations	
		chromosomeSize = max(probePositions)
    annotChunks = split(1:nrow(annotation), trunc(annotation$start/10000000)) #split into 10MB chunks of annotations
    for (i in annotChunks) {
      chunkRange = range(c(annotation$start[i], annotation$end[i]))
      chunkPositions = probePositions[(probePositions>=chunkRange[1])&(probePositions<=chunkRange[2])]-chunkRange[1]+1
      chunkAnnot = annotation[i,]
      chunkAnnot$start = chunkAnnot$start-chunkRange[1]+1
      chunkAnnot$end = chunkAnnot$end-chunkRange[1]+1
		  tempAnnot = processChunk(chunkPositions, chunkAnnot)
		  annotProbes$indexes[i] = tempAnnot$indexes
		  annotProbes$offsets[i] = tempAnnot$offsets
	  }
		return(annotProbes)
  }


	if (is.null(annotation$strand)) { #dont bother with strandedness anywhere
		probesStrandChr <- probes$chr
		annotationStrandChr <- annotation$chr
	} else { #create separate plus & minus chromosomes
		if (is.null(probes$strand)) stop("error: if 'annotation' contains strand information, 'probes' must contain strand information as well")
		probesStrandChr <- paste(probes$chr, probes$strand, sep="")
		annotationStrandChr <- paste(annotation$chr, annotation$strand, sep="")
	}
	#split by strand AND chromosome simultaneously
	annotChr = split(1:nrow(annotation), annotationStrandChr)
	annot = list(indexes=vector(mode='list', length=nrow(annotation)), offsets=vector(mode='list', length=nrow(annotation)))
	for (i in annotChr) {
		thisChr = annotationStrandChr[i[1]]
		cat("Processing",thisChr,"\n")

		#Grab the subset of probes on that chromosome
		tempIndex = which(probesStrandChr==thisChr)
		tempProbes = probes$position[tempIndex]
		#Use probeIndex supplied, or assume probes are in order
		if (is.null(probeIndex)) names(tempProbes) <- tempIndex else names(tempProbes) <- probeIndex[tempIndex]

		#Process the chromosome
		tempAnnot = processChromosome(tempProbes, annotation[i,])
		annot$indexes[i] = tempAnnot$indexes
		annot$offsets[i] = tempAnnot$offsets
	}
	if (!is.null(rownames(annotation))) {
		names(annot$indexes) <- annotation$name
		names(annot$offsets) <- annotation$name
	}
	return(annot)
	#returns $indexes = a list for each annotation entry with the indexes of the probes within the block
	#	 $offsets = a list for each annotation entry with the offsets from the beginning of the block
	
}

annotationLookup <- function(probes, annotation, bpUp, bpDown, probeIndex=NULL) {
#probes = dataframe of $chr and $position
#annotation = dataframe of $chr, $position, $strand ("+" or "-") and $name or rownames = annotation name
#if annotation has no strand, assume are + strand
	if (is.null(annotation$strand)) annotation$strand <- "+"
	annotationTemp <- data.frame(chr=annotation$chr, 
                                     start=annotation$position+ifelse(annotation$strand=="+",-bpUp, +bpUp),
                                     end=annotation$position+ifelse(annotation$strand=="+",+bpDown, -bpDown),
                                     name=rownames(annotation), stringsAsFactors=F)
	annot <- annotationBlocksLookup(probes, annotationTemp, probeIndex)

	#adjust offset by bpUp & bpDown
	adjustOffset <- function(offsets, strand, bpUp, bpDown) {
		if (strand=="+") return(offsets-bpUp) else return(-1*offsets+bpDown)
	}
#	annot$offsets = mapply(adjustOffset, annot$offsets, annotation$strand, MoreArgs=list(bpUp=bpUp, bpDown=bpDown))
	annot$offsets = lapply(annot$offsets, function(x,bpUp) {return(x-bpUp-1)}, bpUp)
	if (!is.null(rownames(annotation))) {
		names(annot$indexes) <- annotation$name
		names(annot$offsets) <- annotation$name
	}
	return(annot)
}

