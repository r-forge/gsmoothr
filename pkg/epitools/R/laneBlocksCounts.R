laneBlocksCounts <- function(lanesPath, TSSDataTable, upstream = 1000, downstream = 1000, removeDuplications=FALSE, doQuality=FALSE, aligner="Bowtie")
{
	require(epitools)
	require(IRanges)
	require(ShortRead)
	
	if(!all(c("chr", "name", "start", "end", "strand")  %in% colnames(TSSDataTable)))
		stop("Incorrect column headings for TSSDataTable. Check documentation for details.")

	fileNames <- dir(lanesPath, pattern="fastq")
	fileNames <- gsub(".fastq", "", fileNames)
	
	startPositions <- as.numeric(apply(TSSDataTable, 1, function(x) ifelse(x["strand"]=="+", as.numeric(x["start"])-upstream, as.numeric(x["end"])-downstream)))
	endPositions <- as.numeric(apply(TSSDataTable, 1, function(x) ifelse(x["strand"]=="+", as.numeric(x["start"])+downstream, as.numeric(x["end"])+upstream)))
	
	TSSranges <- mapply(IRanges, start=split(startPositions, TSSDataTable$chr), end=split(endPositions, TSSDataTable$chr), names=split(as.character(TSSDataTable$name), TSSDataTable$chr))
	namesOnEachChromosome <- split(as.character(TSSDataTable$name), TSSDataTable$chr)
	
	TSSL <- do.call(RangesList, TSSranges)
	readCounts <- matrix(nrow=nrow(TSSDataTable), ncol = length(fileNames))
	rownames(readCounts) <- TSSDataTable$name
	colnames(readCounts) <- fileNames
	if(removeDuplications == TRUE)
	{
		readCountsNoDuplications <- matrix(nrow = nrow(TSSDataTable), ncol = length(fileNames))
		rownames(readCountsNoDuplications) <- TSSDataTable$name
		colnames(readCountsNoDuplications) <- fileNames
	}
		
	if(aligner == "Bowtie" || "MAQ")
	{
		extension <- ".map"
	} else { # is SolexaExport
		extension <- ".txt"
	}
 

	for(index in 1:length(fileNames))
	{
		# If directory corresponding to file doesn't exist.
		if(!file.exists(paste(lanesPath, fileNames[index], sep="/")))
		{
			stop("Can't find your alignment files.")
		}
		
		writeAccess <- file.access(paste(lanesPath,"/",fileNames[index]), 2)
		cat("Warning : Unable to write to directory ", paste(lanesPath, fileNames[index], sep = "/"), ". Can't save alignment in .rda file or write out quality reports.\n")
		
		if(doQuality == TRUE && writeAccess == 0)
		{
			qualityAssessment <- qa(lanesPath, paste(fileNames[index],"\\.fastq", sep=""), type="fastq")
			report(qualityAssessment, dest = paste(lanesPath, fileNames[index], "fastqReport", sep="/"))
		}
		
		#no duplicate filter
		filename <- paste(lanesPath,"/",fileNames[index],"/",fileNames[index],"dups",".rda",sep="")
		if (file.exists(filename)) 
			load(filename)
		else
		{
		   readAligned.dups <- readAligned(paste(lanesPath,"/",fileNames[index],sep=""), paste(fileNames[index], extension, sep=""), type=aligner)
		   if(writeAccess == 0)
			save(readAligned.dups, file=paste(lanesPath,"/",fileNames[index],"/",fileNames[index],"dups",".rda",sep=""))
		}
		
		if(doQuality == TRUE && writeAccess == 0)
		{
			qualityAssessment <- qa(paste(lanesPath,"/",fileNames[index],sep=""), paste("\\", extension, sep=""), type=aligner)
			report(qualityAssessment, dest = paste(lanesPath,fileNames[index], "alignmentReport", sep="/"))
		}

		chr <- chromosome(readAligned.dups)
		readsRanges <- mapply(IRanges, start=split(position(readAligned.dups), chr), width=split(width(readAligned.dups), chr))
		alignedL <- do.call(RangesList, readsRanges)

		countTables <- lapply(overlap(alignedL, TSSL), as.table)
		for(chromosomeIndex in 1 : length(namesOnEachChromosome))
		{
			readCounts[namesOnEachChromosome[[chromosomeIndex]],index] <- countTables[[chromosomeIndex]]
		}
		
		rm(readAligned.dups, readsRanges, alignedL)
		gc()
		
		#duplicate filter
		if(removeDuplications == TRUE)
		{
			filename <- paste(fileNames[index],"/",fileNames[index],"nodups",".rda",sep="")
			if (file.exists(filename))
				load(filename)
			else
			{
			   readAligned.nodups <- readAligned(paste(lanesPath,"/",fileNames[index],sep=""), paste(fileNames[index], extension, sep=""), type=aligner, filter=uniqueFilter(withSread = FALSE))
			   if(writeAccess == 0)
				save(readAligned.nodups, file=paste(lanesPath,"/",fileNames[index],"/",fileNames[index],"nodups",".rda",sep=""))
			}

			chr <- chromosome(readAligned.nodups)
			readsRanges <- mapply(IRanges, start=split(position(readAligned.nodups), chr), width=split(width(readAligned.nodups), chr))
			alignedL <- do.call(RangesList, readsRanges)
		
			countTables <- lapply(overlap(alignedL, TSSL), as.table)
			for(chromosomeIndex in 1 : length(namesOnEachChromosome))
			{
				readCountsNoDuplications[namesOnEachChromosome[[chromosomeIndex]],index] <- countTables[[chromosomeIndex]]	
			}
			
			rm(readAligned.nodups, readsRanges, alignedL)
			gc()
		}
	}
	
	if(removeDuplications==TRUE)
		return(list(readCounts=readCounts, readCountsNoDuplications=readCountsNoDuplications))
	else
		return(list(readCounts=readCounts, readCountsNoDuplications=NULL))
}