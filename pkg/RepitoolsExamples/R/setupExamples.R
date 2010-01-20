setupExamples <- function(doCopy=FALSE) {
	dataDirectory <- system.file("exampleData", package="RepitoolsExamples")
	if  (dir.create(paste(dataDirectory,"test", sep="/"), showWarnings=FALSE)) {
		x = file.remove(paste(dataDirectory,"test", sep="/"))
		setwd(dataDirectory)
		cat("Successfully setwd to example data directory\n")
		invisible(NULL)
	} else {
		cat("Example data directory is not writable\n")
		if (doCopy) {
			cat("Copying example data directory to current working directory\n")
			require(R.utils)
			copyDirectory(dataDirectory)
			cat("Successfully setwd to local copy of example data directory\n")
			invisible(NULL)
		} else cat("If you wish to make a local copy, run setupExamples(doCopy=TRUE)\n")
	}
}
