setupExamples <- function(doCopy=FALSE) {
	dataDirectory <- system.file("exampleData", package="RepitoolsExamples")
	if  (file.access(dataDirectory, 2) == 0) {
		setwd(dataDirectory)
		cat("Successfully setwd to example data directory\n")
		invisible(NULL)
	} else {
		cat("Example data directory is not writable\n")
		if (doCopy) {
			cat("Copying example data directory to current working directory\n")	
			require(R.utils)
			file.copy(dataDirectory, getwd(), recursive = TRUE)
			cat("Successfully setwd to local copy of example data directory\n")
			invisible(NULL)
		} else cat("If you wish to make a local copy, run setupExamples(doCopy=TRUE)\n")
	}
}
