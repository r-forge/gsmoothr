gcContentCalc <- function(locations, window = 500, organism, verbose=FALSE) {	
	return(sequenceCalc(locations, window, organism, DNAString("S"), fixed=FALSE, Nmask=TRUE, verbose=verbose)/window)
}
