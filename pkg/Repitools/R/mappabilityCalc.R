mappabilityCalc <- function(locations, window = 500, organism, verbose=FALSE) {	
	return(sequenceCalc(locations, window, organism, DNAString("N"), verbose=verbose)/window)
}
