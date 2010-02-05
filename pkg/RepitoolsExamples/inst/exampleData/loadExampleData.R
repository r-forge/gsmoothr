require(aroma.affymetrix)

designMatrix <- matrix(c(-1, 0.5, 0.5, 1, -0.5, -0.5), ncol=1)
colnames(designMatrix) <- c("H2A LNCaP - PrEC")
cdfFile <- AffymetrixCdfFile$byChipType("Hs_PromPR_v02")
celSetTables <- AffymetrixCelSet$byName("FaH2A", cdf = cdfFile)
MATNormalise <- MatNormalization(celSetTables, numChunks = 15)
celSetMATNormalised <- process(MATNormalise, verbose=-10)
cdfTableUniquePositions <- getUniqueCdf(cdfFile)
celSetUniquePositions <- convertToUnique(celSetTables)
celSetMATNormalisedUniquePositions <- convertToUnique(celSetMATNormalised)
experimentNames <- toupper(getNames(celSetUniquePositions))

celset  <- celSetUniquePositions
celsetNormalised <- celSetMATNormalisedUniquePositions

MATSmoothingProcess <- MatSmoothing(celSetMATNormalisedUniquePositions, design = designMatrix, probeWindow = 300, tag = "manual", nProbes = 10)
MATSmoothedScores <- process(MATSmoothingProcess, units = NULL, verbose = -10)
smoothedTilingScores <- MATSmoothedScores


cdfTiling <- AffymetrixCdfFile$byChipType("Hs_PromPR_v02", verbose = verbose)
cdfTableUniquePositions <- getUniqueCdf(cdfTiling)
cdfExpression <- AffymetrixCdfFile$byChipType("HuGene-1_0-st-v1", verbose = verbose)
cdfTableUniquePositions <- getUniqueCdf(cdfExpression)
celsetExpression <- AffymetrixCelSet$byName("geneExpression", cdf = cdfExpression, verbose = verbose)
celsetK27me3Tiling <- AffymetrixCelSet$byName("K27me3Tiling", cdf = cdfTiling, verbose = verbose)
celsetK9AcTiling <- AffymetrixCelSet$byName("K9AcTiling", cdf = cdfTiling, verbose = verbose)
cat("\nReading human genome annotation.")
genePositions <- read.csv("annotationData/humanGenomeAnnotation.csv")
cat("\nAnnotation reading done.")
cat("\nReading sequencing data.")
load("rawData/sequencing/seq_data.Rdata")
cat("\nLoading raw data complete.\n")

