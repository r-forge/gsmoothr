\name{laneBlocksCounts}
\alias{laneBlocksCounts}

\title{Finds number of reads around TSS sites from sequencing lane data.}
\description{
Reads in all the .fastq files in the lanesPath directory and finds the count of reads in an upstream and downstream range from every gene's TSS.
}
\usage{
laneBlocksCounts(lanesPath, TSSDataTable, upstream = 1000, downstream = 1000, removeDuplications=FALSE, doQuality=FALSE, aligner)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lanesPath}{Path to where the .fastq files of reads are.}
  \item{TSSannotationPath}{A data frame containing the TSS details of the organism of interest. Column names must be name, chr, start, stop, strand.}
  \item{upstream}{Bases upstream from TSS to start the regions.}
  \item{downstream}{Bases downstream from TSS to end the regions.}
  \item{removeDuplications}{Boolean. Should alignments also be processed, where duplicate reads are firstly eliminated?}
  \item{doQuality}{Boolean. Should .html reports be generated about the quality of reads and of alignments.}
  \item{aligner}{Either "Bowtie", "MAQ", or "SolexaExport". The program that created the alignments.}
}
\details{
	The function makes a list of IRanges of the TSS starts of genes. Then, IRanges of singly-mapped reads, for each chromosome, are read. The alignments for each lane are to be located in a sub-directory of\code{lanesPath} with the sub-directory name being the name of the particular read filename (without the extension). For each TSS, the count of mappable reads overlapping with the TSS regions are made into a table, and this is written to the \code{lanesPath} directory. If \code{doQuality} is FALSE, one file called "counts.txt" is written. If \code{doQuality} is TRUE, the file "countsNoDuplication.txt" is additionally written.
	
	The .html quality reports for each lane are written to the directory for that lane.
}
\value{
A list with elements
\item{readCounts}{Matrix of counts. Rows are gene TSS regions. Columns are sequencing lanes.}
\item{readConutsNoDuplications}{Matrix of counts for non-duplicated reads. Rows are gene TSS regions. Columns are sequencing lanes.}
}
\author{Dario Strbenac}
\seealso{}
\examples{
	laneBlocksCounts(lanes, humanTSSTable, aligner="Bowtie")
}
