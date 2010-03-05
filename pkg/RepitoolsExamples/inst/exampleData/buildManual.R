
library(RepitoolsExamples)
setupExamples()
options(prompt = " ", continue = " ")
# might change this
Sweave("~/R/x86_64-pc-linux-gnu-library/2.10/RepitoolsExamples/exampleData/RepitoolsManual.Rnw", quiet=TRUE)
options(prompt = "> ", continue = "+ ")

