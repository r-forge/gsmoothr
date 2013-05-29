############################################################
##
## file: init_fns.R
##
## Mark Robinson (mrobinson@wehi.edu.au)
## Created 07 Oct 2007
##
##
############################################################

.First.lib <- function(libname, pkgname) {
  s <- search() 
  library.dynam("FIRMAGene",pkgname,libname,now=FALSE)
}
