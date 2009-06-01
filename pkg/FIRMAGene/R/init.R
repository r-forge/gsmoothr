.onLoad <- function(libname, pkgname) {
  s <- search() 
  library.dynam("FIRMAGene",pkgname,libname,now=FALSE)
}

