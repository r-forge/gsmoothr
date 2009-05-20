.onLoad <- function(libname, pkgname) {
  s <- search() 
  library.dynam("epitools",pkgname,libname,now=FALSE)
}

