# Environment, package-version and shipped-input I/O only; no statistical computation.
optional_library <- Sys.getenv('ENDOMETRIAL_R_LIBRARY',unset='')
if(nzchar(optional_library)).libPaths(c(optional_library,.libPaths()))
expected <- c(limma='3.68.5',statmod='1.5.2',jsonlite='2.0.0',digest='0.6.39')
for(p in names(expected)){
 if(!requireNamespace(p,quietly=TRUE))stop('Install required package: ',p)
 if(as.character(packageVersion(p))!=expected[[p]])stop('Reference run requires ',p,' ',expected[[p]],'; found ',packageVersion(p))
}
if(as.character(getRversion())!='4.6.1')warning('Reference R version is4.6.1; validate all numerical comparisons on this R version.')
verify_inputs <- function(){
 manifest <- read.csv(file.path(root,'SEMANTIC_INPUT_MANIFEST.csv'),colClasses='character')
 for(i in seq_len(nrow(manifest))){
  path <- file.path(root,manifest$path[i])
  if(!file.exists(path)||file.info(path)$size!=as.numeric(manifest$bytes[i])||digest::digest(file=path,algo='sha256')!=manifest$sha256[i])stop('Shipped input integrity failure: ',manifest$path[i])
 }
 invisible(TRUE)
}
