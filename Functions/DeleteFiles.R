#' Function to delete output files from ADMB (if they exist) and perform cleanup.
#' @name DeleteFiles
#' @title Delete ADMB output files
#' @param ADMBDir An object of class \code{character} that specifies the 
#'   location of the ADMB executable.
#' @return \code{NULL}
#' @author Adrian Hordyk
#' @seealso \code{\link{}} 
#' @export
#' @examples      
#' \dontrun{
#' 
#' } 
#' 
DeleteFiles <- function(ADMBDir) {
  files <- list.files(ADMBDir)
  delfiles <- grep("admodel", files)
  delfiles <- append(delfiles, grep(".bar", files))
  delfiles <- append(delfiles, grep(".cor", files))
  delfiles <- append(delfiles, grep(".eva", files))
  delfiles <- append(delfiles, grep(".out", files))
  delfiles <- append(delfiles, grep("variance", files))
  delfiles <- append(delfiles, grep(".par", files))
  delfiles <- append(delfiles, grep(".rep", files))
  delfiles <- append(delfiles, grep(".std", files))
  delfiles <- append(delfiles, grep(".b01", files))
  delfiles <- append(delfiles, grep(".p01", files))
  delfiles <- append(delfiles, grep(".r01", files))
  delfiles <- append(delfiles, grep(".log", files))
  delfiles <- append(delfiles, grep(".pin", files))
  delfiles <- append(delfiles, grep(".rpt", files))
  Names <- paste(ADMBDir, files[delfiles], sep="/")
  if (length(delfiles) >= 1) file.remove(Names)
  
}