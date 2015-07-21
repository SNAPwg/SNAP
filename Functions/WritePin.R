#' A function to write the .pin value (initials values) for ADMB.
#' @name WritePin
#' @title Write .pin file
#' @param ExName An object of class \code{character} that specifies the name of 
#'   the ADMB exectuable (without extension).
#' @param ADMBDir An object of class \code{character} that specifies the 
#'   location of the ADMB executable.
#' @param Vals An object of class \code{numeric} of length 3 that specifies the
#'   starting values.
#' @return \code{NULL} if successfully exectuted and .pin file written. Otherwise error message.
#' @author Adrian Hordyk 
#' @seealso \code{\link{}}
#' @export
#' @examples      
#' \dontrun{
#' 
#' } 
#' 
WritePin <- function(ExName, ADMBDir, Vals) {
  ChkFile <- file.exists(paste0(ADMBDir, "/", ExName, ".pin"))
  if (ChkFile) file.remove(paste0(ADMBDir, "/", ExName, ".pin"))
  
  if (length(Vals) != 3) {
    Error <- "Starting parameter vector of incorrect length"
    end (Error)
  }
  
  con=file(paste0(ADMBDir, "/", ExName, ".pin"),open="wt")
  write("#logSelL50",con)
  write(log(Vals[1]),con)
  write("#logDelta",con)
  write(log(Vals[2]),con)
  write("# logFM",con)
  write(log(Vals[3]),con)
  close(con)
}
