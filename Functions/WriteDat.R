#' A function to write the.dat file required by the ADMB code.
#' @name WriteDat
#' @title Write data file (.dat) for ADMB
#' @param AssessPars An object of class \code{list} that contains the neccessary
#'   parameters. Created by \code{LoadAssessPars} function.
#' @param LengthMids A vector containing the mid-points of the length bins.
#'   Created by \code{MakeLengthComp}.
#' @param LenFreq A vector containing the observed number in each length bins.
#'   Created by \code{MakeLengthComp}.
#' @param ADMBDir An object of class \code{character} that specifies the
#'   location of the ADMB executable.
#' @param ExName An object of class \code{character} that specifies the name of
#'   the ADMB exectuable (without extension).
#' @return \code{NULL}
#' @author Adrian Hordyk
#' @seealso \code{\link{}}
#' @export
#' @examples      
#' \dontrun{
#' 
#' }  
#' 
WriteDat <- function(AssessPars, LengthMids, LenFreq, ADMBDir, ExName="lbspr") {

  with(AssessPars, {  			
    ChkFile <- file.exists(paste0(ADMBDir, "/", ExName, ".dat"))
    if (ChkFile) file.remove(paste0(ADMBDir, "/", ExName, ".dat"))
    
    # Error handling
    if (any(!is.numeric(c(MK, Linf, CVLinf, L50, L95, 
                          Walpha, Wbeta)))) {
      Error <- "Error - Input parameters not numeric"
      stop(Error)				 
    }
    
    CheckBy <- NULL
    for (X in 2:length(LengthMids)) {
      CheckBy[X-1] <- round(LengthMids[X] - LengthMids[X-1],2)
    }
    chk <- sum(CheckBy) == length(CheckBy)
    if (chk) By <- CheckBy[1]
    if (chk == FALSE ) {
      Error <- "Error - Length classes not equidistant"
      stop(Error)
    }	  
    
    # Preliminary calculations.
    Linc <- By
    NLenMids <- length(LengthMids)
    LenBins <- seq(from=LengthMids[1]-0.5*Linc, by=Linc, length=length(LengthMids)+1)
    
    # To add - check if length bins start at zero 
    SDLinf <- CVLinf * Linf
    NLenBins <- length(LengthMids)
    
    # Write out data file
    con=file(paste0(ADMBDir, "/", ExName, ".dat"), open="wt")
    write(as.character("#M/K"),con)
    write(MK,con)
    write(as.character("#Linf"),con)
    write(Linf,con)
    write(as.character("#CVLinf"),con)
    write(CVLinf,con)
    write(as.character("#NGTG"),con)
    write(NGTG,con)
    write(as.character("#MaxSD"),con)
    write(MaxSD,con)    
    write(as.character("#NLenMids"),con)
    write(NLenMids,con)  
    write(as.character("#LenMids"),con)
    write(LengthMids,con)  
    write(as.character("#ObsLength"),con)
    write(LenFreq,con)
    write(as.character("#LenBins"),con)
    write(LenBins,con)  	
    write(as.character("#L50"),con)
    write(L50,con) 
    write(as.character("#L95"),con)
    write(L95,con) 	
    write(as.character("#FecB"),con)
    write(FecB,con) 	
    write(as.character("#Mpow"),con)
    write(Mpow,con)
    write(as.character("#logSL50Min"),con)
    write(log(SL50Min),con)	
    write(as.character("#logSL50Max"),con)
    write(log(SL50Max),con)
    write(as.character("#logDeltaMin"),con)
    write(log(DeltaMin),con)	
    write(as.character("#logDeltaMax"),con)
    write(log(DeltaMax),con)
    write(as.character("#kslope"),con)
    write(kslope,con)	
    close(con)
  })
} 
