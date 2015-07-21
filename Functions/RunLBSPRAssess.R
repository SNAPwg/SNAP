#' This function runs the ADMB assessment code and returns of list of parameter 
#' estimates, model fits, and some other information.
#' @name RunLBSPRAssessFun
#' @title A wrapper function to execute the ADMB optimization routine
#' @param AssessPars An object of class \code{list} that contains the parameters
#'   required for the LB-SPR assessment. Use \code{LoadAssessPars} function to 
#'   create this object.
#' @param LenFreq A \code{numeric} vector that contains the counts of the length
#'   observed data corresponding to the \code{LenMids} length bins Must be the 
#'   same length as \code{LenMids}.
#' @param LenMids A \code{numeric} vector that contains the mid-points of the 
#'   length classes corresponding to the \code{LenFreq} data. Must be the same 
#'   length as \code{LenFreq}.
#' @param ADMBDir An object of class \code{character} that specifies the 
#'   location of the ADMB executable file.
#' @param ExName An object of class \code{character} that specifies the name of 
#'   the ADMB executalbe file. Default is \code{lbspr} and should not include 
#'   extension.
#' @param showOutput An object of class \code{logical} that determines if the 
#'   output of the ADMB will be displayed in the R Console. Apparently only 
#'   works on Windows OS.
#' @param MaxCount An object of class \code{numeric} that specifies the maximum 
#'   number of times the ADMB code is executed with different starting 
#'   parameters before the determining that the model failed to converge. 
#'   Default is 5.
#' @param ADMBRead ADD DETAILS HERE
#' @return A \code{list} that contains: 
#' \itemize{ 
#'   \item{\code{ModelFailed}}{ A \code{logical} that determines if the model 
#'   failed to converge. If #'   \code{TRUE} then other objects in list are 
#'   likely to be \code{NULL} or  missing completely.} 
#'   \item{\code{Estimates}}{ A \code{matrix} containing
#'   the parameter estimates.} 
#'   \item{\code{Pred}}{ A \code{numeric} vector
#'   containing the predicted relative size composition.} 
#'   \item{\code{Obs}}{ A \code{numeric} vector
#'   containing the observed size composition.}
#'   \item{\code{Bins}}{ A \code{numeric} vector containing the length class 
#'   corresonding to observed and predicted length data.} 
#'   \item{\code{ObjFunVal}}{ A \code{numeric} value specifying the value 
#'   of the objective function of the ADMB model.} 
#'   \item{\code{Grad}}{ The final gradient value.} 
#'   \item{\code{Cormat}}{ A \code{matrix} containing the estimated 
#'   variance-covariance matrix.} }
#' @author Adrian Hordyk
#' @seealso \code{\link{}}
#' @export
#' @examples      
#' \dontrun{
#' 
#' }  
#' 

RunLBSPRAssess <- function(AssessPars, LenFreq, LenMids, ADMBDir, ExName="lbspr", showOutput=FALSE, MaxCount=5, ADMBRead=NULL) {

  GetWD <- getwd()
  Count <- 0 
  Output <- NULL 
  if (length(ADMBRead) < 1) ADMBRead <- ADMBDir
  setwd(ADMBRead)
  
  ADMBFile <- file.path(ADMBRead,  ExName)
  DeleteFiles(ADMBDir) # Delete old files if they exist
  WriteDat(AssessPars, LenMids, LenFreq, ADMBRead, ExName) # Write Data file 
  # Determine starting values for optimizer
  Ind <- min(which(cumsum(LenFreq)/max(cumsum(LenFreq)) > 0.4))
  startSL50 <- LenMids[Ind]
  startDelta <- 0.1*AssessPars$Linf
  startFM <- 1
  InitVals <- c(startSL50, startDelta, startFM)
  WritePin(ExName, ADMBDir, InitVals) # Write Pin file
  
  # Run ADMB optimizer  
  ModelFailed <- FALSE
  Sys.chmod(as.character(ADMBFile), mode="7777", use_umask =FALSE)
  
  ADMBCode <- system2(ADMBFile) #, show.output.on.console=showOutput)
  if (ADMBCode > 0) ModelFailed <- TRUE
  TryRead <- try(read.table(paste0(ADMBRead, "/", ExName, ".std"), skip=1)[4:7, 2:4])
  if (class(TryRead) == "try-error")  
  { 
    ModelFailed <- TRUE
  }
  
  while (ModelFailed & Count <=MaxCount) {
    Ind <- min(which(cumsum(LenFreq)/max(cumsum(LenFreq)) > runif(1)))
    startSL50 <- LenMids[Ind]
    startDelta <- runif(1, min=0.1*AssessPars$Linf, max=0.6*AssessPars$Linf)
    startFM <- runif(1, min=0.1, max=5)
    WritePin(ExName, ADMBDir, InitVals) # Write Pin file
    ModelFailed <- FALSE
    ADMBCode <- system(ADMBFile)
    if (ADMBCode > 0) ModelFailed <- TRUE
    TryRead <- try(read.table(paste0(ADMBRead, "/", ExName, ".std"), skip=1)[4:7, 2:4])
    if (class(TryRead) == "try-error") {
      Count <- Count + 1  
      ModelFailed <- TRUE
    }  
    print(Count)
  }
  
  if (ModelFailed) {
    setwd(GetWD)
    Output$ModelFailed <- TRUE
    return(Output)
  }	
  
  if (ModelFailed == FALSE) {
    Estimates <- TryRead
    colnames(Estimates) <- c("Par", "Est", "SD")
    
    # Read in Report file
    modoutput <- scan(paste0(ADMBRead, "/", ExName, ".rep"), quiet=TRUE) 
    NLBins <- modoutput[1] 
    Pred <- modoutput[2:(NLBins+1)]
    Obs <- modoutput[(NLBins+2):(NLBins+1+NLBins)]
    Bins <- modoutput[(NLBins+2+NLBins):(NLBins+NLBins+NLBins+1)]
    Unfished <- modoutput[(NLBins+2+NLBins+NLBins):(NLBins+NLBins+NLBins+NLBins+1)] 
    ObjFunVal <- modoutput[length(modoutput)-1]
    Grad <- modoutput[length(modoutput)]
    # Read in Correlation
    # Borrowed from R2admb package and slightly modified
    AllStd <- read.table(paste0(ADMBRead, "/", ExName, ".std"), skip=1)
    Names <- AllStd[,2]
    nsdpar <- nrow(AllStd)
    if (file.exists(paste(ExName, "cor", sep = "."))) {
      ncorpar <- length(readLines(paste(ExName, "cor", sep = "."))) -   2
      cor_dat <- read.table(paste0(ExName, ".cor"), skip = 2, fill = TRUE, as.is = TRUE, col.names = paste("X", 1:(4 + ncorpar), sep = ""))
      cormat <- as.matrix(cor_dat[4:nsdpar, 4 + (4:nsdpar)])
      cormat[upper.tri(cormat)] <- t(cormat)[upper.tri(cormat)]
      colnames(cormat) <- Names[4:nsdpar]
      rownames(cormat) <- Names[4:nsdpar]
    } else {
      warning(".cor file not found")
      cormat <- NA
    } 
    
    setwd(GetWD)
    Output$ModelFailed <- ModelFailed
    Output$ADMBCode <- ADMBCode
    Output$Estimates <- Estimates
    Output$Pred <- Pred
    Output$Obs <- Obs
    Output$Bins <- Bins
    Output$Unfished <- Unfished
    Output$ObjFunVal <- ObjFunVal
    Output$Grad <- Grad
    Output$Cormat <- cormat
    return(Output)
  }	
}
