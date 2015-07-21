#' Function to simulate length-structured growth-type-group (GTG) model to 
#' generate size equilibrium composition of population and catch, as well as SPR
#' of stock and relative yield. This model is follows the same length-structured GTG 
#' approach of \code{SimMod_LHR} but uses the rate parameters directly (i.e M, K, F)
#' rather than the ratios. The results from the two models should be identical if the 
#' life history and fishing ratios are the same.
#' @name SimMod_AgeStruc
#' @title Generate size structure using GTG model and calculate relative YPR and SPR
#' @param SimPars An object of class \code{list} that contains all parameters
#'   required to run GTG model.  Full description of model to be added at later
#'   date.
#' @param Mpar The natural mortality rate. As this model is age-structured it requires
#'   a value for the mean M of the stock. The other rate parameters (K and F) are calculated
#'   from the ratios provided in \code{SimPars}.
#' @param kslope An object of class \code{numeric} which determines the slope of the 
#'   the natural mortality for each GTG to approximate equal fitness across GTG. 
#'   A value of 0 means all GTG have the same natural mortality. The value must be 
#'   quite small to ensure that M is not negative for any groups. See manuscript for 
#'   more details.
#' @return  To add details later.
#' @author Adrian Hordyk
#' @export
 

SimMod_AgeStruc <- function(SimPars, Mpar,...) {
  with(SimPars, {
  Fpar <- Mpar * FM
  kpar <- Mpar/ MK
  
  SDLinf <- CVLinf * Linf # Standard Deviation of Pop. Linf 
  
  # Growth-Type-Group Model Setup 
  DiffLinfs <- seq(from=Linf - MaxSD * SDLinf, to=Linf + MaxSD * SDLinf, length=NGTG)
  GTGLinfdL <- DiffLinfs[2] - DiffLinfs[1]
  RecProbs <- dnorm(DiffLinfs, Linf, SDLinf)/sum(dnorm(DiffLinfs, Linf, SDLinf)) # Recruits normally distributed across GTGs
  
  # Set up Length Bins of Population Model 
  LenBins <- seq(from=0, by=Linc, to=Linf + MaxSD * SDLinf)
  LenMids <- seq(from=LenBins[1] + 0.5*Linc, by=Linc, length=length(LenBins)-1)
  
  MatLen <- 1.0/(1+exp(-log(19)*(LenMids-L50)/(L95-L50))) # Maturity Schedule 
  Weight <- Walpha * LenMids^Wbeta # Weight-at-length 
  FecLen <- MatLen * LenMids^FecB # Relative Fecundity at Length 
  FecLen <- FecLen/max(FecLen) # Divide by maximum to make it more useful for plotting. Can add extra parameter(s) so that fecundity is absolute rather than relative (makes no difference to SPR)
  
  VulLen <- 1.0/(1+exp(-log(19)*(LenBins-SL50)/(SL95-SL50))) # Selectivity-at-Length
  SelLen <- VulLen
  # # Minimum Legal Length 
  # if (length(MLL) > 0 & MLL > 0 ) {
  # Legal <- rep(0, length(LenBins))
  # Legal[LenBins >= MLL] <- 1 
  # SelLen <- VulLen * Legal
  # } 
  
  MVec <- Mpar * (Linf/(LenBins+0.5*Linc))^Mpow # Vector of natural mortality for mean GTG
  MMat <- matrix(MVec, nrow=length(MVec), ncol=NGTG) # Matrix of M for each GTG
  tempFun <- function(X) MVec + kslope*(DiffLinfs[X] - Linf)
  MMat <- sapply(seq_along(DiffLinfs), function (X) tempFun(X))
  
  FVec <- Fpar * SelLen # Vector of F for each size class
  ZMat <- MMat + FVec # Matrix of total mortality for each GTG
  
  
  # Set Up Empty Matrices 
  NPRFished <- NPRUnfished <- matrix(0, nrow=length(LenBins), ncol=NGTG) # number-per-recruit at length
  NatLUnFishedPop <- NatLFishedPop <- NatLUnFishedCatch <- NatLFishedCatch <- FecGTGUnfished <- matrix(0, nrow=length(LenMids), ncol=NGTG) # number per GTG in each length class 
    
  EPR_GTG <- matrix(NA, nrow=NGTG, ncol=2)
  YPR <- rep(NA, NGTG)
 
 # Loop over Growth-Type-Groups 	
  for (GTG in 1:NGTG) {
    NPRUnfished[1, GTG] <- NPRFished[1, GTG] <- RecProbs[GTG]
	options(warn=-1) # turn off warnings
    ALen <- -(log(1-LenBins/DiffLinfs[GTG]))/kpar # Calculate Age for length
	options(warn=0) # turn back on
    ALen[is.na(ALen >= 0) ] <- NA # Get rid of NaNs 
	HighAge <- which.max(ALen)
	LastInc <- ALen[HighAge-1] - ALen[HighAge-2]
	if (ALen[HighAge] != Inf) {
	  LastInc <- ALen[HighAge] - ALen[HighAge-1]
      ALen[HighAge+1] <- ALen[HighAge] + 2 * LastInc # Make a maximum age for last length class 
	}    
    if (ALen[HighAge] == Inf) {
	  LastInc <- ALen[HighAge-1] - ALen[HighAge-2]
	  ALen[HighAge] <- ALen[HighAge-1] + 2 * LastInc # Make a maximum age for last length class 
	}
	for (L in 2:(HighAge)) {
      if (LenBins[L] < DiffLinfs[GTG]) {
        NPRUnfished[L, GTG] <- NPRUnfished[L-1, GTG] * exp(-MMat[L-1, GTG] * (ALen[L] - ALen[L-1]))
        NPRFished[L, GTG] <- NPRFished[L-1, GTG] * exp(-ZMat[L-1, GTG] * (ALen[L] - ALen[L-1]))
      }
    }
    for (L in 1:length(LenMids)) {
      NatLUnFishedPop[L, GTG] <- (NPRUnfished[L,GTG] - NPRUnfished[L+1,GTG])/MMat[L, GTG]
      NatLFishedPop[L, GTG] <- (NPRFished[L,GTG] - NPRFished[L+1,GTG])/ZMat[L, GTG]  
	  FecGTGUnfished[L, GTG] <- NatLUnFishedPop[L, GTG] * FecLen[L]
    }	

    NatLUnFishedCatch[,GTG] <- NatLUnFishedPop[, GTG] * 1.0/(1+exp(-log(19)*(LenMids-SL50)/(SL95-SL50))) 
    NatLFishedCatch[,GTG]  <- NatLFishedPop[, GTG] * 1.0/(1+exp(-log(19)*(LenMids-SL50)/(SL95-SL50)))
	   
    # Eggs-per-recruit for each GTG 
    EPR_GTG[GTG,1] <- sum(NatLUnFishedPop[, GTG] * FecLen)
    EPR_GTG[GTG,2] <- sum(NatLFishedPop[, GTG] * FecLen)
      
    # YPR 
    YPR[GTG] <- sum(NatLFishedPop[, GTG]  * Weight * 1.0/(1+exp(-log(19)*(LenMids-SL50)/(SL95-SL50)))) * FM
  }	
	
  # Calc Unfished Fitness 
  Fit <- apply(FecGTGUnfished, 2, sum, na.rm=TRUE) # Total Fecundity per Group
  FitPR <- Fit/RecProbs # Fitness per-recruit
  ObjFun <- sum((FitPR - median(FitPR, na.rm=TRUE))^2, na.rm=TRUE) # This needs to be minimised to make fitness approximately equal across GTG - by adjusting kslope
  Pen <- 0; if (min(MMat) <= 0 ) Pen <- (1/abs(min(MMat)))^2 * 1E5 # Penalty for optimising kslope   
  ObjFun <- ObjFun + Pen
  # print(cbind(kslope, ObjFun, Pen))
  
  # Calc SPR
  EPR0 <- apply(EPR_GTG,2,sum)[1]
  EPRf <- apply(EPR_GTG,2,sum)[2]
  SPR <- EPRf/EPR0  
  
  # Equilibrium Relative Recruitment 
  reca <- recK/EPR0
  recb <- (reca * EPR0 - 1)/(R0*EPR0)
  RelRec <- max(0, (reca * EPRf-1)/(recb*EPRf))
    
  Yield <- sum(YPR) * RelRec

  # Expected Length Structure of Catch 
  ExpectedLenCatchFished <- apply(NatLFishedCatch, 1, sum)/sum(apply(NatLFishedCatch, 1, sum))
  ExpectedLenPopFished <- apply(NatLFishedPop, 1, sum)/sum(apply(NatLFishedPop, 1, sum))
  ExpectedLenCatchUnfished <- apply(NatLUnFishedCatch, 1, sum)/sum(apply(NatLUnFishedCatch, 1, sum))
  ExpectedLenPopUnfished <- apply(NatLUnFishedPop, 1, sum)/sum(apply(NatLUnFishedPop, 1, sum))
  
  Output <- NULL 
  Output$SPR <- SPR
  Output$Yield <- Yield 
  Output$ExpLenCatchFished <- ExpectedLenCatchFished
  Output$ExpLenPopFished <- ExpectedLenPopFished
  Output$ExpLenCatchUnfished <- ExpectedLenCatchUnfished
  Output$ExpLenPopUnfished <- ExpectedLenPopUnfished
  Output$LenBins <- LenBins
  Output$LenMids <- LenMids
  Output$NGTG <- NGTG
  Output$GTGdL <- GTGLinfdL
  Output$DiffLinfs <- DiffLinfs
  Output$RecProbs <- RecProbs
  Output$Weight <- Weight
  Output$FecLen <- FecLen 
  Output$MatLen <- MatLen 
  Output$SelLen <- SelLen
  Output$MMat <- MMat 
  Output$FVec <- FVec 
  Output$ZMat <- ZMat
  Output$ObjFun <- ObjFun 
  Output$Pen <- Pen
  Output$Mpar <- Mpar 
  Output$kpar <- kpar 
  Output$FitPR <- FitPR
  Output$Diff <- range(FitPR)[2] - range(FitPR)[1]  
  return(Output)
  })
}
