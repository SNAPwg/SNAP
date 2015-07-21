#' Function to simulate length-structured growth-type-group (GTG) model to 
#' generate size equilibrium composition of population and catch, as well as SPR
#' of stock and relative yield. This model is parameterised with the life-history ratios.
#' @name SimMod_LHR
#' @title Generate size structure using GTG model and calculate relative YPR and SPR
#' @param SimPars An object of class \code{list} that contains all parameters
#'   required to run GTG model.  Full description of model to be added at later
#'   date.
#' @param kslope An object of class \code{numeric} which determines the slope of the 
#'   the natural mortality for each GTG to approximate equal fitness across GTG. 
#'   A value of 0 means all GTG have the same natural mortality. The value must be 
#'   quite small to ensure that M is not negative for any groups. See manuscript for 
#'   more details.
#' @return  To add details later.
#' @author Adrian Hordyk
#' @seealso \code{\link{}}
#' @export
 

SimMod_LHR <- function(SimPars, ...) {
  with(SimPars, {
    SDLinf <- CVLinf * Linf # Standard Deviation of Pop. Linf 
    
    # Error Catches #
    
    # if (CVLinf * Linf * MaxSD * MK * ) add all parameters and check that greater than zero
    # if (length(stuff) == 0)
    
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
    
    # Life-History Ratios  
    MKL <- MK * (Linf/(LenBins+0.5*Linc))^Mpow # M/K ratio for each length class
    MKMat <- matrix(MKL, nrow=length(MKL), ncol=NGTG) # Matrix of MK for each GTG
	tempFun <- function(X) MKL + kslope*(DiffLinfs[X] - Linf)
    MKMat <- sapply(seq_along(DiffLinfs), function (X) tempFun(X))
	
    FK <- FM * MK # F/K ratio 
    FKL <- FK * SelLen # F/K ratio for each length class   
    # FkL[Legal == 0] <- FkL[Legal == 0] * DiscardMortFrac 
    ZKLMat <- MKMat + FKL # Z/K ratio (total mortality) for each GTG
    
	
    # Set Up Empty Matrices 
    NPRFished <- NPRUnfished <- matrix(0, nrow=length(LenBins), ncol=NGTG) # number-per-recruit at length
    NatLUnFishedPop <- NatLFishedPop <- NatLUnFishedCatch <- NatLFishedCatch <- FecGTGUnfished <- matrix(0, nrow=length(LenMids), ncol=NGTG) # number per GTG in each length class 
    
    EPR_GTG <- matrix(NA, nrow=NGTG, ncol=2)
    YPR <- rep(NA, NGTG)
    
    # Loop over Growth-Type-Groups 
    for (GTG in 1:NGTG) {
      NPRFished[1, GTG] <- NPRUnfished[1, GTG] <- RecProbs[GTG] # Distribute Recruits into first length class 
      GTGLinf <- DiffLinfs[GTG] # Linf for specific GTG 
      for (L in 2:length(LenBins)) {
        if (LenBins[L] < GTGLinf) {
          NPRUnfished[L, GTG] <- NPRUnfished[L-1, GTG] * ((GTGLinf-LenBins[L])/(GTGLinf-LenBins[L-1]))^MKMat[L-1, GTG]
          NPRFished[L, GTG] <- NPRFished[L-1, GTG] * ((GTGLinf-LenBins[L])/(GTGLinf-LenBins[L-1]))^ZKLMat[L-1, GTG]
        }
      }      
      for (L in 1:length(LenMids)) {
        NatLUnFishedPop[L, GTG] <- (NPRUnfished[L,GTG] - NPRUnfished[L+1,GTG])/MKMat[L, GTG]
        NatLFishedPop[L, GTG] <- (NPRFished[L,GTG] - NPRFished[L+1,GTG])/ZKLMat[L, GTG]  
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
	Pen <- 0; if (min(MKMat) <= 0 ) Pen <- (1/abs(min(MKMat)))^2 * 1E5 # Penalty for optimising kslope   
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
	Output$MKL <- MKL
    Output$MKMat <- MKMat 
    Output$FKL <- FKL 
    Output$ZKLMat <- ZKLMat 
	Output$ObjFun <- ObjFun 
	Output$Pen <- Pen
	Output$FitPR <- FitPR
	Output$Diff <- range(FitPR)[2] - range(FitPR)[1]
    
    return(Output)
  })
}

