BuildLBSPRPars<- function(Fish,Mpow,NGTG,MaxSD,FecB,SL50Min,SL50Max,DeltaMin,DeltaMax)
{
  with(Fish,{
  MK     <- MvK
  Linf   <- Linf
  CVLinf <- 0.1
  L50    <- Mat50
  L95    <- Mat95
  Walpha <- WeightA
  Wbeta  <- WeightB
  FecB   <- FecB
  Mpow   <- Mpow
  NGTG   <- NGTG
  SDLinf <- CVLinf * Linf
  MaxSD  <- MaxSD
  GTGLinfdL <- ((Linf + MaxSD * SDLinf) - (Linf - MaxSD * SDLinf))/(NGTG-1) 
  
  SL50Min <- SL50Min
  if (length(SL50Min) < 1 | is.na(SL50Min)) SL50Min <- 1
  SL50Max <- SL50Max
  if (length(SL50Max) < 1| is.na(SL50Max)) SL50Max <- 0.95 * Linf
  DeltaMin <- DeltaMin
  if (length(DeltaMin) < 1| is.na(DeltaMin)) DeltaMin <- 0.01
  DeltaMax <- DeltaMax
  if (length(DeltaMax) < 1| is.na(DeltaMax)) DeltaMax <- 0.5 * Linf
  
  AssessPars <- list(MK=MK, Linf=Linf, CVLinf=CVLinf, L50=L50, L95=L95, 
                     Walpha=Walpha, Wbeta=Wbeta, FecB=FecB, Mpow=Mpow, 
                     NGTG=NGTG, GTGLinfdL=GTGLinfdL, MaxSD=MaxSD, SL50Min=SL50Min, 
                     SL50Max=SL50Max, DeltaMin=DeltaMin, DeltaMax=DeltaMax)
  AssessPars$kslope <- as.numeric(PredictKSlope(AssessPars))
  return(AssessPars)
})  
}