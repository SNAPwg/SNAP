CreateFish <- function()
{

  Fish<- list()

  Fish$Alpha<- 0.5

  Fish$PeakTol<- 0.001 #Wiggle room parameter for finding peaks in catch curve

  Fish$vbk<- NA

  Fish$AgeSD<- 1

  Fish$Linf<- NA

  Fish$t0<- NA

  Fish$WeightA <- NA#Weight at length parameters taken from Kelp Rockfish from Lea et al. 2002/2003, as documented in MLPA model

  Fish$WeightB<- NA

  Fish$MaxAge<- NA

  Fish$M<- NA

  Fish$MvK<- NA

  Fish$AgeMat50<- NA

  Fish$Mat50<- NA

  Fish$Mat95<- NA

  Fish$LengthMatRatio<- NA

  Fish$PLD<- NA

  Fish$VBErrorSlope<- NA

  Fish$GrowthType <- NA

  return(Fish)
}
