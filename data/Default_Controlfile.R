#### ControlFile ###

# Assessments<- c('CatchCurve')

### Life History ###

Fish<- NULL

Fish$Alpha<- 0.5

Fish$PeakTol<- 0.001 #Wiggle room parameter for finding peaks in catch curve

Fish$vbk<- NA

# Fish$LengthError<- DefaultSD

Fish$AgeSD<- 1


Fish$Linf<- NA #ish, 120mm SL

Fish$t0<- NA

Fish$WeightA <- 9.37e-6 #Weight at length parameters taken from Kelp Rockfish from Lea et al. 2002/2003, as documented in MLPA model

Fish$WeightB<- 3.172

Fish$MaxAge<- NA

Fish$M<- 0.2

Fish$MvK<- NA

# Fish$MaxAge<- log(.01)/-Fish$M #Max age is age at thich only 1% of population are left. Can also be set manually if needed

# Fish$MaxAge<- log(.01)/-Fish$M #Max age is age at thich only 1% of population are left. Can also be set manually if needed

# Fish$MortalityError<- DefaultSD

Fish$AgeMat50<- NA

Fish$Mat50<- NA

Fish$Mat95<- NA

Fish$LengthMatRatio<- NA

Fish$PLD<- 31

# Fish$VBSD<- DefaultSD

Fish$VBErrorSlope<- 0.1


### Fleet Parameters ###

Fleet<- NULL

Fleet$MinSizeCaught<- 20

Fleet$MaxSizeCaught<- Fish$Linf-0.5