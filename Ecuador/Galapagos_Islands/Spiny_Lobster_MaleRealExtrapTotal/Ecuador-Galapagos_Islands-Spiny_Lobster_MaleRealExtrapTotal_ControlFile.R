#### ControlFile ###



AvailableData<- c("LengthData",'CatchData')

 Assessments<- c('DBSRA','CatchMSY','LBAR','CatchCurve','LBSPR')

 # Assessments<- c('LBSPR')

# Assessments<- c('CatchCurve')

### Life History ###

Fish<- NULL

Fish$SciName<- 'P.Gracilatis'
Fish$CommName<- 'Spiny Lobster'


DefaultSD<- 0.05

Fish$LHITol<- 0.1 # The average % deviation from LHI allowed

Fish$vbk<- 0.21625                        ##according to Informe final   0.21625    error 0.004M and 0.06H

Fish$LengthError<- 0.05

Fish$Linf<- 41.9                        ##41.9M and 40.14H averaged from Norahs values (used averages = 41.02) used our equation to determine Linf for tail length (average = 25.15 and H = 24.57)

Fish$t0<- 0

Fish$WeightA <- 0.0340/1000 #From D. Viana's summary of lit, "Lobster Life History..."

Fish$WeightB<- 3 #From D. Viana's summary of lit, "Lobster Life History..."

# Fish$MaxAge<- 25

#Fish$M<- 0.17                   ## average mortality according to Matt Kay and literature review, Hearn 2008 use 0.348M and 0.378H (average = 0.363)

 Fish$M<- 0.348  

# Fish$M<- 0.166  


Fish$MvK<- 1.6 #Ratio of M versus K

Fish$MaxAge<- log(.01)/-Fish$M #Max age is age at thich only 1% of population are left. Can also be set manually if needed

Fish$MortalityError<- 0.1                 ## test error, can change to accomodate more error

Fish$Mat50<- 24.05                        ## 24.05 according to Toral and Hearn (tail only = 13.96)

Fish$Mat95<- 26                              ## 26 according to Toral and Hearn (tail only = 15.24)

Fish$PLD<- 31

Fish$VBSD<- 0.05

Fish$VBErrorSlope<- 0.1

Fish$res<- 'High'

### Fleet Parameters ###

Fleet<- NULL

Fleet$MinSizeCaught<- 26                        ##26 cm minimum landing size for Galapagos      (was 20 before?) (tail = 15.24) 

Fleet$MaxSizeCaught<- Fish$Linf-0.5