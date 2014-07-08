Recruitment<-function(SpawningBiomass,steepness,R0,VirSpBio,detRec,sigmaR)
{
 Recruits<-(0.8*R0*steepness*SpawningBiomass) / ((0.2*VirSpBio*(1-steepness))+ (steepness-0.2)*SpawningBiomass)
 if(detRec==0)
  Recruits<-Recruits*exp(0,sigmaR)
 return(Recruits)
}