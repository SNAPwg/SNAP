AssessLBSPRResiduals<- function(Estimates,Fish,Year)
{
#  Estimates<- runMod
  
 Residuals<- Estimates$Pred - Estimates$Obs
  
 Ages<- floor(AgeAtLength(Estimates$Bins,Fish,0))
  CohortDeviates<- as.data.frame(matrix(NA,nrow=length(Ages),ncol=2))
  
  AgeDeviates<- as.data.frame(matrix(NA,nrow=length(Ages),ncol=2))
  
  for (a in 1:length(Ages))
  {
    Where<- (Ages==Ages[a])
    
    CohortDeviates[a,]<- data.frame(Year-Ages[a],sum(Residuals[Where]))
    
    AgeDeviates[a,]<- data.frame(Ages[a],sum(Residuals[Where]))
    
  }
  
  colnames(CohortDeviates)<- c('Cohort','Residuals')
  
  colnames(AgeDeviates)<- c('Age','Residuals')
  
  return(list(CohortDeviates=CohortDeviates,AgeDeviates=AgeDeviates))
}