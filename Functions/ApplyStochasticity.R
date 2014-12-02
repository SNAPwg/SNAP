ApplyStochasticity<- function(Life,Fleet,SimCTL)
{
  
  WhereError<- grepl('Error',SimCTL[,2])
 
  ErrorStrings<- strsplit(as.character(SimCTL[WhereError,2])," ") #Grab the second string in each variable containing the term "error"
  
  ErrorsToApply<- ldply(ErrorStrings, function(Err) Err[2])
  
  for (d in 1:dim(ErrorsToApply)[1])
  {
    
    ErrorTerm<- strsplit(as.character(ErrorsToApply[d]),"_")[[1]][1]
    
    ErrorSize<- SimCTL[grepl(ErrorsToApply[d],SimCTL[,2]),1]
    
    IsLifeError<- grepl(ErrorTerm,Life[,2])

    IsFleetError<- grepl(ErrorTerm,Fleets[,2])
    
    if (sum(IsLifeError)>0)
    {
      
      Life[IsLifeError,1]<-  rlnorm(1,log(Life[IsLifeError,1]),ErrorSize)
      
    }
    if (sum(IsFleetError)>0)
    {
      Fleets[IsFleetError,1]<-  rlnorm(1,log(Fleets[IsFleetError,1]),ErrorSize)
    }
     
  }
  
  return(list(Life=Life,Fleets=Fleets))
  
}