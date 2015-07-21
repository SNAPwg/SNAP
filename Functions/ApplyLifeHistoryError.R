ApplyLifeHistoryError<- function()
{
  
  #   LHI_Error <- 1.1
  
  #   while(LHI_Error>Fish$LHITol) # Only accept values that don't violate LHI too much
  #   {
  NewFish<- Fish
  
  NewFish$vbk<- Fish$vbk*rlnorm(1,0,Fish$LengthError)
  
  NewFish$MvK<- runif(1,NewFish$MinMvK,NewFish$MaxMvK)
  
  NewFish$M<- NewFish$vbk * NewFish$MvK
  
  #     NewFish$M<-rlnorm(1,log(Fish$M),(Fish$MortalityError))
  
  NewFish$Linf<- Fish$Linf*rlnorm(1,0,Fish$LengthError)
  
  NewFish$LengthMatRatio<- runif(1,NewFish$MinLengthMatRatio,NewFish$MaxLengthMatRatio)
  
  NewFish$Mat50<- NewFish$Linf * NewFish$LengthMatRatio
  
  NewFish$Mat95<- NewFish$Mat50*1.05    
  
  #     if (Fish$t0!=0)
  #     {
  #       NewFish$t0<- rnorm(1,(Fish$t0),(Fish$LengthError))
  #     }
  #       
  
  NewFish$MaxAge<- -log(.01)/Fish$M
  
  #     Mtm_Error<- abs((NewFish$M*AgeAtLength(NewFish$Mat95,Fish,0))/1.65-1)
  #     
  #     MvK_Error<- abs((NewFish$M/NewFish$vbk)/1.6-1)
  #     
  #     LmvLinf_Error<- abs((mean(c(NewFish$Mat50,NewFish$Mat95))/NewFish$Linf)/0.67-1)
  #     
  #     LHI_Error<- (mean(c(Mtm_Error,MvK_Error, LmvLinf_Error)))
  
  #   }
  
  return(NewFish)
}