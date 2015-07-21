MaturityAtAge <- function(Length,Fish)
{
  
  
  s50<- Fish$Mat50
  
  s95<- Fish$Mat95
  
  mature<- ((1/(1+exp(-log(19)*((Length-s50)/(s95-s50))))))  
}
