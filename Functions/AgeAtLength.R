AgeAtLength<- function(vbk,linf,t0,length)
{
  Age<-  ((log(1-length/linf)/-vbk)+t0)
  
  return(Age)
}