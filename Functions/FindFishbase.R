######################################
#UseFishbase--------------------------------------------------
#This code assigns fishbase based life history variables to any data missing life history variables
######################################

FindFishbase<- function(Data)
{
    
#    Data<- GFD
  load('Data/mpack.Rdata')
  
 Fishbase<- mpack$lh

 Fishbase[is.na(Fishbase)]= NA
 rm(mpack)
 
 spnames<- unique((Data$SciName[is.na(Data$SciName)==F]))
 
for (i in 1:length(spnames)){
  
  show(paste(100*(round(i/length(spnames),2)),'% Done With Fishbase Matching',sep=''))
  
  lh_match=grepl(spnames[i],Fishbase$sname) #find if species name name matches FB data, locate row
 
  if (sum(lh_match)>0)
  {
  
    WhereSpecies<- spnames[i]==Data$SciName
#   WhereVb<- spnames[i]==Data$SciName & is.na(Data$VonBertK)
  WhereVb<- WhereSpecies & is.na(Data$vbk)
  
  WhereAgeMat<-  WhereSpecies & is.na(Data$AgeMat50)
  WhereLinf<- WhereSpecies & is.na(Data$Linf)
#   WhereTemp<- WhereSpecies & is.na(Data$Temp)
  
    lh_data=Fishbase[lh_match,] #pull out matching data
#     spcat=lh_data$spcat #identify the species category
    
    #Add in life history
    Data$Linf[WhereLinf]<- lh_data$maxl[1]
    Data$LinfSource[WhereLinf]<-"FishBase"
    
     Data$vbk[WhereVb]<-lh_data$vbk[1]
    Data$vbkSource[WhereVb]<-"FishBase"
    Data$AgeMat50[WhereAgeMat]<-lh_data$agem[1]
    Data$AgeMatSource[WhereAgeMat]<-"FishBase"
 } #close if loop

} #close species loop

return(Data)

}