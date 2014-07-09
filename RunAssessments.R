
# Assessment and Analysis Modules -----------------------------------------

# This section processes data from the research modules and runs available assessment and/or analysis techniques to guide the management module

library(lattice)
source("AssessmentModules.R") #Pull in assessment modules
source("SubFunctions.R") #Pull in helper functions for assessment modules


#Prepare Data for Assessment --------------------------------------------------

#Format available data from each fleet, species, for assessment

#Figure out which assessments can be run for each unique conbination of fleet and species.Using ths assessment data needs csv, figure out what data are needed
#Then, check to see whether those data exist. If yes, that assessment gets a 1, if not a 0, etc. 

#For each fleet, for each species, for each potential assessments

Fishery <- paste(Country,Site,Species,sep='-')

Directory<- paste(Country, "/", Site,'/',Species,'/',sep='')

source(paste(Country, "/", Site,'/',Species,'/',Fishery,"_ControlFile.R", sep = ""))

dir.create(paste(Directory,'Results',sep=''))

dir.create(paste(Directory,'Figures',sep=''))

FigureFolder<- paste(Directory,'Figures/',sep='')

ResultFolder<- paste(Directory,'Results/',sep='')

for (d in 1:length(AvailableData)) #Read in available data
{
  eval(parse(text=paste(AvailableData[d],'<- read.csv(',"'",Directory,Fishery,'_',AvailableData[d],'.csv',"'",')',sep='')))
  
  eval(parse(text=paste('Plot',AvailableData[d],'(',AvailableData[d],')',sep='')))
  
}

AssessmentResults<- as.data.frame(matrix(NA,nrow=length(Assessments)*10,ncol=9))

colnames(AssessmentResults)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')

AssessmentResults$Year<- as.numeric(AssessmentResults$Year)

Count<-0

Fish$LHITol<- 0.6


# Run Assessments ---------------------------------------------------------


for (f in 1:length(Fleets))
{
  
  FleetSpecies<- #Unique species caught by that fleet
  for (s in 1:length(FleetSpecies))
  {
    AvailableAssessments<- #Find assessments that you can run on that species in that fleet
      
      for (a in 1:length(AvailableAssessments))
      {
        
        #Pull out the data needed for this assessment
        
        #Run that assessment
        
        if (AvailableAssessments[a]=='LBAR') #Run LBAR assessment
        {
          
          Temp<- LBAR(LengthData,1,0.2,0,2007,NA,100,1,1,NA)$Output		
          # Temp2<- OldLBAR(LengthData,1,0.2,0,100,1,1)$Output
          
          DataLength<- dim(Temp)[1]
          
          AssessmentResults[(Count+1):(Count+DataLength),]<- Temp	
          
          Count<- Count+DataLength	
        }
        
        if (AvailableAssessments[a]=='CatchCurve') #Run Catch Curve analysis
        {
          
          Temp<- CatchCurve(LengthData,'AgeBased',1,2007,NA,1,100,1,1,1)$Output
          
          DataLength<- dim(Temp)[1]
          
          AssessmentResults[(Count+1):(Count+DataLength),]<- Temp	
          
          Count<- Count+DataLength	
        }
        
        
        if (AvailableAssessments[a]=='DensityRatio') #Run density ratio analysis 
        {
          Temp<- DensityRatio(DensityData,1,0.2,'Biomass',1000,1)$Output
          
          DataLength<- dim(Temp)[1]
          
          AssessmentResults[(Count+1):(Count+DataLength),]<- Temp	
          
          Count<- Count+DataLength
        }
        
        if (AvailableAssessments[a]=='CatchMSY')
        {
          Temp2<- CatchMSY(CatchData,3000,0.05,0,0,1,0,0,1,NA,c(0.75,0.99),NA,NA,c(0.25,0.65))
          
          Temp<- Temp2$Output
          
          DataLength<- dim(Temp)[1]
          
          AssessmentResults[(Count+1):(Count+DataLength),]<- Temp		
          
          Count<- Count+DataLength
          
        }
        
        if (AvailableAssessments[a]=='LBSPR') #Run LBSPR Assessment
        {
          Temp2<- LBSPR(LengthData,0,3,1,1,1)
          
          Temp<- Temp2$Output
          
          DataLength<- dim(Temp)[1]
          
          AssessmentResults[(Count+1):(Count+DataLength),]<- Temp		
          
          Count<- Count+DataLength
        }
        
        
        if (AvailableAssessments[a]=='DBSRA') #Run DBSRA Assessment
        {
          
          DCAC.start.yr <- CatchData$Year[1] #start of the catch period
          DCAC.end.yr<- CatchData$Year[length(CatchData$Year)] #end of the catch period
          delta.yr<- CatchData$Year[length(CatchData$Year)] #Year that current depletion is fit to
          DBSRA.OFL.yr<- CatchData$Year[length(CatchData$Year)] #Year to calculate DBSRA OFL outputs
          FMSYtoMratio <- 0.8 #ratio of Fmsy to M
          SD.FMSYtoMratio<- 0.05
          Delta<- 0.7
          SD.Delta<- 0.1
          DeltaLowerBound<- 0.5
          DeltaUpperBound<- 0.9
          BMSYtoB0ratio <- 0.3
          SD.BMSYtoB0ratio<- 0.1
          BMSYtoB0LowerBound<- 0.2
          BMSYtoB0UpperBound<- 0.5
          CatchInterp<-1
          NIter<- 500	
          
          
          Temp2<- DBSRA(CatchData, DCAC.start.yr, DCAC.end.yr, delta.yr, DBSRA.OFL.yr, FMSYtoMratio, SD.FMSYtoMratio, Delta, SD.Delta, DeltaLowerBound, DeltaUpperBound, BMSYtoB0ratio, SD.BMSYtoB0ratio, BMSYtoB0LowerBound, BMSYtoB0UpperBound, CatchInterp, NIter)
          
          Temp<- Temp2$Output
          
          DataLength<- dim(Temp)[1]
          
          AssessmentResults[(Count+1):(Count+DataLength),]<- Temp		
          
          Count<- Count+DataLength
        }
        
        
        show(paste('Finished ',Assessments[a],'-',round(100*a/length(Assessments),2),'% Done',sep=''))
        
        
        
        #store basic results, as well as list object results for each Fleet-Species-Combo
        
        
      } #close fleet-species-assessment assessment loop     
  } #close fleet-species assessment loop
} # close fleet assessment loop


AssessmentResults<- AssessmentResults[is.na(AssessmentResults$Year)==F,]
AssessmentResults$Year<- as.numeric(AssessmentResults$Year)
AssessmentResults$Value<- as.numeric(AssessmentResults$Value)
AssessmentResults$LowerCI<- as.numeric(AssessmentResults$LowerCI)
AssessmentResults$UpperCI<- as.numeric(AssessmentResults$UpperCI)
AssessmentResults$SD<- as.numeric(AssessmentResults$SD)

AssessmentResults[,4:7]<- round(AssessmentResults[,4:7],2)

show(AssessmentResults)
save.image(file=paste(ResultFolder,Fishery,'_Settings.RData',sep='')) #Save settings used to produce current results
write.csv(file=paste(ResultFolder,Fishery,'_Results.csv',sep=''),AssessmentResults) #Save current results