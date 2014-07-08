
# SNAP Operating Model Controlfile ----------------------------------------

#Description: This is the controlfile for the SNAP management strategy evaluation operating model. 
#Use this file to set basic operating parameters for running code locally on your machine. 


# Run Parameters ----------------------------------------------------------

RunName<- 'TimeStamp' #Batch name where run results will be stores. 'TimeStamp' creates a date stamped batch folder

Species<- 'All' #Which species you want to run.

Fleets<- 'All' #Which fleets you want to run

RunTime<- 30 #The number of time steps each model run will be evaluated for

#....


FontSize<- 14 #Master fontsize for figures
Font<- 'Helvetica' #Master font for figures

# Create Folders ----------------------------------------------------------

InputFolder<- 'Inputs/' #Folder where model data (.csv etc.) are stored 

if (RunName=='TimeStamp')
{
  SeedFolder<- paste('Created-',Sys.Date(),'/',sep="") #Folder to store outputs with a timestamp
}
if (RunName!='TimeStamp')
{
  SeedFolder<- paste(StoreRun,'/',sep='')
}
dir.create(SeedFolder)

FigureFolder<- paste(SeedFolder,'Figures/',sep='')
ResultFolder<- paste(SeedFolder,'Data/',sep='')

dir.create(FigureFolder)
dir.create(ResultFolder)


pdf(file=paste(FigureFolder,'Test.pdf'),family=Font,pointsize=FontSize)
plot(1:10,pch=18)
dev.off()


