#####################################
#Sampling Function for SNAP working group
# 
#Collects and Aggregates (optional) 10 different kinds of data. 
#######################################################################


CollectData = function(timeStep, simTime, burn, FleetN, DataParams, NoTakeZone, SpaceNumAtAgeT,
                       SpaceCatAgeByFisher, SpaceEffort, wgtAtAge, lenAtAge,lenSD, PortX, PortY,kmax){
  #Create Storage Structures Based on Number of Fleets and Aggregation
  if(DataParams$Aggregate == 0){
    histCatchDatOUTPUT <-  array(NA,dim=c(length((DataParams$histStartYr):(timeStep-burn)), dim(NoTakeZone)[1], dim(NoTakeZone)[2],FleetN))
    histCPUEDatOUTPUT <- array(NA,dim=c(length((DataParams$histStartYr):(timeStep-burn)), dim(NoTakeZone)[1], dim(NoTakeZone)[2],FleetN))
    FDCatchOUTPUT <- array(NA,dim=c(dim(NoTakeZone)[1], dim(NoTakeZone)[2], FleetN))
    FDCPUEOUTPUT <- array(NA,dim=c(dim(NoTakeZone)[1], dim(NoTakeZone)[2], FleetN))
    AgeFDMatOUTPUT <- array(NA,dim=c(DataParams$nAgeFD,dim(NoTakeZone)[1], dim(NoTakeZone)[2], FleetN))
    SizeFDMatOUTPUT <- array(NA,dim=c(DataParams$nSizeFD,dim(NoTakeZone)[1], dim(NoTakeZone)[2], FleetN))
  } else {
    histCatchDatOUTPUT <- array(NA,dim=c(length((DataParams$histStartYr):(timeStep-burn)), FleetN))
    histCPUEDatOUTPUT <- array(NA,dim=c(length((DataParams$histStartYr):(timeStep-burn)), FleetN))
    FDCatchOUTPUT <- rep(NA,length = FleetN)
    FDCPUEOUTPUT <- rep(NA,length = FleetN)
    AgeFDMatOUTPUT <- array(NA,dim=c(DataParams$nAgeFD*dim(NoTakeZone)[1]*dim(NoTakeZone)[2], FleetN))
    SizeFDMatOUTPUT <- array(NA,dim=c(DataParams$nSizeFD*dim(NoTakeZone)[1]*dim(NoTakeZone)[2], FleetN))
  }

  ### Conduct FD sampling for each fleet
  for(fl in 1:FleetN){
    #Determine what and where to sample
    if(DataParams$FDSurvType == 0){FDSurvMat <- read.csv("FDSurvMat.csv") }
    if(DataParams$FDSurvType == 1){FDSurvMat <- matrix(1,nrow=dim(NoTakeZone)[1],ncol=dim(NoTakeZone)[2]) }
    if(DataParams$FDSurvType == 2){FDSurvMat <- matrix(0,nrow=dim(NoTakeZone)[1],ncol=dim(NoTakeZone)[2]) }
    if(DataParams$FDSurvType == 3){ 
      TotPatches <- dim(NoTakeZone)[1]*dim(NoTakeZone)[2]
      numberSampled <- round(TotPatches*0.2)
      randPatchNums <- sample(1:TotPatches, numberSampled, replace=FALSE)
      FDSurvMat <- matrix(0,nrow=dim(NoTakeZone)[1],ncol=dim(NoTakeZone)[2])
      FDSurvMat <- replace(FDSurvMat, randPatchNums, 1)
    }
    if(DataParams$FDSurvType == 4){
      # Sample 10% of patches close to port, 10% of patches far from port
      TotPatches <- dim(NoTakeZone)[1]*dim(NoTakeZone)[2]
      numberSampled <- round(TotPatches*0.1)
      #Calculate distance from port
      Distance<-sqrt((row(SpaceUse)-PortX[fl])^2 + (col(SpaceUse)-PortY[fl])^2)
      nearSites <- order(Distance)[1:numberSampled]
      farSites <- order(Distance)[(TotPatches-numberSampled):TotPatches]
      FDSurvMat <- matrix(0,nrow=dim(NoTakeZone)[1],ncol=dim(NoTakeZone)[2]) 
      FDSurvMat <- replace(FDSurvMat, c(nearSites,farSites), 1)
    }
    
    collCatchFD <- FDSurvMat * DataParams$Catch_FD
    collEffortFD <- FDSurvMat * DataParams$Effort_FD
    collAgeFD <- FDSurvMat * DataParams$Ages_FD
    collSizeFD <- FDSurvMat * DataParams$Sizes_FD
    
    ####Collect historical data for current year
    histCatchDat <- CollectHistCatchData((timeStep-burn),burn,DataParams$SampStartYr,SpaceCatAgeByFisher,DataParams$histCatchFD,DataParams$histStartYr,DataParams$Aggregate, DataParams$sigHistCatch, DataParams$catchBias, fl,wgtAtAge)
    histCPUEDat <- CollectHistCPUEData(SpaceEffort,DataParams$histEffortFD,DataParams$histStartYr,DataParams$Aggregate,(timeStep-burn),burn,DataParams$SampStartYr,SpaceCatAgeByFisher,DataParams$sigHistCatch,DataParams$catchBias,fl,wgtAtAge)
    
    #### Collect FD data from these patches
    FDCatch <- CollectFisheryCatch((timeStep-burn+1),SpaceCatAgeByFisher,DataParams$Aggregate,fl,wgtAtAge,DataParams$sigHistCatch,collCatchFD)
    FDCPUE <- CollectFisheryCPUE((timeStep-burn+1),SpaceCatAgeByFisher,SpaceEffort,fl,collEffortFD,DataParams$Aggregate,wgtAtAge,DataParams$sigHistCatch)
    AgeFDMat <- CollectFisheryAges((timeStep-burn+1),SpaceCatAgeByFisher,DataParams$nAgeFD,collAgeFD,kmax,fl,DataParams$Aggregate)
    SizeFDMat <- CollectFisherySizes((timeStep-burn+1),SpaceCatAgeByFisher,DataParams$nSizeFD,collSizeFD,lenSD,lenAtAge,DataParams$Aggregate,fl)
    
    #### Assign to Ouput Arrays
    if(DataParams$Aggregate == 0){
      histCatchDatOUTPUT[,,,fl] <- histCatchDat
      histCPUEDatOUTPUT[,,,fl] <- histCPUEDat
      FDCatchOUTPUT[,,fl] <- FDCatch
      FDCPUEOUTPUT[,,fl] <- FDCPUE
      AgeFDMatOUTPUT[,,,fl] <- AgeFDMat
      SizeFDMatOUTPUT[,,,fl] <- SizeFDMat
    } else {
      histCatchDatOUTPUT[,fl] <- histCatchDat
      histCPUEDatOUTPUT[,fl] <- histCPUEDat
      FDCatchOUTPUT[fl] <- FDCatch
      FDCPUEOUTPUT[fl] <- FDCPUE
      AgeFDMatOUTPUT[,fl] <- AgeFDMat
      SizeFDMatOUTPUT[,fl] <- SizeFDMat
    }
  } # for FleetN
  
  ####Return list of collected data
  DataOut<-NULL
  DataOut$histCatchDat <- histCatchDatOUTPUT
  DataOut$histCPUEDat <- histCPUEDatOUTPUT
  DataOut$FDCatch <- FDCatchOUTPUT
  DataOut$FDCPUE<-FDCPUEOUTPUT 
  DataOut$AgeFD <- AgeFDMatOUTPUT
  DataOut$SizeFD <- SizeFDMatOUTPUT   
#   DataOut$CatchFI <- CatchFI
#   DataOut$CPUEFI <- CPUEFI
#   DataOut$numFI <- numFI
#   DataOut$AgeFI <- AgeFI
#   DataOut$SizeFI <- SizeFI
#   DataOut$VBDat <- VBDat
  
  ### output Single object with 
  
  return(DataOut)
}

###################################################
# Below this is FI functuions to be updated soon.
####################################################

#     ####Survey population
#     Survey <- GoSurvey(timeStep,nyrs,fishPop,survq,survEffort,SurvPatch,Sel50FI,Sel95FI,FFleet,minLen,maxLen,LengthBins)
#     
#   #### Identify patches that have a catch to sample
#   FDsurvPatch<-rep(0,length(collAgeFI))
#   FDsurvPatch[colSums(as.matrix(FFleet$CatchatAge[timeStep,,])) > 0] <- 1
#   
#   #### Identify patches that have a survey to sample
#   FIsurvPatch <- rep(0,length(collAgeFI))
#   FIsurvPatch[colSums(as.matrix(Survey$CatchatAge_FI)) > 0] <- 1
#   
#   #### Collect FI data from these patches
#   CatchFI <- CollectCatchFI(timeStep,nyrs,Survey$CatchatSize_FI,collCatchFI,FFleet,Aggregate,minLen, maxLen, LengthBins)
#   CPUEFI <- CollectSurvCPUE(Survey$CatchatSize_FI,survEffort,SurvPatch,collCatchFI,FFleet,Aggregate,minLen, maxLen, LengthBins)
#   numFI <- CollectSurvnum(Survey$CatchatSize_FI,survEffort,SurvPatch,collCatchFI,FFleet,Aggregate)
#   AgeFIMat <- CollectSurvAges(timeStep,collAgeFI,fishPop,nAgeFI,Survey,FIsurvPatch)
#   SizeFIMat <- CollectSurvSize(collSizeFI,FIsurvPatch,minLen,mxLen,Sel50FI,Sel95FI,fishPop,LengthBins,Survey,nSizeFI)
#   VBDat <- CollectVBData(timeStep,VBbins,collSizeFI,fishPop,survq,SurvPatch,Sel50VB,Sel95VB,minLen,maxLen,FIsurvPatch,numVB)
#   


##############################################################################
# Historical Fishery Dependent Data Collection Functions
##############################################################################
                        
CollectHistCatchData = function(time,burn,SampStartYr,SpaceCatAgeByFisher,histCatchFD,histStartYr,Aggregate, sigHistCatch, catchBias, fl,wgtAtAge){  
  if (time == (SampStartYr) & histCatchFD == 1){
    #Catch at age is converted to catch at weight.
    CatchAtWeight<-sweep(SpaceCatAgeByFisher[,,,,,fl],MARGIN=4,wgtAtAge,`*`)
    
    # True catches at age is sampled -- arrays- indices you want to keep?
    SCATByFleet <- apply(CatchAtWeight[,,,,],1:4,sum,na.rm=TRUE)  #sums over fisher
    SCTByFleet <- apply(SCATByFleet[(histStartYr):(time),,,],1:3,sum,na.rm=TRUE)  #sums over ages
    SCTByFleetNew <- replace(SCTByFleet[,,],which(SCTByFleet[,,]==0),NA)
    
    #Calculate Historical observation error in catch: If not reporting error, sigHistCatch = 0
    #obsError <- rnorm(n=length(histStartYr:time),mean = -(sigHistCatch)^2/2,sd = sigHistCatch) ##### Error in reported catches is based on Bousquet et al. 2009
    obsError <- rnorm(n=length((histStartYr):(time)),mean = 0,sd = sigHistCatch)
    obsevervedCatches <- sweep(SCTByFleetNew,MARGIN=1,obsError,'+')
    
    # Add in a potential reporting bias
    biasedCatches <- obsevervedCatches + obsevervedCatches*catchBias #catchBias is (-1,1),with 0 = no Bias
    CatchHist_FD <-  biasedCatches
    
    # Aggregate Historical Catches?
    if(Aggregate == 0){
      Final <- CatchHist_FD
    } else {   #sum across all patches for all years
      Final <- apply(CatchHist_FD,1,sum,na.rm=TRUE)
    } 
  } else {
    Final <- NA
  }
  return(Final) 
}

CollectHistCPUEData = function(SpaceEffort,histEffortFD,histStartYr,Aggregate,time,burn,SampStartYr,SpaceCatAgeByFisher,sigHistCatch,catchBias,fl,wgtAtAge){  
  if (time == SampStartYr & histEffortFD == 1){
    #Catch at age is converted to catch at weight.
    CatchAtWeight<-sweep(SpaceCatAgeByFisher[,,,,,fl],MARGIN=4,wgtAtAge,`*`)
    
    # True catches at age is sampled -- arrays- indices you want to keep?
    SCATByFleet <- apply(CatchAtWeight[,,,,],1:4,sum,na.rm=TRUE)  #sums over fisher
    SCTByFleet <- apply(SCATByFleet[(histStartYr):(time),,,],1:3,sum,na.rm=TRUE)  #sums over ages
    SCTByFleetNew <- replace(SCTByFleet[,,],which(SCTByFleet[,,]==0),NA)
    
    #Calculate Historical observation error in catch: If not reporting error, sigHistCatch = 0
    #obsError <- rnorm(n=length(histStartYr:time),mean = -(sigHistCatch)^2/2,sd = sigHistCatch) ##### Error in reported catches is based on Bousquet et al. 2009
    obsError <- rnorm(n=length((histStartYr):(time)),mean = 0,sd = sigHistCatch)
    obsevervedCatches <- sweep(SCTByFleetNew,MARGIN=1,obsError,'+')
    
    # Add in a potential reporting bias
    biasedCatches <- obsevervedCatches + obsevervedCatches*catchBias #catchBias is (-1,1),with 0 = no Bias
    CatchHist_FD <-  biasedCatches
    
    # Eliminate Zeros in Effort to avoid dividing by zero
    HistSpaceEffort <- SpaceEffort[(histStartYr):(time),,,fl]
    NewSpaceEff <- replace(HistSpaceEffort[,,],which(HistSpaceEffort[,,]==0),NA)
    
    # Divide effort per patch
    if(Aggregate == 0){
      Final <- CatchHist_FD / NewSpaceEff
    } else {   #sum across all patches for all years
      AggCPUE <- apply(NewSpaceEff,1,sum,na.rm=TRUE)
      AggCPUENew <- replace(AggCPUE, which(AggCPUE == 0), NA)
      Final <- apply(CatchHist_FD,1,sum,na.rm=TRUE)/AggCPUENew 
    }
  } else {
    Final <- NA
  }
  return(Final)
}

##############################################################################
# Current Fishery Dependent Data Collection Functions
##############################################################################
CollectFisheryCatch <- function(time,SpaceCatAgeByFisher,Aggregate,fl,wgtAtAge,sigHistCatch,collCatchFD){
  #Catch at age is converted to catch at weight.
  CatchAtWeight<-sweep(SpaceCatAgeByFisher[time,,,,,fl],MARGIN=3,wgtAtAge,`*`)
  
  # True catches at age is sampled -- arrays- indices you want to keep?
  SCATByFleet <- apply(CatchAtWeight[,,,],1:3,sum,na.rm=TRUE)  #sums over fisher
  SCTByFleet <- apply(SCATByFleet[,,],1:2,sum,na.rm=TRUE)  #sums over ages
  SCTByFleetNew <- replace(SCTByFleet[,],which(SCTByFleet[,]==0),NA)
  
  #Calculate Historical observation error in catch: If not reporting error, sigHistCatch = 0
  #obsError <- rnorm(n=length(histStartYr:time),mean = -(sigHistCatch)^2/2,sd = sigHistCatch) ##### Error in reported catches is based on Bousquet et al. 2009
  obsError <- rnorm(n=1, mean = 0,sd = sigHistCatch)
  obsevervedCatches <- obsError + SCTByFleetNew
  FisheryCatch <- obsevervedCatches 
  
  if (Aggregate == 0){
    Final <- FisheryCatch
  } else {
    Final <- sum(FisheryCatch,na.rm=TRUE)
  }
  return(Final)
}

CollectFisheryCPUE <- function(time,SpaceCatAgeByFisher,SpaceEffort,fl,collEffortFD,Aggregate,wgtAtAge, sigHistCatch){
  
  #Catch at age is converted to catch at weight.
  CatchAtWeight<-sweep(SpaceCatAgeByFisher[time,,,,,fl],MARGIN=3,wgtAtAge,`*`)
  
  # True catches at age is sampled -- arrays- indices you want to keep?
  SCATByFleet <- apply(CatchAtWeight[,,,],1:3,sum,na.rm=TRUE)  #sums over fisher
  SCTByFleet <- apply(SCATByFleet[,,],1:2,sum,na.rm=TRUE)  #sums over ages
  SCTByFleetNew <- replace(SCTByFleet[,],which(SCTByFleet[,]==0),NA)
  
  #Calculate Historical observation error in catch: If not reporting error, sigHistCatch = 0
  #obsError <- rnorm(n=length(histStartYr:time),mean = -(sigHistCatch)^2/2,sd = sigHistCatch) ##### Error in reported catches is based on Bousquet et al. 2009
  obsError <- rnorm(n=1, mean = 0,sd = sigHistCatch)
  obsevervedCatches <- obsError + SCTByFleetNew
  FisheryCatch <- obsevervedCatches 
  
  ## effort data collected from each patch?
  Eff_Temp <-SpaceEffort[time,,,fl]*collEffortFD ## should end up with effort in eahc patch, with 0s where there is no data
  Effort_FD <- replace(Eff_Temp, Eff_Temp == 0, NA)  #avoiding divide by zero problem in unsampled patches
  if (Aggregate == FALSE){
    Final <- FisheryCatch/Effort_FD
  } else {
    Final <- sum(FisheryCatch,na.rm=TRUE)/sum(Effort_FD,na.rm=TRUE)
  }
  return(Final)
}

#Collect Age Data
CollectFisheryAges <- function(time,SpaceCatAgeByFisher,nAgeFD,collAgeFD, kmax, fl, Aggregate){
  #Sum over fishers for this fleet and this year
  SCAT <- apply(SpaceCatAgeByFisher[time,,,,,fl],1:3,sum,na.rm=TRUE)  #sums over fisher
  AgeDatHolder <- array(NA,dim=c(max(nAgeFD,1),dim(collAgeFD)[1],dim(collAgeFD)[2]))
  
  #Figure out if there was any catch in each patch
  CatchPerPatch <- apply(SpaceCatAgeByFisher[time,,,,,fl],1:2,sum,na.rm=TRUE)  #sums over fisher
  
  #Loop through possible patches
  for(rows in 1:dim(collAgeFD)[1]){
    for(cols in 1:dim(collAgeFD)[1]){
      if (collAgeFD[rows,cols] == 1 & CatchPerPatch[rows,cols] > 0){ 
        set.seed(fl*rows*cols*time)
        AgeDatHolder[,rows,cols]<-sample(1:kmax,nAgeFD,replace=TRUE,prob=SCAT[rows,cols,])
      }
    }
  }
  
  # Aggregate if necessary
  if (Aggregate == FALSE){
    AgeDatHolderFINAL <- AgeDatHolder
  } else {
    AgeDatHolderFINAL <- as.vector(AgeDatHolder)
  }
  return(AgeDatHolderFINAL)
}

#Collect FD Size Data 
CollectFisherySizes <- function(time,SpaceCatAgeByFisher,nSizeFD,collSizeFD,lenSD,lenAtAge,Aggregate,fl){
  ## Collect Ages using the same procedure as above
  #Sum over fishers for this fleet and this year
  SCAT <- apply(SpaceCatAgeByFisher[time,,,,,fl],1:3,sum,na.rm=TRUE)  #sums over fisher
  AgeDatHolder <- array(NA,dim=c(max(nSizeFD,1),dim(collSizeFD)[1],dim(collSizeFD)[2]))
  
  #Figure out if there was any catch in each patch
  CatchPerPatch <- apply(SpaceCatAgeByFisher[time,,,,,fl],1:2,sum,na.rm=TRUE)  #sums over fisher
  
  #Loop through possible patches
  for(rows in 1:dim(collSizeFD)[1]){
    for(cols in 1:dim(collSizeFD)[2]){
      if (collSizeFD[rows,cols] == 1 & CatchPerPatch[rows,cols] > 0){ 
        set.seed(fl*rows*cols*time)
        AgeDatHolder[,rows,cols]<-sample(1:kmax,nSizeFD,replace=TRUE,prob=SCAT[rows,cols,])
      }
    }
  }
  
  #Storage for size data
  SizeDatHolder <- array(NA,dim=c(max(nSizeFD,1),dim(collSizeFD)[1],dim(collSizeFD)[2]))
  ## Assign sizes to each age
  
  LengthVec = seq(0, ceiling(Linf)*1.4, by= 1)   #changed to make 1cm length bins
  LenMids = seq(LengthVec[1] +((LengthVec[2]-LengthVec[1])/2), by=(LengthVec[2]-LengthVec[1]), length=length(LengthVec)-1)
  AgeLenProbMat = matrix(0, nrow=length(LengthVec)-1, ncol=kmax) 
  for (i in seq_along(1:kmax)) 
  {
    for (lengthclass in seq_along(LengthVec)[-1])
    {  
      length_cat = lengthclass-1   #Calculate the probability of falling into each length bin
      if(length_cat==1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], lenAtAge[i],  lenSD) 
      if(length_cat>1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], lenAtAge[i],  lenSD) - pnorm(LengthVec[lengthclass-1], lenAtAge[i], lenSD)       
      
    }
    # For each Age normalize across lengths
    AgeLenProbMat[,i] = AgeLenProbMat[,i]/sum(AgeLenProbMat[,i])
  }
  
  #Loop through possible patches
  for(rows in 1:dim(collSizeFD)[1]){
    for(cols in 1:dim(collSizeFD)[2]){
      if (collSizeFD[rows,cols] == 1 & CatchPerPatch[rows,cols] > 0){ 
        for(a in 1:nSizeFD){
          SizeDatHolder[a,rows,cols]<-sample(LenMids,1,replace=TRUE,prob=AgeLenProbMat[,AgeDatHolder[a,rows,cols]])
        }
      }
    }
  }  
  
  # Aggregate if necessary
  if (Aggregate == FALSE){
    SizeDatHolderFINAL <- SizeDatHolder
  } else {
    SizeDatHolderFINAL <- as.vector(SizeDatHolder)
  }
  return(SizeDatHolderFINAL)
}
##############################################################################
# Survey/ Fishery Independent Data Collection Functions
##############################################################################

# #Determine Catch at age from FI survey
# GoSurvey = function(timeStep,nyrs,fishPop,survq,survEffort,SurvPatch,sel50,sel95,FFleet,minLen,maxLen,LengthBins){
#   #Calculate survey selectivity at age
#   selFI<-ConvertSel(sel50,sel95,fishPop$fish$lengthSD,fishPop)
#   Fsurv <- survq*((survEffort/length(SurvPatch))*SurvPatch)
#   Harv <- matrix((1-exp(-Fsurv)),nrow = dim(selFI)[1],ncol = dim(selFI)[2],byrow=T)
#   CatchatAge_FI <- pmin(fishPop$N[timeStep,,], fishPop$N[timeStep,,]*(selFI*Harv))
#   
#   #Calculate the probability of falling into length bins
#   LengthVec = seq(minLen, maxLen, by= LengthBins)   #changed to make 1cm length bins
#   LenMids = seq(LengthVec[1] +((LengthVec[2]-LengthVec[1])/2), by=(LengthVec[2]-LengthVec[1]), length=length(LengthVec)-1)
#   AgeLenProbMat = matrix(0, nrow=length(LengthVec)-1, ncol=fishPop$fish$maxAge) 
#   for (i in seq_along(1:fishPop$fish$maxAge)) 
#   {
#     for (lengthclass in seq_along(LengthVec)[-1])
#     {  
#       length_cat = lengthclass-1
#       if(length_cat==1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], fishPop$fish$length[i], 
#                                                             fishPop$fish$lengthSD[i]) 
#       if(length_cat>1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], fishPop$fish$length[i], 
#                                                            fishPop$fish$lengthSD[i]) - pnorm(LengthVec[lengthclass-1], fishPop$fish$length[i], fishPop$fish$lengthSD[i])       
#     }
#     AgeLenProbMat[,i] = AgeLenProbMat[,i]/sum(AgeLenProbMat[,i])
#   }
#   
#   ##Calculate the selectivity at length for each length bin
#   SelLengthVec = plogis(LenMids, location=sel50, scale=(sel95-sel50)/log(19))
#   #SelLengthVec = rep(1,56)
#   
#   ###Create matrix of normalized probabilities of being above the size and selected.
#   Prob2Mat = AgeLenProbMat* SelLengthVec    #SelLengthMat
#   NormProbMat = Prob2Mat
#   ColSums = apply(Prob2Mat, 2, sum)
#   for (i in 1:fishPop$fish$maxAge) {
#     NormProbMat[,i] = Prob2Mat[,i]/ColSums[i] # Normalise probabilities
#   }
#   
#   #Calculate the Catch at Size
#   CatchatSize_FI <- matrix(NA,nrow=length(LenMids),ncol=fishPop$nPatch)
#   for(i in 1:fishPop$nPatch){
#     CatchLengthMat <-  CatchatAge_FI[,i] * t(NormProbMat)
#     CatchatSize_FI[,i] <- apply(CatchLengthMat, 2, sum)
#   }
#   
#   Outs<-NULL
#   Outs$CatchatAge_FI<-CatchatAge_FI
#   Outs$CatchatSize_FI <- CatchatSize_FI
#   return(Outs)
# }
# 
# CollectCatchFI <- function(timeStep,firstDatYr,CatchatSizeFI,collCatchFI,FFleet,Aggregate,minLen, maxLen, LengthBins){
#   LengthVec = seq(minLen, maxLen, by= LengthBins)   #changed to make 1cm length bins
#   LenMids = seq(LengthVec[1] +((LengthVec[2]-LengthVec[1])/2), by=(LengthVec[2]-LengthVec[1]), length=length(LengthVec)-1)
#   WeightAtLength = FFleet$fishPop$fish$wa*LenMids^FFleet$fishPop$fish$wb
#   #weight <- FFleet$fishPop$fish$weightAtAge()
#   if (length(collCatchFI) > 1){
#     Catch_FI <- colSums(as.matrix(CatchatSizeFI[1:(length(LengthVec)-1),]*WeightAtLength))
#   }else{
#     Catch_FI <- sum(CatchatSizeFI[1:(length(LengthVec)-1)]*WeightAtLength)
#   }
#   if (Aggregate == FALSE){
#     Final <- Catch_FI
#   } else {
#     Final <- sum(Catch_FI,na.rm=TRUE)
#   }
#   return(Final)
# }
# 
# #Determine survey CPUE
# CollectSurvCPUE <- function(CatchatSizeFI,survEffort,SurvPatch,collCatchFI,FFleet,Aggregate,minLen, maxLen, LengthBins){
#   LengthVec = seq(minLen, maxLen, by= LengthBins)   #changed to make 1cm length bins
#   LenMids = seq(LengthVec[1] +((LengthVec[2]-LengthVec[1])/2), by=(LengthVec[2]-LengthVec[1]), length=length(LengthVec)-1)
#   WeightAtLength = FFleet$fishPop$fish$wa*LenMids^FFleet$fishPop$fish$wb
#   #SurvEffTemp <- replace(((survEffort/length(SurvPatch))*SurvPatch), survEffort*SurvPatch == 0, NA) #avoiding dividing by zero
#   SurvEffTemp <- replace(((survEffort/sum(SurvPatch))*SurvPatch), survEffort*SurvPatch == 0, NA) 
#   if (length(collCatchFI) >1){
#     SurvCatch <- colSums(as.matrix(CatchatSizeFI[1:(length(LengthVec)-1),]*WeightAtLength))
#   }else{
#     SurvCatch <- sum(CatchatSizeFI[1:(length(LengthVec)-1)]*WeightAtLength) 
#   }
#   CPUE_FI <- SurvCatch/SurvEffTemp
#   if (Aggregate == FALSE){
#     Final <- CPUE_FI
#   } else {
#     Final <- sum(SurvCatch,na.rm=TRUE)/sum(SurvEffTemp,na.rm=TRUE)
#   }
#   return(Final)
# } 
# 
# #collect CPUE in numbers by size -- needed in dtree
# CollectSurvnum <- function(CatchatSizeFI,survEffort,SurvPatch,collCatchFI,FFleet,Aggregate){
#   SurvEffTemp <- replace(((survEffort/sum(SurvPatch))*SurvPatch), survEffort*SurvPatch == 0, NA) #avoiding dividing by zero, and effort is only allocated to patches that are surveyed
#   # SurvEffTemp <- replace(((survEffort/length(SurvPatch))*SurvPatch), survEffort*SurvPatch == 0, NA) #avoiding dividing by zero
#   if (length(collCatchFI) >1){
#     #SurvCatchNums <- as.matrix(CatchatSizeFI[1:FFleet$NumSizes,which(collCatchFI==1)])
#     SurvCatchNums <- as.matrix(CatchatSizeFI[1:FFleet$NumSizes,])
#   }else{
#     #SurvCatchNums <- CatchatSizeFI[1:FFleet$NumSizes,]*collCatchFI  
#     SurvCatchNums <- CatchatSizeFI[1:FFleet$NumSizes,]
#   }
#   numCPUE_FI <- round(SurvCatchNums,0)   ##Converts to whole numbers   
#   
#   SurvEffMat <- matrix(SurvEffTemp,nrow=dim(numCPUE_FI)[1],ncol=dim(numCPUE_FI)[2],byrow=T)
#   if (Aggregate == FALSE){
#     Final <- numCPUE_FI/SurvEffMat
#   } else {
#     Final <- rowSums(numCPUE_FI,na.rm=TRUE)/sum(SurvEffTemp,na.rm=TRUE)
#   }
#   return(Final)
# }
# 
# #Collect Age Data -- assign ages to size data
# CollectSurvAges <- function(timeStep,collAgeFI,fishPop,nAgeFI,Survey,FIsurvPatch){
#   CatchatAge_FI <-  as.matrix(round(Survey$CatchatAge))  #turns the catch at size into whole numbers
#   AgeDatHolder <- matrix(NA,nrow=max(nAgeFI,1),ncol=length(collAgeFI))
#   for(patch in 1:length(collAgeFI)){
#     if ((FIsurvPatch * collAgeFI)[patch] == 1){ 
#       AgeDatHolder[,patch] <- sample(1:fishPop$fish$maxAge,nAgeFI,replace=TRUE,prob=(CatchatAge_FI[,patch]/sum(CatchatAge_FI[,patch])))
#       
#     }
#   }
#   return(AgeDatHolder)
# }
# 
# #Collect VB data -- age and size are paired. 
# CollectVBData <- function(timeStep,VBbins,collSizeFI,fishPop,survq,SurvPatch,Sel50VB,Sel95VB,minLen,maxLen,FIsurvPatch,numVB){
#   #Calculate survey selectivity at age
#   survEffort <- 25
#   selFI<-ConvertSel(Sel50VB,Sel95VB,fishPop$fish$lengthSD,fishPop)
#   Fsurv <- survq*((survEffort/length(SurvPatch))*SurvPatch)
#   Harv <- matrix((1-exp(-Fsurv)),nrow = dim(selFI)[1],ncol = dim(selFI)[2],byrow=T)
#   CatchatAge_FI <- round(pmin(fishPop$N[timeStep,,], fishPop$N[timeStep,,]*(selFI*Harv)))
#   
#   MatSize <- max(colSums(CatchatAge_FI))
#   AgeDatHolder <- matrix(NA,nrow=MatSize,ncol=length(collSizeFI))
#   for(patch in 1:length(collSizeFI)){
#     if ((FIsurvPatch * collSizeFI)[patch] == 1){ 
#       vec <- NULL
#       for(a in 1:fishPop$fish$maxAge){
#         vec <- append(vec,rep(a,CatchatAge_FI[a,patch]))
#       }
#       AgeDatHolder[1:length(vec),patch]<-vec
#     }
#   }
#   
#   # Create age-length probability matrix
#   LengthVec = seq(minLen, maxLen, by= VBbins)   #changed to make 1cm length bins
#   LenMids = seq(LengthVec[1] +((LengthVec[2]-LengthVec[1])/2), by=(LengthVec[2]-LengthVec[1]), length=length(LengthVec)-1)
#   AgeLenProbMat = matrix(0, nrow=length(LengthVec)-1, ncol=fishPop$fish$maxAge) 
#   for (i in seq_along(1:fishPop$fish$maxAge)) 
#   {
#     for (lengthclass in seq_along(LengthVec)[-1])
#     {  
#       length_cat = lengthclass-1   #Calculate the probability of falling into each length bin
#       if(length_cat==1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], fishPop$fish$length[i],  fishPop$fish$lengthSD[i]) 
#       if(length_cat>1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], fishPop$fish$length[i],  fishPop$fish$lengthSD[i]) - pnorm(LengthVec[lengthclass-1], fishPop$fish$length[i], fishPop$fish$lengthSD[i])       
#     }
#     # For each Age normalize across lengths
#     AgeLenProbMat[,i] = AgeLenProbMat[,i]/sum(AgeLenProbMat[,i])
#   }
#   
#   ##Calculate the selectivity at length for each length bin
#   SelLengthVec = plogis(LenMids, location=Sel50VB, scale=(Sel95VB-Sel50VB)/log(19))
#   #SelLengthVec = rep(1,34)
#   
#   ###Create matrix of normalized probabilities of being above the size and selected.
#   Prob2Mat = AgeLenProbMat* SelLengthVec    #SelLengthMat
#   
#   ## Need to re-normalized across all lengths after selectivity was accounted for 
#   NormProbMat = Prob2Mat
#   ColSums = apply(Prob2Mat, 2, sum)
#   for (i in 1:fishPop$fish$maxAge) {
#     NormProbMat[,i] = Prob2Mat[,i]/ColSums[i] # Normalise probabilities
#   }
#   
#   SizeDatHolder<-  matrix(NA,nrow=dim(AgeDatHolder)[1],ncol=dim(AgeDatHolder)[2])
#   for(patch in 1:length(collSizeFI)){
#     if ((FIsurvPatch * collSizeFI)[patch] == 1){
#       TempSize <- NULL
#       for(a in 1:fishPop$fish$maxAge){
#         Na <- length(AgeDatHolder[which(AgeDatHolder[,patch]==a),patch])
#         TempSize <- append(TempSize,sample(LenMids, Na, replace=TRUE, prob=NormProbMat[,a]))
#       }
#     } else {
#       TempSize <- rep(NA,dim(AgeDatHolder)[1])
#     }
#     SizeDatHolder[1:length(TempSize),patch] <- TempSize
#   }
#   
#   # Sample VB data
#   SampInd <- sample(1:MatSize,NumVB,replace=TRUE)
#   AgeDatHolderVB <- matrix(NA,nrow=MatSize,ncol=length(collSizeFI))
#   SizeDatHolderVB <- matrix(NA,nrow=MatSize,ncol=length(collSizeFI))
#   AgeDatHolderVB <- AgeDatHolder[SampInd,]
#   SizeDatHolderVB <- SizeDatHolder[SampInd,]
#   
#   Outs<-NULL
#   Outs$Age <- AgeDatHolderVB
#   Outs$Size <- SizeDatHolderVB
#   return(Outs)
# }
# 
# #Collect FI Size Data
# CollectSurvSize <- function(collSizeFI,FIsurvPatch,minLen,mxLen,Sel50FI,Sel95FI,fishPop,LengthBins,Survey,nSizeFI){
#   CatchatAge_FI <-  as.matrix(round(Survey$CatchatAge))  #turns the catch at size into whole numbers
#   AgeDatHolder <- matrix(NA,nrow=max(nSizeFI,1),ncol=length(collSizeFI))
#   for(patch in 1:length(collSizeFI)){
#     if ((FIsurvPatch * collSizeFI)[patch] == 1){ 
#       AgeDatHolder[,patch] <- sample(1:fishPop$fish$maxAge,nSizeFI,replace=TRUE,prob=(CatchatAge_FI[,patch]/sum(CatchatAge_FI[,patch])))
#       
#     }
#   }
#   
#   # Create age-length probability matrix
#   LengthVec = seq(minLen, maxLen, by= LengthBins)   #changed to make 1cm length bins
#   LenMids = seq(LengthVec[1] +((LengthVec[2]-LengthVec[1])/2), by=(LengthVec[2]-LengthVec[1]), length=length(LengthVec)-1)
#   AgeLenProbMat = matrix(0, nrow=length(LengthVec)-1, ncol=fishPop$fish$maxAge) 
#   for (i in seq_along(1:fishPop$fish$maxAge)) 
#   {
#     for (lengthclass in seq_along(LengthVec)[-1])
#     {  
#       length_cat = lengthclass-1   #Calculate the probability of falling into each length bin
#       #if(length_cat==1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], fishPop$fish$length[i],  10) 
#       if(length_cat==1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], fishPop$fish$length[i],  fishPop$fish$lengthSD[i]) 
#       #if(length_cat>1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], fishPop$fish$length[i],  10) - pnorm(LengthVec[lengthclass-1], fishPop$fish$length[i], 12)       
#       if(length_cat>1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], fishPop$fish$length[i],  fishPop$fish$lengthSD[i]) - pnorm(LengthVec[lengthclass-1], fishPop$fish$length[i], fishPop$fish$lengthSD[i])       
#     }
#     # For each Age normalize across lengths
#     AgeLenProbMat[,i] = AgeLenProbMat[,i]/sum(AgeLenProbMat[,i])
#   }
#   
#   ##Calculate the selectivity at length for each length bin
#   SelLengthVec = plogis(LenMids, location=Sel50FI, scale=(Sel95FI-Sel50FI)/log(19))
#   #SelLengthVec = rep(1,34)
#   
#   ###Create matrix of normalized probabilities of being above the size and selected.
#   Prob2Mat = AgeLenProbMat* SelLengthVec    #SelLengthMat
#   
#   ## Need to re-normalized across all lengths after selectivity was accounted for 
#   NormProbMat = Prob2Mat
#   ColSums = apply(Prob2Mat, 2, sum)
#   for (i in 1:fishPop$fish$maxAge) {
#     NormProbMat[,i] = Prob2Mat[,i]/ColSums[i] # Normalise probabilities
#   }
#   
#   SizeDatHolder<- matrix(NA,nrow=dim(AgeDatHolder)[1],ncol=dim(AgeDatHolder)[2])
#   for(patch in 1:length(collSizeFI)){
#     if ((FIsurvPatch * collSizeFI)[patch] == 1){
#       TempSizeFreq <- rep(NA,length(1:fishPop$fish$maxAge))
#       TempMat <-  Prob2Mat 
#       for(a in 1:fishPop$fish$maxAge){
#         Na <- AgeDatHolder[which(AgeDatHolder[,patch]==a),patch]
#         TempMat[,a] <- NormProbMat[,a] * length(Na)
#       }
#       TempSizeFreq <- rowSums(TempMat)
#       TempSizeDat <- NULL
#       for(i in 1:length(LenMids)){
#         TempSizeDat <- append( TempSizeDat, rep(LenMids[i],round(TempSizeFreq[i])))
#       }
#       if(length(TempSizeDat) > dim(SizeDatHolder)[1]){TempSizeDat <- TempSizeDat[1:dim(SizeDatHolder)[1]]}
#       SizeDatHolder[1:length(TempSizeDat),patch] <- TempSizeDat
#     }
#   }
#   return(SizeDatHolder)
# }
# 
# 
