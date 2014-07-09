#####################################
#Sampling Function for SNAP working group
# 
#Collects and Aggregates (optional) 10 different kinds of data. 
#######################################################################
SpaceNumAtAgeT<-array(dim=c((simTime+burn),SpaceR,SpaceC,kmax))
# Sample this in survey

CatchByFisher<-matrix(ncol=simTime-burn,nrow=Fishers)
totCatch <- sums catch by fisher

CollectData = function(timeStep, simTime, burn, FFleetvec, Aggregate, fish, NoTakeZone, SpaceC, SpaceR){
  
  # Set some stuff up
  nyrs <- simTime-burn
  Yr <- timeStep
  
  
  ### Conduct FD sampling for each fleet
  for(FFleet in FFleetvec){
    
    ####Collect historical data for current year
    histCatchDat <- CollectHistCatchData(Yr,nyrs,FFleet,histCatch_FD,histStartYr,Aggregate)
    histCPUEDat <- CollectHistCPUEData(Yr,nyrs,FFleet,histEffort_FD,histStartYr,Aggregate)
    
    #### Collect FD data from these patches
    FDCatch <- CollectFisheryCatch(Yr,FFleet,Aggregate)
    FDCPUE <- CollectFisheryCPUE(Yr,FFleet,collEffortFD,Aggregate)
    AgeFDMat <- CollectFisheryAges(Yr,FFleet,nAgeFD,collAgeFD,fishPop,FDsurvPatch)
    SizeFDMat <- CollectFisherySizes(Yr,FFleet,nSizeFD,collSizeFD,FDsurvPatch,fishPop,minLen,maxLen,LengthBins,Sel50,Sel95)
  }
  
    ####Survey population
    Survey <- GoSurvey(Yr,nyrs,fishPop,survq,survEffort,SurvPatch,Sel50FI,Sel95FI,FFleet,minLen,maxLen,LengthBins)
  
    
    
    
    
  #### Identify patches that have a catch to sample
  FDsurvPatch<-rep(0,length(collAgeFI))
  FDsurvPatch[colSums(as.matrix(FFleet$CatchatAge[Yr,,])) > 0] <- 1
  
 
  #### Identify patches that have a survey to sample
  FIsurvPatch <- rep(0,length(collAgeFI))
  FIsurvPatch[colSums(as.matrix(Survey$CatchatAge_FI)) > 0] <- 1
  
  #### Collect FI data from these patches
  CatchFI <- CollectCatchFI(Yr,nyrs,Survey$CatchatSize_FI,collCatchFI,FFleet,Aggregate,minLen, maxLen, LengthBins)
  CPUEFI <- CollectSurvCPUE(Survey$CatchatSize_FI,survEffort,SurvPatch,collCatchFI,FFleet,Aggregate,minLen, maxLen, LengthBins)
  numFI <- CollectSurvnum(Survey$CatchatSize_FI,survEffort,SurvPatch,collCatchFI,FFleet,Aggregate)
  AgeFIMat <- CollectSurvAges(Yr,collAgeFI,fishPop,nAgeFI,Survey,FIsurvPatch)
  SizeFIMat <- CollectSurvSize(collSizeFI,FIsurvPatch,minLen,mxLen,Sel50FI,Sel95FI,fishPop,LengthBins,Survey,nSizeFI)
  VBDat <- CollectVBData(Yr,VBbins,collSizeFI,fishPop,survq,SurvPatch,Sel50VB,Sel95VB,minLen,maxLen,FIsurvPatch,numVB)
  
  if (Aggregate == FALSE){
    AgeFI <- AgeFIMat
    SizeFI <- SizeFIMat
    AgeFD <- AgeFDMat
    SizeFD <- SizeFDMat
  } else {
    TempAgeFI <- as.vector(AgeFIMat)
    AgeFI <- TempAgeFI[!is.na(TempAgeFI)]
    TempSizeFI <- as.vector(SizeFIMat)
    SizeFI <- TempSizeFI[!is.na(TempSizeFI)]
    TempAgeFD <- as.vector(AgeFDMat)
    AgeFD <- TempAgeFD[!is.na(TempAgeFD)]
    TempSizeFD <- as.vector(SizeFDMat)
    SizeFD <- TempSizeFD[!is.na(TempSizeFD)]
  }
  
  ####Return list of collected data
  DataOut<-NULL
  DataOut$histCatchDat <- histCatchDat
  DataOut$histCPUEDat <- histCPUEDat
  DataOut$FDCatch <- FDCatch
  DataOut$FDCPUE<-FDCPUE 
  DataOut$AgeFD <- AgeFD
  DataOut$SizeFD <- SizeFD    
  DataOut$CatchFI <- CatchFI
  DataOut$CPUEFI <- CPUEFI
  DataOut$numFI <- numFI
  DataOut$AgeFI <- AgeFI
  DataOut$SizeFI <- SizeFI
  DataOut$VBDat <- VBDat
  
  
  return(DataOut)
}


CollectHistCatchData = function(Yr,nyrs,FFleet,histCatch_FD,histStartYr,Aggregate){  
  dYr <- Yr-nyrs
  if (dYr == 1 & histCatch_FD == 1){
    CatchHist_FD <- FFleet$totCatch[histStartYr:Yr,]  
    if(Aggregate == FALSE){
      Final <- CatchHist_FD
    } else {   #sum across all patches for all years
      Final <- rowSums(CatchHist_FD,na.rm=TRUE )
    }
  } else {Final <- NA}
  return(Final) 
}

CollectHistCPUEData = function(Yr,nyrs,FFleet,histEffort_FD,histStartYr,Aggregate){  
  dYr <- Yr-nyrs  
  if (dYr == 1 & histEffort_FD == 1){
    CPUEHist_FD <- FFleet$Catch[histStartYr:Yr,]/FFleet$effort[histStartYr:Yr,]
    if(Aggregate == FALSE){
      Final <- CPUEHist_FD
    } else {   #sum across all patches for all years
      Final <- rowSums(FFleet$Catch[histStartYr:Yr,],na.rm=TRUE)/ rowSums(FFleet$effort[histStartYr:Yr,],na.rm=TRUE)
    }
  } else {Final <- NA}
  return(Final)
}

#Determine Catch at age from FI survey
GoSurvey = function(Yr,nyrs,fishPop,survq,survEffort,SurvPatch,sel50,sel95,FFleet,minLen,maxLen,LengthBins){
  #Calculate survey selectivity at age
  selFI<-ConvertSel(sel50,sel95,fishPop$fish$lengthSD,fishPop)
  Fsurv <- survq*((survEffort/length(SurvPatch))*SurvPatch)
  Harv <- matrix((1-exp(-Fsurv)),nrow = dim(selFI)[1],ncol = dim(selFI)[2],byrow=T)
  CatchatAge_FI <- pmin(fishPop$N[Yr,,], fishPop$N[Yr,,]*(selFI*Harv))
  
  #Calculate the probability of falling into length bins
  LengthVec = seq(minLen, maxLen, by= LengthBins)   #changed to make 1cm length bins
  LenMids = seq(LengthVec[1] +((LengthVec[2]-LengthVec[1])/2), by=(LengthVec[2]-LengthVec[1]), length=length(LengthVec)-1)
  AgeLenProbMat = matrix(0, nrow=length(LengthVec)-1, ncol=fishPop$fish$maxAge) 
  for (i in seq_along(1:fishPop$fish$maxAge)) 
  {
    for (lengthclass in seq_along(LengthVec)[-1])
    {  
      length_cat = lengthclass-1
      if(length_cat==1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], fishPop$fish$length[i], 
                                                            fishPop$fish$lengthSD[i]) 
      if(length_cat>1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], fishPop$fish$length[i], 
                                                           fishPop$fish$lengthSD[i]) - pnorm(LengthVec[lengthclass-1], fishPop$fish$length[i], fishPop$fish$lengthSD[i])       
    }
    AgeLenProbMat[,i] = AgeLenProbMat[,i]/sum(AgeLenProbMat[,i])
  }
  
  ##Calculate the selectivity at length for each length bin
  SelLengthVec = plogis(LenMids, location=sel50, scale=(sel95-sel50)/log(19))
  #SelLengthVec = rep(1,56)
  
  ###Create matrix of normalized probabilities of being above the size and selected.
  Prob2Mat = AgeLenProbMat* SelLengthVec    #SelLengthMat
  NormProbMat = Prob2Mat
  ColSums = apply(Prob2Mat, 2, sum)
  for (i in 1:fishPop$fish$maxAge) {
    NormProbMat[,i] = Prob2Mat[,i]/ColSums[i] # Normalise probabilities
  }
  
  #Calculate the Catch at Size
  CatchatSize_FI <- matrix(NA,nrow=length(LenMids),ncol=fishPop$nPatch)
  for(i in 1:fishPop$nPatch){
    CatchLengthMat <-  CatchatAge_FI[,i] * t(NormProbMat)
    CatchatSize_FI[,i] <- apply(CatchLengthMat, 2, sum)
  }
  
  Outs<-NULL
  Outs$CatchatAge_FI<-CatchatAge_FI
  Outs$CatchatSize_FI <- CatchatSize_FI
  return(Outs)
}

CollectCatchFI <- function(Yr,firstDatYr,CatchatSizeFI,collCatchFI,FFleet,Aggregate,minLen, maxLen, LengthBins){
  LengthVec = seq(minLen, maxLen, by= LengthBins)   #changed to make 1cm length bins
  LenMids = seq(LengthVec[1] +((LengthVec[2]-LengthVec[1])/2), by=(LengthVec[2]-LengthVec[1]), length=length(LengthVec)-1)
  WeightAtLength = FFleet$fishPop$fish$wa*LenMids^FFleet$fishPop$fish$wb
  #weight <- FFleet$fishPop$fish$weightAtAge()
  if (length(collCatchFI) > 1){
    Catch_FI <- colSums(as.matrix(CatchatSizeFI[1:(length(LengthVec)-1),]*WeightAtLength))
  }else{
    Catch_FI <- sum(CatchatSizeFI[1:(length(LengthVec)-1)]*WeightAtLength)
  }
  if (Aggregate == FALSE){
    Final <- Catch_FI
  } else {
    Final <- sum(Catch_FI,na.rm=TRUE)
  }
  return(Final)
}

#Determine survey CPUE
CollectSurvCPUE <- function(CatchatSizeFI,survEffort,SurvPatch,collCatchFI,FFleet,Aggregate,minLen, maxLen, LengthBins){
  LengthVec = seq(minLen, maxLen, by= LengthBins)   #changed to make 1cm length bins
  LenMids = seq(LengthVec[1] +((LengthVec[2]-LengthVec[1])/2), by=(LengthVec[2]-LengthVec[1]), length=length(LengthVec)-1)
  WeightAtLength = FFleet$fishPop$fish$wa*LenMids^FFleet$fishPop$fish$wb
  #SurvEffTemp <- replace(((survEffort/length(SurvPatch))*SurvPatch), survEffort*SurvPatch == 0, NA) #avoiding dividing by zero
  SurvEffTemp <- replace(((survEffort/sum(SurvPatch))*SurvPatch), survEffort*SurvPatch == 0, NA) 
  if (length(collCatchFI) >1){
    SurvCatch <- colSums(as.matrix(CatchatSizeFI[1:(length(LengthVec)-1),]*WeightAtLength))
  }else{
    SurvCatch <- sum(CatchatSizeFI[1:(length(LengthVec)-1)]*WeightAtLength) 
  }
  CPUE_FI <- SurvCatch/SurvEffTemp
  if (Aggregate == FALSE){
    Final <- CPUE_FI
  } else {
    Final <- sum(SurvCatch,na.rm=TRUE)/sum(SurvEffTemp,na.rm=TRUE)
  }
  return(Final)
} 

#collect CPUE in numbers by size -- needed in dtree
CollectSurvnum <- function(CatchatSizeFI,survEffort,SurvPatch,collCatchFI,FFleet,Aggregate){
  SurvEffTemp <- replace(((survEffort/sum(SurvPatch))*SurvPatch), survEffort*SurvPatch == 0, NA) #avoiding dividing by zero, and effort is only allocated to patches that are surveyed
  # SurvEffTemp <- replace(((survEffort/length(SurvPatch))*SurvPatch), survEffort*SurvPatch == 0, NA) #avoiding dividing by zero
  if (length(collCatchFI) >1){
    #SurvCatchNums <- as.matrix(CatchatSizeFI[1:FFleet$NumSizes,which(collCatchFI==1)])
    SurvCatchNums <- as.matrix(CatchatSizeFI[1:FFleet$NumSizes,])
  }else{
    #SurvCatchNums <- CatchatSizeFI[1:FFleet$NumSizes,]*collCatchFI  
    SurvCatchNums <- CatchatSizeFI[1:FFleet$NumSizes,]
  }
  numCPUE_FI <- round(SurvCatchNums,0)   ##Converts to whole numbers   
  
  SurvEffMat <- matrix(SurvEffTemp,nrow=dim(numCPUE_FI)[1],ncol=dim(numCPUE_FI)[2],byrow=T)
  if (Aggregate == FALSE){
    Final <- numCPUE_FI/SurvEffMat
  } else {
    Final <- rowSums(numCPUE_FI,na.rm=TRUE)/sum(SurvEffTemp,na.rm=TRUE)
  }
  return(Final)
}

#Collect Age Data -- assign ages to size data
CollectSurvAges <- function(Yr,collAgeFI,fishPop,nAgeFI,Survey,FIsurvPatch){
  CatchatAge_FI <-  as.matrix(round(Survey$CatchatAge))  #turns the catch at size into whole numbers
  AgeDatHolder <- matrix(NA,nrow=max(nAgeFI,1),ncol=length(collAgeFI))
  for(patch in 1:length(collAgeFI)){
    if ((FIsurvPatch * collAgeFI)[patch] == 1){ 
      AgeDatHolder[,patch] <- sample(1:fishPop$fish$maxAge,nAgeFI,replace=TRUE,prob=(CatchatAge_FI[,patch]/sum(CatchatAge_FI[,patch])))
      
    }
  }
  return(AgeDatHolder)
}

#Collect VB data -- age and size are paired. 
CollectVBData <- function(Yr,VBbins,collSizeFI,fishPop,survq,SurvPatch,Sel50VB,Sel95VB,minLen,maxLen,FIsurvPatch,numVB){
  #Calculate survey selectivity at age
  survEffort <- 25
  selFI<-ConvertSel(Sel50VB,Sel95VB,fishPop$fish$lengthSD,fishPop)
  Fsurv <- survq*((survEffort/length(SurvPatch))*SurvPatch)
  Harv <- matrix((1-exp(-Fsurv)),nrow = dim(selFI)[1],ncol = dim(selFI)[2],byrow=T)
  CatchatAge_FI <- round(pmin(fishPop$N[Yr,,], fishPop$N[Yr,,]*(selFI*Harv)))
  
  MatSize <- max(colSums(CatchatAge_FI))
  AgeDatHolder <- matrix(NA,nrow=MatSize,ncol=length(collSizeFI))
  for(patch in 1:length(collSizeFI)){
    if ((FIsurvPatch * collSizeFI)[patch] == 1){ 
      vec <- NULL
      for(a in 1:fishPop$fish$maxAge){
        vec <- append(vec,rep(a,CatchatAge_FI[a,patch]))
      }
      AgeDatHolder[1:length(vec),patch]<-vec
    }
  }
  
  # Create age-length probability matrix
  LengthVec = seq(minLen, maxLen, by= VBbins)   #changed to make 1cm length bins
  LenMids = seq(LengthVec[1] +((LengthVec[2]-LengthVec[1])/2), by=(LengthVec[2]-LengthVec[1]), length=length(LengthVec)-1)
  AgeLenProbMat = matrix(0, nrow=length(LengthVec)-1, ncol=fishPop$fish$maxAge) 
  for (i in seq_along(1:fishPop$fish$maxAge)) 
  {
    for (lengthclass in seq_along(LengthVec)[-1])
    {  
      length_cat = lengthclass-1   #Calculate the probability of falling into each length bin
      if(length_cat==1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], fishPop$fish$length[i],  fishPop$fish$lengthSD[i]) 
      if(length_cat>1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], fishPop$fish$length[i],  fishPop$fish$lengthSD[i]) - pnorm(LengthVec[lengthclass-1], fishPop$fish$length[i], fishPop$fish$lengthSD[i])       
    }
    # For each Age normalize across lengths
    AgeLenProbMat[,i] = AgeLenProbMat[,i]/sum(AgeLenProbMat[,i])
  }
  
  ##Calculate the selectivity at length for each length bin
  SelLengthVec = plogis(LenMids, location=Sel50VB, scale=(Sel95VB-Sel50VB)/log(19))
  #SelLengthVec = rep(1,34)
  
  ###Create matrix of normalized probabilities of being above the size and selected.
  Prob2Mat = AgeLenProbMat* SelLengthVec    #SelLengthMat
  
  ## Need to re-normalized across all lengths after selectivity was accounted for 
  NormProbMat = Prob2Mat
  ColSums = apply(Prob2Mat, 2, sum)
  for (i in 1:fishPop$fish$maxAge) {
    NormProbMat[,i] = Prob2Mat[,i]/ColSums[i] # Normalise probabilities
  }
  
  SizeDatHolder<-  matrix(NA,nrow=dim(AgeDatHolder)[1],ncol=dim(AgeDatHolder)[2])
  for(patch in 1:length(collSizeFI)){
    if ((FIsurvPatch * collSizeFI)[patch] == 1){
      TempSize <- NULL
      for(a in 1:fishPop$fish$maxAge){
        Na <- length(AgeDatHolder[which(AgeDatHolder[,patch]==a),patch])
        TempSize <- append(TempSize,sample(LenMids, Na, replace=TRUE, prob=NormProbMat[,a]))
      }
    } else {
      TempSize <- rep(NA,dim(AgeDatHolder)[1])
    }
    SizeDatHolder[1:length(TempSize),patch] <- TempSize
  }
  
  # Sample VB data
  SampInd <- sample(1:MatSize,NumVB,replace=TRUE)
  AgeDatHolderVB <- matrix(NA,nrow=MatSize,ncol=length(collSizeFI))
  SizeDatHolderVB <- matrix(NA,nrow=MatSize,ncol=length(collSizeFI))
  AgeDatHolderVB <- AgeDatHolder[SampInd,]
  SizeDatHolderVB <- SizeDatHolder[SampInd,]
  
  Outs<-NULL
  Outs$Age <- AgeDatHolderVB
  Outs$Size <- SizeDatHolderVB
  return(Outs)
}

#Collect FI Size Data
CollectSurvSize <- function(collSizeFI,FIsurvPatch,minLen,mxLen,Sel50FI,Sel95FI,fishPop,LengthBins,Survey,nSizeFI){
  CatchatAge_FI <-  as.matrix(round(Survey$CatchatAge))  #turns the catch at size into whole numbers
  AgeDatHolder <- matrix(NA,nrow=max(nSizeFI,1),ncol=length(collSizeFI))
  for(patch in 1:length(collSizeFI)){
    if ((FIsurvPatch * collSizeFI)[patch] == 1){ 
      AgeDatHolder[,patch] <- sample(1:fishPop$fish$maxAge,nSizeFI,replace=TRUE,prob=(CatchatAge_FI[,patch]/sum(CatchatAge_FI[,patch])))
      
    }
  }
  
  # Create age-length probability matrix
  LengthVec = seq(minLen, maxLen, by= LengthBins)   #changed to make 1cm length bins
  LenMids = seq(LengthVec[1] +((LengthVec[2]-LengthVec[1])/2), by=(LengthVec[2]-LengthVec[1]), length=length(LengthVec)-1)
  AgeLenProbMat = matrix(0, nrow=length(LengthVec)-1, ncol=fishPop$fish$maxAge) 
  for (i in seq_along(1:fishPop$fish$maxAge)) 
  {
    for (lengthclass in seq_along(LengthVec)[-1])
    {  
      length_cat = lengthclass-1   #Calculate the probability of falling into each length bin
      #if(length_cat==1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], fishPop$fish$length[i],  10) 
      if(length_cat==1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], fishPop$fish$length[i],  fishPop$fish$lengthSD[i]) 
      #if(length_cat>1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], fishPop$fish$length[i],  10) - pnorm(LengthVec[lengthclass-1], fishPop$fish$length[i], 12)       
      if(length_cat>1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], fishPop$fish$length[i],  fishPop$fish$lengthSD[i]) - pnorm(LengthVec[lengthclass-1], fishPop$fish$length[i], fishPop$fish$lengthSD[i])       
    }
    # For each Age normalize across lengths
    AgeLenProbMat[,i] = AgeLenProbMat[,i]/sum(AgeLenProbMat[,i])
  }
  
  ##Calculate the selectivity at length for each length bin
  SelLengthVec = plogis(LenMids, location=Sel50FI, scale=(Sel95FI-Sel50FI)/log(19))
  #SelLengthVec = rep(1,34)
  
  ###Create matrix of normalized probabilities of being above the size and selected.
  Prob2Mat = AgeLenProbMat* SelLengthVec    #SelLengthMat
  
  ## Need to re-normalized across all lengths after selectivity was accounted for 
  NormProbMat = Prob2Mat
  ColSums = apply(Prob2Mat, 2, sum)
  for (i in 1:fishPop$fish$maxAge) {
    NormProbMat[,i] = Prob2Mat[,i]/ColSums[i] # Normalise probabilities
  }
  
  SizeDatHolder<- matrix(NA,nrow=dim(AgeDatHolder)[1],ncol=dim(AgeDatHolder)[2])
  for(patch in 1:length(collSizeFI)){
    if ((FIsurvPatch * collSizeFI)[patch] == 1){
      TempSizeFreq <- rep(NA,length(1:fishPop$fish$maxAge))
      TempMat <-  Prob2Mat 
      for(a in 1:fishPop$fish$maxAge){
        Na <- AgeDatHolder[which(AgeDatHolder[,patch]==a),patch]
        TempMat[,a] <- NormProbMat[,a] * length(Na)
      }
      TempSizeFreq <- rowSums(TempMat)
      TempSizeDat <- NULL
      for(i in 1:length(LenMids)){
        TempSizeDat <- append( TempSizeDat, rep(LenMids[i],round(TempSizeFreq[i])))
      }
      if(length(TempSizeDat) > dim(SizeDatHolder)[1]){TempSizeDat <- TempSizeDat[1:dim(SizeDatHolder)[1]]}
      SizeDatHolder[1:length(TempSizeDat),patch] <- TempSizeDat
    }
  }
  return(SizeDatHolder)
}


CollectFisheryCatch <- function(Yr,FFleet,Aggregate){
  FisheryCatch <-FFleet$Catch[Yr,]
  if (Aggregate == FALSE){
    Final <- FisheryCatch
  } else {
    Final <- sum(FisheryCatch,na.rm=TRUE)
  }
  return(Final)
}

CollectFisheryCPUE <- function(Yr,FFleet,collEffortFD,Aggregate){
  FisheryCatch <- FFleet$Catch[Yr,]
  Eff_Temp <- FFleet$effort[Yr,]*collEffortFD
  Effort_FD <- replace(Eff_Temp, Eff_Temp == 0, NA)  #avoiding divide by zero problem in unsampled patches
  CPUE_FD <- FisheryCatch/Effort_FD   ##Collect CPUE if effort is collected in patch. 
  if (Aggregate == FALSE){
    Final <- CPUE_FD
  } else {
    Final <- sum(FisheryCatch,na.rm=TRUE)/sum(Effort_FD,na.rm=TRUE)
  }
  return(Final)
}

#Collect Age Data
CollectFisheryAges <- function(Yr,FFleet,nAgeFD,collAgeFD,fishPop,FDsurvPatch){
  CatchatAge_FD <-  as.matrix(round(FFleet$CatchatAge[Yr,,]))  #turns the catch at age into whole numbers
  AgeDatHolder <- matrix(NA,nrow=max(nAgeFD,1),ncol=length(collAgeFD))
  for(patch in 1:length(collAgeFD)){
    if ((FDsurvPatch * collAgeFD)[patch] == 1){ 
      AgeDatHolder[,patch]<-sample(1:fishPop$fish$maxAge,nAgeFD,replace=TRUE,prob=(CatchatAge_FD[,patch]/sum(CatchatAge_FD[,patch])))
    }
  }
  return(AgeDatHolder)
}

#Collect FD Size Data 
CollectFisherySizes <- function(Yr,FFleet,nSizeFD,collSizeFD,FDsurvPatch,fishPop,minLen,maxLen,LengthBins,sel50,sel95){
  CatchatAge_FD <-  as.matrix(round(FFleet$CatchatAge[Yr,,]))  #turns the catch at age into whole numbers
  AgeDatHolder <- matrix(NA,nrow=max(nSizeFD,1),ncol=length(collSizeFD))
  for(patch in 1:length(collSizeFD)){
    if ((FDsurvPatch * collSizeFD)[patch] == 1){ 
      AgeDatHolder[,patch]<-sample(1:fishPop$fish$maxAge,nSizeFD,replace=TRUE,prob=(CatchatAge_FD[,patch]/sum(CatchatAge_FD[,patch])))
    }
  }
  
  # Create age-length probability matrix
  LengthVec = seq(minLen, maxLen, by= VBbins)   #changed to make 1cm length bins
  LenMids = seq(LengthVec[1] +((LengthVec[2]-LengthVec[1])/2), by=(LengthVec[2]-LengthVec[1]), length=length(LengthVec)-1)
  AgeLenProbMat = matrix(0, nrow=length(LengthVec)-1, ncol=fishPop$fish$maxAge) 
  for (i in seq_along(1:fishPop$fish$maxAge)) 
  {
    for (lengthclass in seq_along(LengthVec)[-1])
    {  
      length_cat = lengthclass-1   #Calculate the probability of falling into each length bin
      if(length_cat==1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], fishPop$fish$length[i],  fishPop$fish$lengthSD[i])     
      if(length_cat>1) AgeLenProbMat[length_cat,i] = pnorm(LengthVec[lengthclass], fishPop$fish$length[i],  fishPop$fish$lengthSD[i]) - pnorm(LengthVec[lengthclass-1], fishPop$fish$length[i], fishPop$fish$lengthSD[i])       
    }
    # For each Age normalize across lengths
    AgeLenProbMat[,i] = AgeLenProbMat[,i]/sum(AgeLenProbMat[,i])
  }
  
  ##Calculate the selectivity at length for each length bin
  SelLengthVec = plogis(LenMids, location=sel50, scale=(sel95-sel50)/log(19))
  #SelLengthVec = rep(1,34)
  
  ###Create matrix of normalized probabilities of being above the size and selected.
  Prob2Mat = AgeLenProbMat* SelLengthVec    #SelLengthMat
  
  ## Need to re-normalized across all lengths after selectivity was accounted for 
  NormProbMat = Prob2Mat
  ColSums = apply(Prob2Mat, 2, sum)
  for (i in 1:fishPop$fish$maxAge) {
    NormProbMat[,i] = Prob2Mat[,i]/ColSums[i] # Normalise probabilities
  }
  
  SizeDatHolder<- matrix(NA,nrow=dim(AgeDatHolder)[1],ncol=dim(AgeDatHolder)[2])
  for(patch in 1:length(collSizeFD)){
    if ((FDsurvPatch * collSizeFD)[patch] == 1){
      TempSizeFreq <- rep(NA,length(1:fishPop$fish$maxAge))
      TempMat <-  Prob2Mat 
      for(a in 1:fishPop$fish$maxAge){
        Na <- AgeDatHolder[which(AgeDatHolder[,patch]==a),patch]
        TempMat[,a] <- NormProbMat[,a] * length(Na)
      }
      TempSizeFreq <- rowSums(TempMat)
      TempSizeDat <- NULL
      for(i in 1:length(LenMids)){
        TempSizeDat <- append( TempSizeDat, rep(LenMids[i],round(TempSizeFreq[i])))
      }
      if(length(TempSizeDat) > dim(SizeDatHolder)[1]){TempSizeDat <- TempSizeDat[1:dim(SizeDatHolder)[1]] }
      SizeDatHolder[1:length(TempSizeDat),patch] <- TempSizeDat
    }
  }
  return(SizeDatHolder)
}

