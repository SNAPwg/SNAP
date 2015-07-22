#==this function runs the main population and economic dynamics for simulations given several inputs (usually in the form of csvs)

Master<-function(Life,SimCTL,Fleets,Species, Taxa, season,Samp,ManageStrats,Management,NoTakeZoneInit,NoTakeZoneImp,habitat,Graphs=F,GraphsFish=F,PrintLifeHistory=F)
{

  BaseLife<- Life

  BaseFleets<- Fleets

  BaseSeason<- season

  DataParams<- ShapeDatParams(Samp) #didn't realize these were included in Fishery below...

  Results<- list()

  Assessments <- list()

  Iterations<- SimCTL[grepl('NumIterations',SimCTL[,2]),1]

  for (s in 1:(Iterations))
  {


    show(paste(round(100*(s/Iterations),2),'% Done with Iterations',sep=''))
    if (s>1)
    {

      Life<- BaseLife

      Fleets<- BaseFleets

      season<-   BaseSeason

      StochTerms<- ApplyStochasticity(Life,Fleets,SimCTL)

      Life<- StochTerms$Life

      Fleets<- StochTerms$Fleets
    }

    Fishery<- ShapeFishery(Life = Life,Species = Species, Taxa = Taxa, Fleets = Fleets,SimCTL = SimCTL,season = season, NoTakeZone = NoTakeZoneInit,Samp = Samp)


    Fishery$ManagementPlan<- ManageStrats$Plan

    Fishery$Iteration<- s

    Things<- names(Fishery)

    for (t in 1:length(Things))
    {
      eval(parse(text=paste(Things[t],'<- ','Fishery$',Things[t],sep='')))
    }


    if(DataParams$Aggregate == 0){
      FDCatchData <- array(dim=c(simTime-burn+1, SpaceR, SpaceC, FleetN))
      FDCPUEData <- array(dim=c(simTime-burn+1, SpaceR, SpaceC, FleetN))
      FDAgeData <- array(NA,dim=c(simTime-burn+1, DataParams$nAgeFD,SpaceR, SpaceC, FleetN))
      FDSizeData <- array(NA,dim=c(simTime-burn+1, DataParams$nSizeFD,SpaceR, SpaceC, FleetN))
      FICatchData <- array(dim=c(simTime-burn+1, SpaceR, SpaceC))
      FICPUEData <- array(dim=c(simTime-burn+1, SpaceR, SpaceC))
    } else {
      FDCatchData <- array(dim=c(simTime-burn+1, FleetN))
      FDCPUEData <- array(dim=c(simTime-burn+1, FleetN))
      FDAgeData <- array(NA,dim=c(simTime-burn+1, DataParams$nAgeFD*SpaceR*SpaceC, FleetN))
      FDSizeData <- array(NA,dim=c(simTime-burn+1, DataParams$nSizeFD*SpaceR*SpaceC, FleetN))
      FICatchData <- array(dim=c(simTime-burn+1))
      FICPUEData <- array(dim=c(simTime-burn+1))
    }


    PortLc<-c(PortX,PortY) #Does this get used anywhere?  I DON'T KNOW.


    #==dummy spatial matrix
    SpaceUse<-matrix(0,ncol=SpaceC,nrow=SpaceR)

    if(PrintLifeHistory==T)
    {
      #==print out maturity, selectivity, growth, etc.
      pdf(file=paste(FigureFolder,'/LifeHistory ',ManageStrats$Plan,' Iteration ',s,'.pdf',sep=''))
      par(mfrow=c(2,1),mar=c(0.1,4,.1,4),oma=c(4,0,0,0))
      plot(MatAge,type="l",xaxt='n',las=2,ylab='')
      mtext(side=2,"Probability",line=3)
      lines(FishSel,lty=2,col=2)
      lines(MoveProb,lty=3,col=3)
      legend("bottomright",lty=c(1,2,3),col=c(1,2,3),legend=
               c("Maturity","Selectivity","Movement Probability"),bty='n')
      plot(lenAtAge,type="l",las=1,ylab='',ylim=c(0,max(lenAtAge)))
      par(new=T)
      plot(wgtAtAge,type="l",lty=2,col=2,yaxt='n',ylab='',ylim=c(0,max(wgtAtAge)))
      legend("bottomright",lty=c(1,2),col=c(1,2),legend=
               c("Length","Weight"),bty='n')
      axis(side=4,las=2)
      mtext(side=4,"Weight (kg)",line=2)
      mtext(side=2,"Length (cm)",line=2)
      mtext(side=1,"Age (yr)",line=2)

      filled.contour(matrix(as.numeric(unlist(habitat)),ncol=SpaceC),x=seq(1,SpaceR),y=seq(1,SpaceC))
      mtext(side=3,"Habitat quality")
      filled.contour(matrix(as.numeric(unlist(NoTakeZone)),ncol=SpaceC),x=seq(1,SpaceR),y=seq(1,SpaceC))
      mtext(side=3,"No take zones (0)")
      dev.off()
    }

    #===========================================================
    #== Begin simulations====================================
    #===========================================================

    #==population dynamics array tracking number at age in a patch by year
    SpaceNumAtAgeT      <-array(dim=c((simTime),SpaceR,SpaceC,kmax))
    SpaceCatAgeByFisher <-array(dim=c(simTime-burn+1,SpaceR,SpaceC,kmax,max(Fishers),FleetN)) #in numbers
    SpaceEffort         <-array(0,dim=c(simTime-burn+1,SpaceR,SpaceC,FleetN))

    CatchByFisher       <-array(dim=c(max(Fishers),(simTime-burn+1),FleetN))
    ProfitByFisher      <-array(dim=c(max(Fishers),(simTime-burn+1),FleetN))
    CostByFisher        <-array(dim=c(max(Fishers),simTime-burn+1,FleetN))
    SpawningBiomass     <-rep(0,simTime-burn)
    NumbersExplt        <-rep(0,simTime-burn)
    CostOfManagement    <-rep(0,simTime-burn)
    InsideMPAspbio      <-rep(0,simTime-burn)
    OutsideMPAspbio     <-rep(0,simTime-burn)

    #==data.frames for data from assessment
    CatchData           <-data.frame(Species=character(0),TimeStep=numeric(0),Catch=numeric(0),Units=character(0))
    IndexData           <-data.frame(Species=character(0),TimeStep=numeric(0),Index=numeric(0),Units=character(0),MPA=numeric(0),Inside=numeric(0),PropMPA=numeric(0))
    LengthData          <-data.frame(Species=character(0),TimeStep=numeric(0),Length=numeric(0),Units=character(0),MPA=numeric(0),Inside=numeric(0),PropMPA=numeric(0))

    #==index for spatial area==========================
    coords<-NULL
    for(x in 1:SpaceR)
    {
      temp<-cbind(rep(x,SpaceC),seq(1,SpaceC))
      coords<-rbind(coords,temp)
    }

    #==make movement array==============================
    #=should there be two movement arrays?
    #=one for adults and one for recruitment?
    Movement<-movArray(SpaceR,SpaceC,sdx,sdy)

    #==initialize population and allow movement to equilibrate
    SpaceNumAtAgeT[1:burn,,,]<-InitialPop(R0,M,kmax,SpaceR,SpaceC,MoveProb,coords,Movement,Graphs=F,burn)

    #==calculate virgin spawning biomass (1843.584)================
    tempSpBio<-rep(0,kmax)
    for(i in 1:kmax)
      tempSpBio[i]<-sum(SpaceNumAtAgeT[burn,,,i]*MatAge[i]*wgtAtAge[i])

    VirSpBio<-sum(tempSpBio)
    Fishery$VirSpBio<- VirSpBio

    #==relate the cost of fishing to the number of fish in a patch?? (e.g. density is low == must fish harder)

    if(Graphs==T | GraphsFish==T)
    {
      #setwd("C:/Users/Cody/Desktop/spatial gif")
      #dev.new(height=7,width=7)
      # 	    ani.record(reset=TRUE)
      #       ani.options(interval=.01)
    }

    for(timeStep in burn:simTime)
    {
      NoTakeZone<-NoTakeZoneInit
      if(timeStep>=initManage)
      {
        NoTakeZone<-NoTakeZoneImp

        Fishery<- ApplyManagement(Strat=ManageStrats,Management=Management,Fishery=Fishery)

        Things<- names(Fishery)

        for (t in 1:length(Things))
        {
          eval(parse(text=paste(Things[t],'<- ','Fishery$',Things[t],sep='')))
        }

      }
      if(sum(dim(NoTakeZone)- c(SpaceR,SpaceC))!=0)
      {
        print("no take zone is of the wrong dimensions")
        break
      }
      for(flt in 1:FleetN)
      {
        #==fishers decide where to fish============
        #==initial costs
        Distance<-sqrt((row(SpaceUse)-PortX[flt])^2 + (col(SpaceUse)-PortY[flt])^2)

        tempBenefit<-array(dim=c(SpaceR,SpaceC,kmax))
        for(x in 1:nrow(FishSel))
        {
          tempBenefit[,,x]<-SpaceNumAtAgeT[timeStep,,,x]*FishSel[x,flt]*wgtAtAge[x]
        }
        BenefitPatch<-apply(tempBenefit,c(1,2),sum)

        if(timeStep==burn)OrigMaxBen<-max(BenefitPatch)

        #==costs depend both on the distance from port and the density of fish
        #tempCostFish<-costFish[flt]*(OrigMaxBen/BenefitPatch)^costSteep
        tempCostFish<-(1+(costSteep-(costSteep)*(BenefitPatch/OrigMaxBen)))*costFish[flt]
        CostPatch<-Distance*costTrv[flt]+tempCostFish
        #==cap potential benefit at what can be removed after calculating cost of fishing
        BenefitPatch[BenefitPatch>maxCapac[flt]]<-maxCapac[flt]
        #==translate to dollars
        BenefitPatch<-BenefitPatch*price[flt]
        #==in dollars (through selected kg of biomass), this is what a fisher would catch by going to a given patch
        BenefitPatch<-BenefitPatch*NoTakeZone
        NetBenefitPatch<-BenefitPatch-CostPatch
        if(timeStep==burn)OrigMaxNetBen<-max(NetBenefitPatch)

        tempSNALT<-SpaceNumAtAgeT[timeStep,,,]
        if(timeStep==burn)OrigTempSNALT<-max(tempSNALT[,,3])
        #==all fishers select their patch and fish if the season is 'in'

        # Go Fishing -----------------------------------------------------------------------

        if(any(timeStep%%12==(which(!is.na(season[,flt+1]))-1)))
        {
          for(f in 1:Fishers[flt])
          {
            #==find closest, highest value, fishable patch=================
            #==insert option for choosing randomly in a bivariate normal centered on the port of entry
            #==insert method to close patches here=========================
            #==open access, fishers increase until profit == 0 ??
            maxNetBene<-which(NetBenefitPatch==max(NetBenefitPatch),arr.ind=T)
            chosenPatch<-maxNetBene[which(Distance[maxNetBene]==min(Distance[maxNetBene])),]

            if(NetBenefitPatch[chosenPatch[1],chosenPatch[2]]>0)
            {
              #===FISHHHHH===#=><">====><"">====><""">===============================
              potentialCatch<-sum(SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel[,flt]*wgtAtAge)
              if(potentialCatch<=maxCapac[flt])
              {
                tempSNALT[chosenPatch[1],chosenPatch[2],]<-SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],] -
                  SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel[,flt]
                #==spatial catch at age by patch by fisher by fleet
                for(p in 1:kmax)
                  SpaceCatAgeByFisher[timeStep-burn+1,chosenPatch[1],chosenPatch[2],p,f,flt]<-SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],p]*FishSel[p,flt]
                CatchByFisher[f,timeStep-burn+1,flt]<-potentialCatch
                CostByFisher[f,timeStep-burn+1,flt]<-CostPatch[chosenPatch[1],chosenPatch[2]]
                ProfitByFisher[f,timeStep-burn+1,flt]<-CatchByFisher[f,timeStep-burn+1,flt]*price[flt] - CostByFisher[f,timeStep-burn+1,flt]
              }
              #==if they can catch more than capactity, they'll probably change selectivity and throw small ones back....think about that
              if(potentialCatch>maxCapac[flt])
              {
                #===find harvest rate that would put the fisher at max capacity===
                maxHarv<-1
                minHarv<-.000000001
                for(x in 1:25)
                {
                  useHarv<-(maxHarv+minHarv)/2
                  tempCat<-sum(SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel[,flt]*useHarv*wgtAtAge)
                  if(tempCat<maxCapac[flt])
                    minHarv<-useHarv
                  if(tempCat>maxCapac[flt])
                    maxHarv<-useHarv
                }
                #==apply harvest===============================================
                tempSNALT[chosenPatch[1],chosenPatch[2],]<-SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],] -
                  SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel[,flt]*useHarv
                CatchByFisher[f,timeStep-burn+1,flt]<-sum(SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],]*FishSel[,flt]*useHarv*wgtAtAge)
                for(p in 1:kmax)
                  SpaceCatAgeByFisher[timeStep-burn+1,chosenPatch[1],chosenPatch[2],p,f,flt]<-SpaceNumAtAgeT[timeStep,chosenPatch[1],chosenPatch[2],p]*FishSel[p,flt]*useHarv

                CostByFisher[f,timeStep-burn+1,flt]<-CostPatch[chosenPatch[1],chosenPatch[2]]
                ProfitByFisher[f,timeStep-burn+1,flt]<-CatchByFisher[f,timeStep-burn+1,flt]*price[flt] - CostByFisher[f,timeStep-burn+1,flt]

              }

              #==tracks number of fishermen fishing in a patch at time
              SpaceEffort[timeStep-burn+1,chosenPatch[1],chosenPatch[2],flt]<-SpaceEffort[timeStep-burn+1,chosenPatch[1],chosenPatch[2],flt]+1

              #==update benefitPatch for next fisher=======================================
              tempBenefit<-array(dim=c(SpaceR,SpaceC,kmax))
              for(x in 1:nrow(FishSel))
                tempBenefit[,,x]<-tempSNALT[,,x]*FishSel[x,flt]*wgtAtAge[x]
              BenefitPatch<-apply(tempBenefit,c(1,2),sum)

              #    tempCostFish<-costFish[flt]*(OrigMaxBen/BenefitPatch)^costSteep
              tempCostFish<-(1+(costSteep-(costSteep)*(BenefitPatch/OrigMaxBen)))*costFish[flt]
              CostPatch<-Distance*costTrv[flt]+tempCostFish
              BenefitPatch[BenefitPatch>maxCapac[flt]]<-maxCapac[flt]
              BenefitPatch<-BenefitPatch*price[flt]
              BenefitPatch<-BenefitPatch*NoTakeZone
              NetBenefitPatch<-BenefitPatch-CostPatch
              if(Graphs==T)
              {
                filled.contour(as.matrix(NetBenefitPatch),y=seq(1,SpaceC),x=seq(1,SpaceR))
                #,zlim=c(-OrigMaxNetBen,OrigMaxNetBen)
                mtext(side=3,paste("Timestep = ",timeStep))
                # ani.record()
              }

            }
          }	#end fishers
        } # end if for season

        MPApresent<-sum(NoTakeZone==0)>0
        PropMPA<-sum(NoTakeZone==0)/(nrow(NoTakeZone)*ncol(NoTakeZone))
        LengthSampN<-1000

        #==if no data collection int hat time step, add a row without data (this is important for some methods)
        #

        if((timeStep %in% (sampleTimeSteps+burn))==FALSE)
        {
          inCat           <-data.frame(Species=Species,TimeStep=timeStep,Catch=NA,Units="Biomass")
          CatchData       <-rbind(CatchData,inCat)

          inIndex         <-data.frame(Species=Species,TimeStep=timeStep,Index=NA,Units="Biomass",MPA=MPApresent,Inside=NA,PropMPA=PropMPA)
          IndexData       <-rbind(IndexData,inIndex)

          inLenDat        <-data.frame(Species=Species,TimeStep=timeStep,Length=NA,Units="mm",MPA=MPApresent,Inside=NA,PropMPA=PropMPA)
          LengthData      <-rbind(LengthData,inLenDat)
        }

        #==if data collection, then deterministic sampling
        #==make into a function eventually
        #           CollecteDataDet<-function(arguments,LengthSampN=5000)
        #           {
        if(timeStep %in% (sampleTimeSteps+burn))
        {
          #==actual catch
          ObservedCatch   <-apply(CatchByFisher,c(2),sum,na.rm=T)[timeStep-burn+1] # catch in biomass (Fishers,year,fleet)
          inCat           <-data.frame(Species=Species,TimeStep=timeStep,Catch=ObservedCatch,Units="Biomass")
          CatchData       <-rbind(CatchData,inCat)


          if(PropMPA>0)
          {
            # sample for lengths inside MPA
            tempNbyAgeAndSpace       <-SpaceNumAtAgeT[timeStep,,,]
            NatAgeInsideMPA          <-tempNbyAgeAndSpace
            for(d in 1:length(FishSel[,flt]))
              NatAgeInsideMPA[,,d]   <-tempNbyAgeAndSpace[,,d]*FishSel[d,1]*as.matrix(-1*(NoTakeZone-1))

            NatAgeInsideMPA[NatAgeInsideMPA==0]<-NA
            totalLengthStruc         <-apply(NatAgeInsideMPA,c(3),sum,na.rm=T)
            sampledAges              <-sample(seq(1,kmax),LengthSampN,prob=totalLengthStruc,replace=TRUE)
            sampledLengths           <-sampledAges*Fishery$lenAtAge[sampledAges]*rnorm(LengthSampN,1,Fishery$lenSD)
            inLenDat                 <-data.frame(Species=rep(Species,LengthSampN),TimeStep=rep(timeStep,LengthSampN),Length=sampledLengths,
                                                  Units=rep("mm",LengthSampN),MPA=rep(MPApresent,LengthSampN),Inside=rep(1,LengthSampN),
                                                  PropMPA=rep(PropMPA,LengthSampN))
            inLenDat                 <-inLenDat[-as.numeric(inLenDat[,3])<0,] # sometimes you get negatives when adding the uncertainty at small sizes
            LengthData               <-rbind(LengthData,inLenDat)

            ExpBioInsideMPA          <-sum(apply(NatAgeInsideMPA,c(3),sum,na.rm=t)*FishSel[,1]*wgtAtAge)
            inExpBio                 <-data.frame(Species=Species,TimeStep=timeStep,Index=ExpBioInsideMPA,Units="Biomass",MPA=MPApresent,Inside=1,PropMPA=PropMPA)
            IndexData                <-rbind(IndexData,inExpBio)
          }


          # sample for lengths and biomass outside MPA
          tempNbyAgeAndSpace       <-SpaceNumAtAgeT[timeStep,,,]
          NatAgeOutsideMPA          <-tempNbyAgeAndSpace
          for(d in 1:length(FishSel[,flt]))
            NatAgeOutsideMPA[,,d]   <-tempNbyAgeAndSpace[,,d]*FishSel[d,1]*as.matrix(NoTakeZone)

          NatAgeOutsideMPA[NatAgeOutsideMPA==0]<-NA
          totalLengthStruc         <-apply(NatAgeOutsideMPA,c(3),sum,na.rm=T)
          sampledAges              <-sample(seq(1,kmax),LengthSampN,prob=totalLengthStruc,replace=TRUE)
          sampledLengths           <-sampledAges*Fishery$lenAtAge[sampledAges]*rnorm(LengthSampN,1,Fishery$lenSD)
          inLenDat                 <-data.frame(Species=rep(Species,LengthSampN),TimeStep=rep(timeStep,LengthSampN),Length=sampledLengths,
                                                Units=rep("mm",LengthSampN),MPA=rep(MPApresent,LengthSampN),Inside=rep(1,LengthSampN),
                                                PropMPA=rep(PropMPA,LengthSampN))
          inLenDat                 <-inLenDat[-as.numeric(inLenDat[,3])<0,] # sometimes you get negatives
          LengthData               <-rbind(LengthData,inLenDat)

          # Index of exploitable biomass. This isn't 'CPUE', because effort is omniscient.
          ExpBioOutsideMPA        <-sum(apply(NatAgeOutsideMPA,c(3),sum,na.rm=t)*FishSel[,1]*wgtAtAge)
          inExpBio                 <-data.frame(Species=Species,TimeStep=timeStep,Index=ExpBioOutsideMPA,Units="Biomass",MPA=MPApresent,Inside=0,PropMPA=PropMPA)
          IndexData               <-rbind(IndexData,inExpBio)
          #           return()
          #           }
        } # SAMPLING TIME STEP

        # Assess -----------------------------------------------------------------------
        # This should be moved up before fishing later

        Monitor_Data <- list (Monitor_CatchData = CatchData, Monitor_CPUEData = IndexData, Monitor_LengthData = LengthData)

        # This "PossibleAssessmen" function will later produce the
        # matrix of assessments and requirements

        #         if (timeStep == simTime)
        #         {
        what_can_we_do <- Find_Possible_Assessments(Fishery = Fishery, Data = Monitor_Data, time_steps = timeStep)

        PossibleAssessments <- what_can_we_do$possible_assessments
        #           browser()
        #         }

        Output<- as.data.frame(matrix(NA,nrow=0,ncol=9))

        colnames(Output)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')


        Assessments[[timeStep]] <- Output
        #         Assessments <- Assess(Fishery = Fishery, Data = Monitor_Data, Assessments = PossibleAssessments, BatchFolder = BatchFolder, FigureFolder = FigureFolder)

        #==SARAH'S MORE COMPLICATED SAMPLING PROTOCOL FOR FUTURE USE==
        #Sample at the end of fishing if timeStep equals a sampling time step.
        #             if(timeStep %in% (sampleTimeSteps+burn)){
        #               Data <- CollectData(timeStep, simTime, burn, FleetN, DataParams, NoTakeZone, SpaceNumAtAgeT,
        #                                   SpaceCatAgeByFisher, SpaceEffort, wgtAtAge, lenAtAge, lenSD, PortX, PortY,kmax,Linf)
        #
        #               # Assign data to storage
        #               if(DataParams$Aggregate == 0){
        #                 if(timeStep == DataParams$SampStartYr+burn){
        #                   FDCatchData[DataParams$histStartYr:(timeStep-burn),,,] <- Data$histCatchDat
        #                   FDCPUEData[DataParams$histStartYr:(timeStep-burn),,,] <- Data$histCPUEDat
        #                   } else {
        #                   FDCatchData[(timeStep-burn),,,] <- Data$FDCatch
        #                   FDCPUEData[(timeStep-burn),,,] <- Data$FDCPUE
        #                   FICatchData[(timeStep-burn),,] <- Data$FICatch
        #                   FICPUEData[(timeStep-burn),,] <- Data$FICPUE
        #                   }
        #                 FDAgeData[(timeStep-burn),,,,] <- Data$AgeFD
        #                 FDSizeData[(timeStep-burn),,,,] <- Data$SizeFD
        #               } else {
        #                 if(timeStep == DataParams$SampStartYr+burn){
        #                   FDCatchData[DataParams$histStartYr:(timeStep-burn),] <- Data$histCatchDat
        #                   FDCPUEData[DataParams$histStartYr:(timeStep-burn),] <- Data$histCPUEDat
        #                 } else {
        #                   FDCatchData[(timeStep-burn),] <- Data$FDCatch
        #                   FDCPUEData[(timeStep-burn),] <- Data$FDCPUE
        #                   FICatchData[(timeStep-burn)] <- Data$FICatch
        #                   FICPUEData[(timeStep-burn)] <- Data$FICPUE
        #                 }
        #                 FDAgeData[(timeStep-burn),,] <- Data$AgeFD
        #                 FDSizeData[(timeStep-burn),,] <- Data$SizeFD
        #               }
        #
        #               }# if sample time steps
        #
        ######## Should add if statement here with a management loop####
        ####################################################
      } #end fleet

      #==natural mortality (depends on timestep size)
      tempSNALT<-tempSNALT*exp(-(M/yearMark))

      #==recruitment (happens once a year)
      #==calculate recruits based on total spawnign biomass, allocate them according to either the spawning biomass in a patch or evenly.
      #==current SpBio
      tempNum<-rep(0,kmax)
      for(i in 1:kmax)
        tempNum[i]<-sum(tempSNALT[,,i]*FishSel[i,flt])
      NumbersExplt[timeStep]<-sum(tempNum)

      tempSpBio<-rep(0,kmax)
      for(i in 1:kmax)
        tempSpBio[i]<-sum(tempSNALT[,,i]*MatAge[i]*wgtAtAge[i])
      SpawningBiomass[timeStep]<-sum(tempSpBio)

      tempSpBio<-rep(0,kmax)
      for(i in 1:kmax)
        tempSpBio[i]<-sum(-1*(NoTakeZoneImp-1)*tempSNALT[,,i]*MatAge[i]*wgtAtAge[i])
      InsideMPAspbio[timeStep]<-sum(tempSpBio)

      tempSpBio<-rep(0,kmax)
      for(i in 1:kmax)
        tempSpBio[i]<-sum(NoTakeZoneImp*tempSNALT[,,i]*MatAge[i]*wgtAtAge[i])
      OutsideMPAspbio[timeStep]<-sum(tempSpBio)

      if(timeStep%%yearMark == 0)
      {
        #==move to next year class
        tempSNALT[,,2:(kmax-1)]<-tempSNALT[,,1:(kmax-2)]
        tempSNALT[,,kmax]<-tempSNALT[,,kmax]+tempSNALT[,,kmax-1]
        #==recruitment
        Recruits<-Recruitment(SpawningBiomass[timeStep],steepness,R0,VirSpBio,detRec,sigmaR)
        if(RecDist==1) # equally spread recruitment
          tempSNALT[,,1]<-Recruits/(SpaceR*SpaceC)
        if(RecDist==2)
        {
          #==find proportion of total SpBio in each area
          tempSpBioMat<-matrix(0,nrow=SpaceR,ncol=SpaceC)
          for(i in 1:kmax)
            tempSpBioMat<-tempSpBioMat+tempSNALT[,,i]*MatAge[i]
          RecDistMat<-tempSpBioMat/sum(tempSpBioMat)
          tempSNALT[,,1]<-Recruits*RecDistMat
        }
        if(RecDist==3)
        {
          ##===INSERT habitat.csv DISTRIBUTION==
          RecDistMat<-habitat/sum(habitat)
          tempSNALT[,,1]<-(Recruits*as.matrix(RecDistMat))
        }
      }

      #==movement depends on the probability of moving at age

      temp1<-tempSNALT
      for(age in 1:kmax)
      {
        tempAge<-array(dim=c(SpaceR,SpaceC,SpaceR*SpaceC))
        Movers<-MoveProb[age]*tempSNALT[,,age]
        NonMovers<-(1-MoveProb[age])*tempSNALT[,,age]

        for(h in 1:nrow(coords))
          tempAge[,,h]<-Movers[coords[h,1],coords[h,2]]*Movement[,,h]

        Moved<-apply(tempAge,c(1,2),sum)
        tempSNALT[,,age]<-Moved+NonMovers
      }

      if(GraphsFish==T)
      {

        filled.contour(temp1[,,3],main="postfishery, premovement: age 2",zlim=c(0,2*OrigTempSNALT))
        mtext(side=3,paste("Timestep = ",timeStep))
        #       ani.record()
        filled.contour(tempSNALT[,,3],main="postfishery, postmovement: age 2",zlim=c(0,2*OrigTempSNALT))
        mtext(side=3,paste("Timestep = ",timeStep))
        #       ani.record()
      }
      #==update popdym
      if(timeStep<simTime)
        SpaceNumAtAgeT[timeStep+1,,,]<-tempSNALT


      #===========================================================================
      #== Calculate the cost of management for a given scenario
      #===========================================================================

      #==costs of management

      CostOfManagement[timeStep]<- CalculateManagementCosts(Fishery,timeStep,initManage,ManageStrats,Management)

      #       if(sum(1-NoTakeZone)>1 & timeStep==initManage)
      #         CostOfManagement[timeStep]<- CostOfManagement[timeStep] + MPAsunk
      #
      #       if((SizeLimit[flt]>0)& timeStep==initManage)
      #         CostOfManagement[timeStep]<- CostOfManagement[timeStep] + SizeSunk

      #==if there are any months in which fishing is not allowed, add sunk costs
      FleetSeason<-season[,2]
      FleetSeason[is.na(FleetSeason)]<-0
      if(FleetN>1)
        FleetSeason<-apply(!is.na(season[,2:ncol(season)]),1,sum)
      #       if(any(FleetSeason<FleetN,na.rm=T) & timeStep==initManage)
      #         CostOfManagement[timeStep]<- CostOfManagement[timeStep] + SeasonSunk
      #
      #       CostOfManagement[timeStep]<- CostOfManagement[timeStep] + sum(1-NoTakeZone)*MPAcost

      #       if((SizeLimit[flt])>0)
      #         CostOfManagement[timeStep]<- CostOfManagement[timeStep] + SizeCost
      #
      #==if all fisheries are open, there is no cost of enforcing a season
      #       if(FleetSeason[timeStep%%12+1]<FleetN)
      #         CostOfManagement[timeStep]<- CostOfManagement[timeStep] + SeasonCost

    } # end timestep


    Results[[s]]<- (list(CatchByFisher=CatchByFisher,CostByFisher=CostByFisher,ProfitByFisher=ProfitByFisher,CostOfManagement=CostOfManagement,
                         SpawningBiomass=SpawningBiomass,SpaceEffort=SpaceEffort,InsideMPAspbio=InsideMPAspbio,OutsideMPAspbio=OutsideMPAspbio,
                         ExploitableNumbers=NumbersExplt,Fishery=Fishery,IndexData=IndexData,LengthData=LengthData,CatchData=CatchData))
  } #Close stochasticity loop
  return(Results)
}


# outs<-Master(Life,SimCTL,Fleets,season,Samp,NoTakeZone,habitat,Graphs=F,GraphsFish=F,PrintLifeHistory=F)
#
# totCatch<-apply(CatchByFisher,2,sum,na.rm=T)
# totCost<-apply(CostByFisher,2,sum,na.rm=T)
# totProfit<-apply(ProfitByFisher,2,sum,na.rm=T)
#
# par(mfrow=c(4,1),mar=c(.1,4,.1,.1))
# plot(totCatch,type="b",xaxt='n',las=2)
# plot(totCost,lty=2,type="b",xaxt='n',las=2)
# plot(totProfit,lty=2,type="b",xaxt='n',las=2)
# lines(CostOfManagement,lty=2,type="b",col=2)
# plot(SpawningBiomass[burn:simTime],type="b",xaxt='n',las=2,ylab="SpawningBio",ylim=c(0,max(SpawningBiomass[burn:simTime],na.rm=T)))
#
# dim(CatchByFisher)
# dim(CostByFisher)
# dim(ProfitByFisher)
# length(CostOfManagement)
# length(SpawningBiomass)


