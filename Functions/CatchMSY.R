CatchMSY<- function(CatchDat,n,ProcessError,Smooth,Display,CatchJitters,CatchError,CatchBias,InterpCatch,StartYear,StartBio,MidYear,MidBio,FinalBio)
{
  ### Run CatchMSY, adapted form Martell & Froese 2012
  # CatchDat=CatchData
  # n=5000
  # ProcessError=0.05
  # Smooth=0
  # Display=1
  # CatchJitters=1
  # CatchError=0.075
  # CatchBias= 0 #the mean of the log normal errors in catch, 0 means no bias, above 0 positive bias, - negative bias
  # InterpCatch=1
  # StartYear=1994
  # StartBio=c(0.5,0.75)
  # MidYear=NA
  # MidBio=NA
  # FinalBio=c(0.25,0.5)
  # # Fish$res='Low'
  # CatchDat: Time series of catch
  # n: The number of iterations to search over
  # ErrorSize: The maginitude of the gradient of life history terms to search over
  # Smooth: Marks whether or not to smooth catch history
  # Display: Show diagnostics as it goes
  # CatchJitters: Number of times to try modifications of the catch
  # CatchError: Log error magnitude 
  
  
  CatchDat<- if (is.na(StartYear)){CatchDat} else {CatchDat[CatchDat$Year>=StartYear,]}
  
  yr   <- CatchDat$Year
  
  set.seed(999)  ## for same random sequence
  
  Output<- as.data.frame(matrix(NA,nrow=length(yr),ncol=9))
  
  colnames(Output)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')
  
  MCOutput<- as.data.frame(matrix(NA,nrow=length(yr),ncol=7))
  
  colnames(MCOutput)<- c('Iteration','Year','Method','SampleSize','Value','Metric','Flag')
  
  MCDetails<- as.data.frame(matrix(NA,nrow=length(yr)*CatchJitters,ncol=2))
  
  colnames(MCDetails)<- c('Year','FvFmsy')  
  
  
  ## filename <- "RAM_MSY.csv"
  # setwd('/Users/danovando/Desktop/Bren/SFG Work/Will Fish Recovery Increase Food')
  # filename <- "ICESct2.csv"
  # outfile  <- "CatchMSY_Output.csv"
  # outfile2  <- "Clean_CatchMSY_Output.csv"
  
  # FigureFolder<- 'Figures'
  # dir.create(FigureFolder)
  # cdat <- read.csv2(filename, header=T, dec=".")
  # cat("\n", "File", filename, "read successfully","\n")
  # output2<- NULL
  
  
  ## stock_id <- "cod-2224" ## for selecting individual stocks
  
  ## Loop through stocks
  for(c in 1:CatchJitters) {
    
    if (InterpCatch==0)
    {
      ct   <- CatchDat$Catch *rlnorm(length(yr),CatchBias,CatchError)
    }
    if (InterpCatch==1)
    {
      
      InterpCatch<- na.approx(CatchDat$Catch)
      
      ct=InterpCatch * rlnorm(length(yr),CatchBias,CatchError)
    }
    
    if(Smooth==1){ct<- runmed(ct,3)}
    
    res  <- Fish$res
    
    nyr  <- length(yr)    ## number of years in the time series
    
    ## PARAMETER SECTION
    
    ## If resilience is to be used, delete ## in rows 1-4 below and set ## in row 5	below
    for (i in 1)
    {
      start_r  <- if(res == "Very low"){c(0.015, 0.1)}
      else if(res == "Low") {c(0.05,0.5)}
      else if(res == "High") {c(0.6,1.5)}
      else {c(0.2,1)} ## Medium, or default if no res is found
    }	
    ## start_r     <- c(0.5,1.5)  ## disable this line if you use resilience
    start_k     <- c(max(ct),100*max(ct)) ## default for upper k e.g. 100 * max catch
    ## startbio 	<- c(0.8,1)   ## assumed biomass range at start of time series, as fraction of k
    startbio    <- if(ct[1]/max(ct) < 0.5) {c(0.5,0.9)} else {c(0.3,0.6)} ## use for batch processing 
    
    startbio<- if(is.na(StartBio[1])){startbio} else {StartBio}
    
    interyr 	<- if (is.na(MidYear[1])){yr[2]} else {MidYear}   ## interim year within time series for which biomass estimate is available; set to yr[2] if no estimates are available #SUB IN INTERMIN YEAR
    interbio 	<- if (is.na(MidBio[1])){c(0, 1)} else {MidBio} ## biomass range for interim year, as fraction of k; set to 0 and 1 if not available
    
    ## finalbio 	<- c(0.8, 0.9) ## biomass range after last catches, as fraction of k
    
    finalbio    <- if(ct[nyr]/max(ct) > 0.5) {c(0.3,0.7)} else {c(0.01,0.4)} ## use for batch processing #SET TO KNOWN B/BMSY RANGE
    
    finalbio<- if(is.na(FinalBio[1])){finalbio} else {FinalBio}
    
    sigR        <- ProcessError      ## process error; 0 if deterministic model; 0.05 reasonable value? 0.2 is too high
    
    startbt     <- seq(startbio[1], startbio[2], by = 0.05) ## apply range of start biomass in steps of 0.05	
    parbound <- list(r = start_r, k = start_k, lambda = finalbio, sigR)
    
    if (Display==1)
    {
      cat("Last year =",max(yr),", last catch =",ct[nyr],"\n")
      cat("Resilience =",res,"\n")
      cat("Process error =", sigR,"\n")
      cat("Assumed initial biomass (B/k) =", startbio[1],"-", startbio[2], " k","\n")
      cat("Assumed intermediate biomass (B/k) in", interyr, " =", interbio[1],"-",interbio[2]," k","\n")
      cat("Assumed final biomass (B/k) =", parbound$lambda[1],"-",parbound$lambda[2]," k","\n")
      cat("Initial bounds for r =", parbound$r[1], "-", parbound$r[2],"\n")
      cat("Initial bounds for k =", format(parbound$k[1], digits=3), "-", format(parbound$k[2],digits=3),"\n")
    }
    
    
    ## FUNCTIONS
    .schaefer	<- function(theta)
    {
      
      with(as.list(theta), {  ## for all combinations of ri & ki
        bt=vector()
        ell = 0  ## initialize ell
        for (j in startbt)
        {
          if(ell == 0) 
          {
            bt[1]=j*k*exp(rnorm(1,0, sigR))  ## set biomass in first year
            for(i in 1:nyr) ## for all years in the time series
            {
              xt=rnorm(1,0, sigR)
              bt[i+1]=(bt[i]+r*bt[i]*(1-bt[i]/k)-ct[i])*exp(xt) ## calculate biomass as function of previous year's biomass plus net production minus catch
            }
            
            #Bernoulli likelihood, assign 0 or 1 to each combination of r and k
            ell = 0
            if(bt[nyr+1]/k>=lam1 && bt[nyr+1]/k <=lam2 && min(bt) > 0 && max(bt) <=k && bt[which(yr==interyr)]/k>=interbio[1] && bt[which(yr==interyr)]/k<=interbio[2]) 
            {ell = 1
            }
            
          }
        }
        
        return(list(ell=ell,bio=bt))
        
        
      })
    }
    
    sraMSY	<-function(theta, N)
    {
      #This function conducts the stock reduction
      #analysis for N trials
      #args:
      #	theta - a list object containing:
      #		r (lower and upper bounds for r)
      #		k (lower and upper bounds for k)
      #		lambda (limits for current depletion)
      
      
      with(as.list(theta), 
           {
             ri = exp(runif(N, log(r[1]), log(r[2])))  ## get N values between r[1] and r[2], assign to ri
             ki = exp(runif(N, log(k[1]), log(k[2])))  ## get N values between k[1] and k[2], assing to ki
             itheta=cbind(r=ri,k=ki, lam1=lambda[1],lam2=lambda[2], sigR=sigR) ## assign ri, ki, and final biomass range to itheta
             M = apply(itheta,1,.schaefer) ## call Schaefer function with parameters in itheta
             
             i=1:N
             
             ## prototype objective function
             get.ell=function(i) M[[i]]$ell
             ell = sapply(i, get.ell) 
             
             get.bio=function(i) M[[i]]$bio
             bio = sapply(i, get.bio) 
             
             return(list(r=ri,k=ki, ell=ell,bio=bio))	
           })
    }
    
    ## MAIN
    R1 = sraMSY(parbound, n)  
    
    
    ## Get statistics on r, k, MSY and determine new bounds for r and k
    r1 	<- R1$r[R1$ell==1]
    k1 	<- R1$k[R1$ell==1]
    bio1<- R1$bio[,R1$ell==1]
    msy1  <- r1*k1/4
    
    mean_msy1 <- exp(mean(log(msy1))) 
    max_k1a  <- min(k1[r1<1.1*parbound$r[1]]) ## smallest k1 near initial lower bound of r
    max_k1b  <- max(k1[r1*k1/4<mean_msy1]) ## largest k1 that gives mean MSY
    max_k1 <- if(max_k1a < max_k1b) {max_k1a} else {max_k1b}
    
    if(length(r1)<10) {
      cat("Too few (", length(r1), ") possible r-k combinations, check input parameters","\n")
      flush.console()
    }
    
    if(length(r1)>=10) {
      
      ## set new upper bound of r to 1.2 max r1
      parbound$r[2] <- 1.2*max(r1)
      ## set new lower bound for k to 0.9 min k1 and upper bound to max_k1 
      parbound$k 	  <- c(0.9 * min(k1), max_k1)
      
      if (Display==1)
      {
        
        cat("First MSY =", format(1000*mean_msy1, digits=3),"\n")
        cat("First r =", format(exp(mean(log(r1))), digits=3),"\n")
        cat("New upper bound for r =", format(parbound$r[2],digits=2),"\n")	
        cat("New range for k =", format(1000*parbound$k[1], digits=3), "-", format(1000*parbound$k[2],digits=3),"\n")
      }
      ## Repeat analysis with new r-k bounds
      R1 = sraMSY(parbound, n)
      
      ## Get statistics on r, k and msy
      r = R1$r[R1$ell==1]
      k = R1$k[R1$ell==1]
      msy = r * k / 4
      mean_ln_msy = mean(log(msy))
      bio<- R1$bio[,R1$ell==1]
      
      BvK<- bio/(matrix(rep(k,nyr+1),nrow=dim(bio)[1],ncol=dim(bio)[2],byrow=T))
      BvBmsy<- bio/(matrix(rep(k,nyr+1),nrow=dim(bio)[1],ncol=dim(bio)[2],byrow=T)/2)
      ctt<- rep(as.matrix(ct),length(msy))
      CvMSY=matrix(ctt,nrow=length(ct),ncol=length(msy))/matrix(rep(msy,length(ct)),nrow=length(ct),length(msy),byrow=T)
      
      FvFmsy<- CvMSY/BvBmsy[1:nyr,]
      
      F<- FvFmsy*matrix(rep(r/2,nyr),nrow=dim(FvFmsy)[1],ncol=dim(FvFmsy)[2],byrow=T)
      
      pdf(file=paste(FigureFolder,'CatchMSY BvK.pdf',sep=''))
      BvKBox=boxplot(t(BvK),names=c(yr,yr[length(yr)]+1),las=0,ylab='B/K')
      dev.off()
      
      pdf(file=paste(FigureFolder,'CatchMSY BvBmsy.pdf',sep=''))
      BvBmsyBox=boxplot(t(BvBmsy),names=c(yr,yr[length(yr)]+1),las=0,ylab='B/Bmsy')
      dev.off()
      
      pdf(file=paste(FigureFolder,'CatchMSY CatchvMSY.pdf',sep=''))
      CvMSYBox=boxplot(t(CvMSY),names=c(yr),las=0,ylab='Catch/MSY')
      dev.off()
      
      pdf(file=paste(FigureFolder,'CatchMSY FvFmsy.pdf',sep=''))
      FvFmsyBox=boxplot(t(FvFmsy),names=c(yr),las=0,ylab='F/Fmsy',outline=FALSE)
      dev.off()
      
      pdf(file=paste(FigureFolder,'CatchMSY KobePlot.pdf',sep=''))
      plot(BvBmsyBox$stats[3,1:nyr],FvFmsyBox$stats[3,],type='b',xlab='BvBmsy',ylab='FvFmsy')
      text(BvBmsyBox$stats[3,c(1,nyr)],FvFmsyBox$stats[3,c(1,nyr)],labels=c(yr[1],yr[nyr]),pos=4)
      abline(h=1,v=1,lty=2)
      dev.off()
      
      pdf(file=paste(FigureFolder,'Catch MSY Plots.pdf',sep=''))
      ## plot MSY over catch data
      par(mfcol=c(2,3))
      plot(yr, ct, type="l", ylim = c(0, max(ct)), xlab = "Year", ylab = "Catch ")
      abline(h=exp(mean(log(msy))),col="red", lwd=2)
      abline(h=exp(mean_ln_msy - 2 * sd(log(msy))),col="red")
      abline(h=exp(mean_ln_msy + 2 * sd(log(msy))),col="red")
      
      hist(r, probability=T, xlim=c(0, 1.2 * max(r)), main = "")
      abline(v=exp(mean(log(r))),col="red",lwd=2)
      abline(v=exp(mean(log(r))-2*sd(log(r))),col="red")
      abline(v=exp(mean(log(r))+2*sd(log(r))),col="red")
      
      plot(r1, k1, xlim = start_r, ylim = start_k, xlab="r", ylab="k")
      
      hist(k, probability=T, xlim=c(0, 1.2 * max(k)), xlab="k", main = "")
      abline(v=exp(mean(log(k))),col="red", lwd=2)	
      abline(v=exp(mean(log(k))-2*sd(log(k))),col="red")
      abline(v=exp(mean(log(k))+2*sd(log(k))),col="red")
      
      
      plot(log(r), log(k),xlab="ln(r)",ylab="ln(k)")
      abline(v=mean(log(r)))
      abline(h=mean(log(k)))
      abline(mean(log(msy))+log(4),-1, col="red",lwd=2)
      abline(mean(log(msy))-2*sd(log(msy))+log(4),-1, col="red")
      abline(mean(log(msy))+2*sd(log(msy))+log(4),-1, col="red")
      
      hist(msy, probability=T, xlim=c(0, 1.2 * max(msy)), xlab=paste("MSY(",CatchDat$Units[1],')'),main = "")
      abline(v=exp(mean(log(msy))),col="red", lwd=2)
      abline(v=exp(mean_ln_msy - 2 * sd(log(msy))),col="red")
      abline(v=exp(mean_ln_msy + 2 * sd(log(msy))),col="red")
      dev.off()
      
      if (Display==1)
      {
        cat("Possible combinations = ", length(r),"\n")
        cat("geom. mean r =", format(exp(mean(log(r))),digits=3), "\n")
        cat("r +/- 2 SD =", format(exp(mean(log(r))-2*sd(log(r))),digits=3),"-",format(exp(mean(log(r))+2*sd(log(r))),digits=3), "\n")
        cat("geom. mean k =", format(1000*exp(mean(log(k))),digits=3), "\n")
        cat("k +/- 2 SD =", format(1000*exp(mean(log(k))-2*sd(log(k))),digits=3),"-",format(1000*exp(mean(log(k))+2*sd(log(k))),digits=3), "\n")
        cat("geom. mean MSY =", format(1000*exp(mean(log(msy))),digits=3),"\n")
        cat("MSY +/- 2 SD =", format(1000*exp(mean_ln_msy - 2 * sd(log(msy))),digits=3), "-", format(1000*exp(mean_ln_msy + 2 * sd(log(msy))),digits=3), "\n")
        
      }
      Output<- as.data.frame(matrix(NA,nrow=length(yr),ncol=9))
      
      colnames(Output)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')
      
      Output$Year<-yr
      Output$Method<-'CatchMSY'
      Output$SampleSize<-ct
      Output$Metric<- 'Catch/MSY'
      Output$Value<- ct/exp(mean_ln_msy)
      Output$UpperCI<- ct/quantile(msy,0.05)
      Output$LowerCI<- ct/quantile(msy,0.95)
      Output$SD<- sd(log(msy))
      
      
      MCDetails$Year<- yr
      
      MCDetails$FvFmsy<- FvFmsyBox$stats[3,]
      
      MCDetails$LowerFvFmsy<- FvFmsyBox$stats[1,]
      
      MCDetails$UpperFvFmsy<- FvFmsyBox$stats[5,]
      
      
      MCDetails$BvBmsy<- BvBmsyBox$stats[3,1:nyr]
      
      MCDetails$LowerBvBmsy<- BvBmsyBox$stats[1,1:nyr]
      
      MCDetails$UpperBvBmsy<- BvBmsyBox$stats[5,1:nyr]
      
      ## Write results into outfile, in append mode (no header in file, existing files will be continued)
      output = data.frame( sigR, startbio[1], startbio[2], interbio[1], interbio[2], finalbio[1], finalbio[2], min(yr), max(yr), res, max(ct), ct[1], ct[nyr], length(r), exp(mean(log(r))), sd(log(r)), min(r), quantile(r,0.05), quantile(r,0.25), median(r), quantile(r,0.75), quantile(r,0.95), max(r), exp(mean(log(k))), sd(log(k)), min(k), quantile(k, 0.05), quantile(k, 0.25), median(k), quantile(k, 0.75), quantile(k, 0.95), max(k), exp(mean(log(msy))), sd(log(msy)), min(msy), quantile(msy, 0.05), quantile(msy, 0.25), median(msy), quantile(msy, 0.75), quantile(msy, 0.95), max(msy)) 
      
      write.csv(output, file = paste(ResultFolder,'Raw CatchMSY Output.csv',sep=''), row.names = FALSE)
      
    }
  }  ## Catch Jitter loop, get next stock or exit
  
  return(list(Output=Output,Details=MCDetails))
  
}
