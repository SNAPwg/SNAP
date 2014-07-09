
 DBSRA<- function(CatchDat,DCAC.start.yr,DCAC.end.yr,delta.yr,DBSRA.OFL.yr,FMSYtoMRatio,SD.FMSYtoMratio,Delta,SD.Delta, DeltaLowerBound,DeltaUpperBound,BMSYtoB0ratio, SD.BMSYtoB0ratio,BMSYtoB0LowerBound, BMSYtoB0UpperBound ,InterpCatch, Niter)
{

# DCAC and DB-SRA
# Version 4d
# CatchDat<- CatchData
# DCAC.start.yr <- 1963 #start of the catch period
# DCAC.end.yr<- 2011 #end of the catch period
# delta.yr<- 2011 #Year that current depletion is fit to
# DBSRA.OFL.yr<- 2011 #Year to calculate DBSRA OFL outputs
# M.est<- .32 #Natural morality estiamte
# SD.lnM<- .5 #natural mortality error estimate
# FMSYtoMratio <- 0.8 #ratio of Fmsy to M
# SD.FMSYtoMratio<- .05
# Delta<- 0.6
# SD.Delta<- 0.1
# DeltaLowerBound<- 0.5
# DeltaUpperBound<- 0.9
# BMSYtoB0ratio<- 0.4
# SD.BMSYtoB0ratio<- 0.1
# BMSYtoB0LowerBound<- 0.2
# BMSYtoB0UpperBound<- 0.5
# InterpCatch<-1

start.enchilada <- Sys.time()

yr<- unique(CatchDat$Year)

Output<- as.data.frame(matrix(NA,nrow=length(yr),ncol=9))
  
  colnames(Output)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')
  
  # MCOutput<- as.data.frame(matrix(NA,nrow=length(yr),ncol=7))
  
  # colnames(MCOutput)<- c('Iteration','Year','Method','SampleSize','Value','Metric','Flag')
  


# number of CPUs available for optional parallel processing
NumCPUs <- 1

# turn on parallel processing if #CPUs > 1
do.parallel <- ifelse(NumCPUs > 1, TRUE, FALSE)

# toggle M correction term in production model
M.correction <- TRUE

# specify number of random draws for DCAC and delay-difference model
# Niter <- 500

############################
# LOAD DATA AND R PACKAGES #
############################

# load snowfall package (also requires snow)
require(snowfall)

# get life history data

parms.df<- as.data.frame(matrix(NA,nrow=1,ncol=24))
colnames(parms.df)<- c('run','group','sci.name','common.name','species.code','exclude','age.mat','DCAC.start.yr','DCAC.end.yr','delta.yr','DBSRA.OFL.yr','M.est','SD.lnM','FMSYtoMratio','SD.FMSYtoMratio','Delta','SD.Delta','DeltaLowerBound','DeltaUpperBound','BMSYtoB0ratio','SD.BMSYtoB0ratio','BMSYtoB0LowerBound','BMSYtoB0UpperBound','random.seed')

parms.df$run<- 1
parms.df$group<- Fish$Group
parms.df$sci.name<- Fish$SciName
parms.df$common.name<- Fish$CommName
parms.df$common.name<- Fish$CommName
parms.df$species.code<- 'SpnyLob'
parms.df$exclude<- 0
parms.df$age.mat<- round(AgeAtLength(Fish$Mat50,Fish, Fish$LengthError))
parms.df$DCAC.start.yr<- DCAC.start.yr
parms.df$DCAC.end.yr<- DCAC.end.yr
parms.df$delta.yr<- delta.yr
parms.df$DBSRA.OFL.yr <- DBSRA.OFL.yr
parms.df$M.est <- Fish$M
parms.df$SD.lnM <- Fish$MortalityError
parms.df$FMSYtoMratio <- FMSYtoMratio
parms.df$SD.FMSYtoMratio <- SD.FMSYtoMratio
parms.df$Delta <- Delta
parms.df$SD.Delta <- SD.Delta
parms.df$DeltaLowerBound <- DeltaLowerBound
parms.df$DeltaUpperBound <- DeltaUpperBound
parms.df$BMSYtoB0ratio <- BMSYtoB0ratio
parms.df$SD.BMSYtoB0ratio <- SD.BMSYtoB0ratio
parms.df$BMSYtoB0LowerBound <- BMSYtoB0LowerBound
parms.df$BMSYtoB0UpperBound <- BMSYtoB0UpperBound
parms.df$random.seed <- 4000


# vector of species codes
sp.vec <- as.character(parms.df$species.code)
N.spp  <- length(sp.vec)

# import landings data

TempCatch<- NULL
for (y in 1:length(yr))
{
	TempCatch[y]<- sum(CatchDat$Catch[CatchDat$Year==yr[y]])
}


lands.df<- data.frame(CatchDat$Year, sp.vec[1],TempCatch)
colnames(lands.df)<- c('year','sp.code','catch.mt')

	if (InterpCatch==1)
	{
		
		InterpCatch<- na.approx(lands.df$catch.mt)
			
		lands.df$catch.mt<- InterpCatch	
			
	}


# get first and last year of landings for each species and append to parms.df
first.yr <- numeric(N.spp)
last.yr <- numeric(N.spp)
for (i in 1:N.spp)
{
   first.yr[i] <- min(lands.df[lands.df$sp.code==sp.vec[i], "year"])
   last.yr[i]  <- max(lands.df[lands.df$sp.code==sp.vec[i], "year"])
}

parms.df <- cbind.data.frame(parms.df, first.yr, last.yr)
dimnames(parms.df)[[1]] <- sp.vec
rm(i, first.yr, last.yr)

# print warning if DCAC start year or end year is beyond landings 
if (any(parms.df$DCAC.start.yr < parms.df$first.yr)) print("DCAC start year is earlier than landings")
if (any(parms.df$DCAC.end.yr > parms.df$last.yr)) print("DCAC end year is later than landings")

# print warning if any species has fewer catch records than years in time series
for (i in sp.vec)
{
   if ( nrow(lands.df[lands.df$sp.code==i,]) != length(parms.df[i,"first.yr"]:parms.df[i,"last.yr"]) )
   {
       print(paste("Species", i, "has missing catch records."))
   }
}

# sort annual landings data frame by species, year
lands.df <- lands.df[with(lands.df, order(sp.code, year)), ]

# convert time series of catch to list object (convenient format for DB-SRA)
Catch.list <- as.list(1:N.spp)
names(Catch.list) <- sp.vec
for (i in 1:N.spp)
{
   Catch.list[[i]] <- lands.df[lands.df$sp.code==sp.vec[i], "catch.mt"]
   names(Catch.list[[i]]) <- lands.df[lands.df$sp.code==sp.vec[i], "year"]
}

# sum catch for species i between start and target years for DCAC
DCAC.SumC <- numeric(N.spp)
names(DCAC.SumC) <- sp.vec
for (i in 1:N.spp)
{
    catch.i <- lands.df[lands.df$sp.code==sp.vec[i] &
                       lands.df$year>=parms.df[sp.vec[i],"DCAC.start.yr"] &
                       lands.df$year<=parms.df[sp.vec[i],"DCAC.end.yr"], "catch.mt"]
   DCAC.SumC[i] <- sum(catch.i)
}

parms.df <- cbind.data.frame(parms.df, DCAC.SumC)
DCAC.Nyrs <- parms.df$DCAC.end.yr - parms.df$DCAC.start.yr + 1
parms.df <- cbind.data.frame(parms.df, DCAC.Nyrs, AvgCatch=(DCAC.SumC/DCAC.Nyrs))
rm(DCAC.SumC, DCAC.Nyrs, i, catch.i)

###############################################################################
# function to simulate from bounded beta distribution with
# m=mean and s=standard deviation;
# NOTE: the mean is the mean of the truncated dist, but the SD is from a standard (0,1) beta
# the user specifies the upper and lower values [a,b]

rbeta.ab <- function(n, m, s, a, b)
{
   # calculate mean of corresponding standard beta dist
   mu.std <- (m-a)/(b-a)
   
   # calculate parameters of std. beta with mean=mu.std and sd=s
   alpha <- (mu.std^2 - mu.std^3 - mu.std*s^2) / s^2
   beta  <- (mu.std - 2*mu.std^2 + mu.std^3 - s^2 + mu.std*s^2) / s^2
   
   # generate n draws from standard beta
   b.std <- rbeta(n, alpha, beta)
   
   # linear transformation from beta(0,1) to beta(a,b)
   b.out <- (b-a)*b.std + a
   
   return(b.out)
}
###############################################################################

# initialize list object for random draws
rand.list <- as.list(1:N.spp)
names(rand.list) <- sp.vec

# generate random draws from input distributions
for (i in 1:N.spp)
{
   M <- parms.df[sp.vec[i], "M.est"]
   SD.lnM <- parms.df[sp.vec[i], "SD.lnM"]
   FMSYtoMratio <- parms.df[sp.vec[i], "FMSYtoMratio"]
   SD.FMSYtoMratio <- parms.df[sp.vec[i], "SD.FMSYtoMratio"]
   Delta <- parms.df[sp.vec[i], "Delta"]
   SD.Delta <- parms.df[sp.vec[i], "SD.Delta"]
   DeltaLowerBound <- parms.df[sp.vec[i], "DeltaLowerBound"]
   DeltaUpperBound <- parms.df[sp.vec[i], "DeltaUpperBound"]
   BMSYtoB0ratio <- parms.df[sp.vec[i], "BMSYtoB0ratio"]
   SD.BMSYtoB0ratio <- parms.df[sp.vec[i], "SD.BMSYtoB0ratio"]
   BMSYtoB0LowerBound <- parms.df[sp.vec[i], "BMSYtoB0LowerBound"]
   BMSYtoB0UpperBound <- parms.df[sp.vec[i], "BMSYtoB0UpperBound"]
   myseed <- parms.df[sp.vec[i], "random.seed"] + 9167*(0:3)
   
   # lognormal distribution for natural mortality (M)
   set.seed(myseed[1])
   M.vec     <- rlnorm(Niter, meanlog=(log(M)-0.5*SD.lnM^2), sdlog=SD.lnM)
   
   # lognormal distribution for ratio of ratio FMSY/M
   set.seed(myseed[2])
   FMSYtoM.vec  <- rlnorm(Niter, meanlog=(log(FMSYtoMratio)-0.5*SD.FMSYtoMratio^2), sdlog=SD.FMSYtoMratio)
   
   # simulate bounded beta draws for delta
   set.seed(myseed[3])
   Delta.vec <- rbeta.ab(Niter, Delta, SD.Delta,
                                DeltaLowerBound, DeltaUpperBound)
   
   # simulate draws for ratio of BMSY/B0
   set.seed(myseed[4])
   BMSYtoB0.vec  <- rbeta.ab(Niter, BMSYtoB0ratio, SD.BMSYtoB0ratio,
                             BMSYtoB0LowerBound, BMSYtoB0UpperBound)
   
   rand.list[[i]] <- cbind(M.vec, FMSYtoM.vec, Delta.vec, BMSYtoB0.vec)
   
}

# stop if RNGs generarate NAs (happens when delta SD is large relative to mean)
if(any(unlist(lapply(rand.list, FUN=function(x) any(is.na(x))))))
{
   stop("NAs in random draws -- check distribution parameters")
}

###############################################################################
# DCAC function that accepts vector of parameters from input distributions
# parm.vec is a vector of length 4 with values for M, Fmsy/M, Depletion Delta, and Bmsy/B0
# call this function from apply(); faster than using for loop
DCAC.fun <- function(parm.vec, SumOfCatch, NumberOfYears)
{
   
   # assign names to elements of parameter vector (parm.vec) to clarify code
   M        <- parm.vec[1]
   FMSYtoM  <- parm.vec[2]
   Delta    <- parm.vec[3] # proportion that biomass is reduced relative to biomass in start year (NOT "depletion")
   BMSYtoB0 <- parm.vec[4]
   
   # calculate DCAC
   DCAC <- SumOfCatch / (NumberOfYears + (Delta/(BMSYtoB0 * FMSYtoM * M)))
   
   return(DCAC)
}
###############################################################################

# PARALLEL EXECUTION OF DCAC

# initialize "cluster" (use of multiple cores)
sfInit(parallel=do.parallel, cpus=NumCPUs)

# export variables & functions needed by slaves for parallel operation
sfExport(list=c("DCAC.fun","rand.list","parms.df","sp.vec"))

DCAC.start.time <- Sys.time() # record time to see how fast this works

# apply function DCAC.fun to elements of rand.list, parms.df$DCAC.SumC, and parms.df$DCAC.Nyrs
DCAC.list <- sfLapply(1:N.spp, function(i) apply(rand.list[[i]], MARGIN=1, FUN=DCAC.fun,
                                                 SumOfCatch=parms.df[sp.vec[i], "DCAC.SumC"],
                                                 NumberOfYears=parms.df[sp.vec[i], "DCAC.Nyrs"])
                     )

# how long did DCAC calculation take?
# print(Sys.time()-DCAC.start.time)

names(DCAC.list) <- sp.vec
DCAC.df <- as.data.frame(DCAC.list)
rm(DCAC.list)

# print summary (mean and quantiles) of DCAC distributions
DCAC.summary <- t(apply(DCAC.df, MARGIN=2, FUN=quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))
DCAC.summary <- cbind.data.frame(species=parms.df[, "common.name"],
                                 mean.DCAC=apply(DCAC.df, 2, mean),
                                 DCAC.summary)

# stop parallel cluster
sfStop()

# write DCAC output
# write.csv(DCAC.df, "DCAC_results.csv", row.names=FALSE)
# write.csv(DCAC.summary, "DCAC_summary.csv", row.names=TRUE)

####################################################################
# functions to get exponent in Fletcher parameterization of P-T model
# based on input value of BMSY/B0
get.n <- function(BMSYtoB0ratio)
{
   # define half-width of region around 1/e (Fox model) to randomly
   # assign a value that won't blow up but is reasonable approximation
   eps <- 1e-3

   # search for n when BMSY/B0 < B0/e
   if(BMSYtoB0ratio < (1/exp(1) - eps))
   {
      # lower bound of 0.05 for n sets minimum BMSY/B0 just below 5% (about 4.3%); plenty low
      n.out <- optimize(phi.fun, interval=c(0.05, 1-eps^2), phi.true=BMSYtoB0ratio)[[1]]
   }

   # approximate Fox model for values of BMSY/B0 near B0/e
   if( all( BMSYtoB0ratio >= (1/exp(1) - eps),
            BMSYtoB0ratio < (1/exp(1) + eps)))
   {
      n.out <- sample(c(1-eps, 1+eps), size=1)
   }

   # search for n when BMSY/B0 > B0/e
   if(BMSYtoB0ratio >= (1/exp(1) + eps))
   {
      # upper bound of 90 for n allows BMSY to occur at up to approx. 95% of B0
      n.out <- optimize(phi.fun, interval=c(1+eps^2, 90), phi.true=BMSYtoB0ratio)[[1]]
   }

   return(n.out)
}

# objective function to minimize error between input BMSY/B0 (phi.true) and estimate based on current n value
phi.fun <- function(n, phi.true)
{
   phi.out <- (1/n)^(1/(n-1))
   out <- (phi.true - phi.out)^2
   return(out)
}

# estimate lower and upper bounds for B0 for each species
# lower bound = average catch used in DCAC
# upper bound = 110% of sumCatch/min(Delta) (or twice the hake B0, whichever is smaller)
B0.low.vec <- numeric(N.spp)
B0.high.vec <- numeric(N.spp)
for (i in 1:N.spp)
{
   B0.low.vec[i] <- parms.df[sp.vec[i],"AvgCatch"]
   B0.high.vec[i] <- min(3000000, 1.1 * sum(Catch.list[[i]])/min(rand.list[[i]][,"Delta.vec"]))
}

# append B0 bounds to parms.df, so you can check if the optimization step hit the boundaries
parms.df <- cbind.data.frame(parms.df, B0.low=B0.low.vec, B0.high=B0.high.vec)
rm(B0.low.vec, B0.high.vec)

###############################################################################
# Pella-Tomlinson-Fletcher-Schaefer-MacCall model

Delay.Diff.fun <- function(i)
{
   i=1
   
   first.yr  <- parms.df[sp.vec[i],"first.yr"]
   last.yr   <- parms.df[sp.vec[i],"last.yr"]
   
   # "delta.yr" IS THE YEAR THAT WILL BE FIT TO THE CURRENT DEPLETION VALUE ###
   delta.yr <- parms.df[sp.vec[i],"delta.yr"]
   
   catch.vec <- Catch.list[[i]]
   Amat <- parms.df[sp.vec[i], "age.mat"]
   N.yrs <- last.yr - first.yr + 1
   
   # assume B0 is greater than average catch
   B0.low  <- parms.df[sp.vec[i],"B0.low"]
   
   # assume B0 is less than it would be if surplus production were zero,
   # and depletion equaled (B0 - sumCatch)/B0; i.e. some catch is from production
   B0.high <- parms.df[sp.vec[i],"B0.high"]

   # create dataframe to hold results for this species:
   # 4 input distributions, FMSY, EMSY, BMSY, MSY,
   # OFL vector, biomass (B) vector, production (P) vector,
   # Bjoin, Flag.NegBiomass, Flag.MissTarget, Flag.OptimWarn, Flag.NAs
   temp.mat <- matrix(NA, ncol = (N.yrs+1)*3 + 8, nrow=Niter) # project time series one year beyond last landing
   results.df <- cbind.data.frame(rand.list[[i]],temp.mat)
   rm(temp.mat)
   
   # give names to columns in the results data frame
   names(results.df) <- c(names(results.df)[1:4],
                          c("FMSY", "EMSY", "BMSY", "MSY"),
                          paste("OFL",first.yr:(last.yr+1),sep=""),		# note extra year (forecast 1 yr past landings)
                          paste("B",first.yr:(last.yr+1),sep=""),		# note extra year (forecast 1 yr past landings)
                          paste("P",first.yr:(last.yr+1),sep=""),		# note extra year (forecast 1 yr past landings)
                          "Bjoin", "Flag.NegBiomass", "Flag.NAs","Flag.OptimWarn")
   
   results.df[,"Flag.OptimWarn"] <- 0
   results.df[,"Flag.NAs"] <- 0
   
   # calculate distribution of FMSY (instantaneous rate) from input values
   results.df[,"FMSY"] <- results.df[,"M.vec"] * results.df[,"FMSYtoM.vec"]
   
   # calculate distribution of EMSY (exploitation rate) from input values
   results.df[,"EMSY"] <- (1 - exp(-(results.df[,"M.vec"] + results.df[,"FMSY"]))) *
                          ( results.df[,"FMSY"] / (results.df[,"M.vec"] + results.df[,"FMSY"]))
   
   for (j in 1:Niter)
   {  
   	# print(j)
      M        <- results.df[j, "M.vec"]
      FMSYtoM  <- results.df[j, "FMSYtoM.vec"]
      Delta    <- results.df[j, "Delta.vec"] # fraction biomass is reduced relative to biomass in start year (Delta is NOT "depletion")
      BMSYtoB0 <- results.df[j, "BMSYtoB0.vec"]
      FMSY     <- results.df[j, "FMSY"]
      EMSY     <- results.df[j, "EMSY"]
      
      if (any(Delta>=1, Delta<=0))
      {
         stop("Delta must be between 0 and 1 given assumptions of current production model")
      }
      Depletion <- 1-Delta

      n <- get.n(BMSYtoB0)
      g <- (n^(n/(n-1)))/(n-1)
      
      # create vector to hold biomass time series for this set of parameters
      B.vec <- numeric(N.yrs + 1)
      
      # create a vector to hold production time series for this set of parameters
      P.vec <- numeric(N.yrs + 1)
      
      # production functions that minimize difference between obs/pred depletion in target year
      # over a range of B0 values; parm.vec is current draw of c(M, FMSYtoM, Delta, BMSYtoB0, plus FMSY, EMSY)
      prod.fun <- function(B0, parm.vec, C.vec, Amat, detailed.output = TRUE)
      {
         M        <- parm.vec[1]
         FMSYtoM  <- parm.vec[2]
         Delta    <- parm.vec[3]
         BMSYtoB0 <- parm.vec[4]
         FMSY     <- parm.vec[5]
         EMSY     <- parm.vec[6]

         # calculate BMSY for current guess at B0 and current value of BMSYtoB0
         BMSY <- B0 * BMSYtoB0
         MSY <- BMSY * EMSY
         
         # calculate Bjoin and "s" (slope of production/biomass ratio) if BMSYtoB0 <= 0.5
         if (BMSYtoB0 < 0.3)
         {
            Bjoin <- B0 * (0.5*BMSYtoB0)
            s <- (1-n) * g * MSY * (Bjoin^(n-2)) * B0^(-n)
         }
         if (all(BMSYtoB0 >= 0.3, BMSYtoB0 <= 0.5))
         {
            Bjoin <- B0 * (0.75*BMSYtoB0 - 0.075)
            s <- (1-n) * g * MSY * (Bjoin^(n-2)) * B0^(-n)
         }
         
         if (BMSYtoB0 > 0.5)
         {
            Bjoin <- NA
         }
         
         # production function for Pella-Tomlinson-Fletcher model
         PTF.fun <- function(n, g, MSY, Blag, B0)
         {
            if (Blag > 0)
            {
               P <- g * MSY * (Blag/B0) - g * MSY * (Blag/B0)^n
            } else {
               P <- 0
            }
            return(P)
         }
   
         # production function for hybrid PTF-Schaefer model
         # that approximates BHSRR productivity when BMSY/B0 < 0.5
         hybrid.fun <- function(n, g, MSY, Blag, B0, Bjoin, s)
         {
            if (Blag > 0)
            {
               if (Blag >= Bjoin)
               {
                  P <- g * MSY * (Blag/B0) - g * MSY * (Blag/B0)^n
               } else {
                  P <-  Blag * ( ((g*MSY*(Bjoin/B0)-g*MSY*(Bjoin/B0)^n) / Bjoin) + s*(Blag - Bjoin) )
               }
            } else {
               P <- 0
            }
            return(P)
         }

         # assume B0 = B(t=1)
         B.vec[1] <- B0
         
         # production in first year is always zero (starting from B0)
         P.vec[1] <- 0
         
         # PROJECT FORWARD FROM CURRENT B0 USING PRODUCTION MODEL
         # two options: with and without M correction
         if (M.correction)
         {
            for (k in 2:(N.yrs+1))
	    {
	        # production is zero (at B0) until source of spawning output has been harvested
	        if (k <= Amat)
	        {
	           # remove production term (P) because production is based on unfished biomass until k>Amat,
	           # but include M correction and replace B.vec[k-Amat] with B0
	           B.vec[k] <- B.vec[k-1] + (1-exp(-M))*(B0-B.vec[k-1]) - C.vec[k-1]
	        } else {
	        
	        
	           # if BMSY/B0 > 0.5, then use P-T-F model
	           if (BMSYtoB0>0.5)
	           {
	              P.vec[k] <- PTF.fun(n, g, MSY, B.vec[k-Amat], B0)
	           } else {
	           
	              # if BMSY/B0 <= 0.5, use hybrid model
	              P.vec[k] <- hybrid.fun(n, g, MSY, B.vec[k-Amat], B0, Bjoin, s)
	           }
	           # production model with M correction term
	           B.vec[k] <- B.vec[k-1] + P.vec[k] + (1-exp(-M))*(B.vec[k-Amat]-B.vec[k-1]) - C.vec[k-1]
	        }
	    }
         } else { # loop without M correction term
            for (k in 2:(N.yrs+1))
	    {
	        # production is zero (at B0) until source of spawning output has been harvested
	        if (k <= Amat)
	        {
	           # remove production term (P) because production is based on unfished biomass until k>Amat
	           B.vec[k] <- B.vec[k-1] - C.vec[k-1]
	        } else {
	        
	           # if BMSY/B0 > 0.5, then use P-T-F model
	           if (BMSYtoB0>0.5)
	           {
	              P.vec[k] <- PTF.fun(n, g, MSY, B.vec[k-Amat], B0)
	           } else {
	           
	              # if BMSY/B0 <= 0.5, use hybrid model
	              P.vec[k] <- hybrid.fun(n, g, MSY, B.vec[k-Amat], B0, Bjoin, s)
	           }
	           # production model without M correction term
	           B.vec[k] <- B.vec[k-1] + P.vec[k] - C.vec[k-1]
	        }
	    }
	 }

          # add year labels to the time series vectors
          names(B.vec) <- names(P.vec) <- first.yr:(last.yr+1)
          names(C.vec) <- first.yr:last.yr
          
          Btgt.to.B0.ratio <- as.numeric(B.vec[as.character(delta.yr)]/B0)
          obj.fun <- as.numeric(abs(B.vec[as.character(delta.yr)] - Depletion*B0))
          
          if (detailed.output == FALSE)
          {
             return(as.numeric(obj.fun))
          }
          if (detailed.output == TRUE)
          {
             return(list(Species = sp.vec[i],
                         ObjectiveFunction = as.numeric(obj.fun),
                         B0.Bounds = c(B0.low, B0.high),
                         # set catch in last.yr+1 to zero (no effect)
                         TimeSeries = cbind.data.frame(B.vec, P.vec, C.vec=c(C.vec,0)), 
                         BMSY = BMSY,
                         MSY = MSY,
                         Bjoin = Bjoin))
          }
          
      } # end of prod.fun
      
      # global assignment (used in tryCatch)
      tag <<- 0
      
      tryCatch(B0.opt <- optimize(prod.fun,
                                  interval = c(B0.low, B0.high),
                                  ### the following are arguments passed to prod.fun()
                                  parm.vec = as.numeric(results.df[j,1:6]),
                                  C.vec = catch.vec,
                                  Amat = Amat,
                                  detailed.output = FALSE,
                                  tol = 0.0001),
                                  # if optimize() generates warning,
                                  # then tag this iteration and store result in "Flag.OptimWarn" column
                                  warning = function(warn) { tag <<- 1 },
                                  finally = B0.opt <- optimize(prod.fun,
                                                               interval = c(B0.low, B0.high),
                                                               ### the following are arguments passed to prod.fun()
                                                               parm.vec = as.numeric(results.df[j,1:6]),
                                                               C.vec = catch.vec,
                                                               Amat = Amat,
                                                               detailed.output = FALSE,
                                                               tol = 0.0001)
              )

      # generate model results as a list
      prod.out <- prod.fun(B0.opt[[1]],
                           parm.vec=as.numeric(results.df[j,1:6]),
                           catch.vec,
                           Amat,
                           detailed.output = TRUE)
      
      # store model results for this iteration in this species
      Flag.NegBiomass <- any(prod.out[["TimeSeries"]]$B.vec <= 0)
      Flag.NAs <- any(is.na(prod.out[["TimeSeries"]]))
      result.vec <- as.numeric( c(prod.out[["BMSY"]],
                                  prod.out[["MSY"]],
                                  results.df[j,"EMSY"] * prod.out[["TimeSeries"]][,"B.vec"], # time series of OFL (EMSY*B(t))
                                  prod.out[["TimeSeries"]]$B.vec, # time series of biomass
                                  prod.out[["TimeSeries"]]$P.vec, # time series of production
                                  prod.out[["Bjoin"]],
                                  Flag.NegBiomass,
                                  Flag.NAs,
                                  tag)
                              )
      results.df[j,7:ncol(results.df)] <- result.vec
      
   } # end of iteration loop
   
   # add results.df to list containing results for all species
   return(results.df)
   
} # end of species loop (function)

### PARALLEL EXECUTION OF DB-SRA

# initialize "cluster" (use of multiple cores)
sfInit(parallel=do.parallel, cpus=NumCPUs)

# export variables & functions needed by slaves for parallel operation
sfExport(list=c("Catch.list",
                "M.correction",
                "rand.list",
                "parms.df",
                "sp.vec",
                "Niter",
                "get.n",
                "phi.fun")
        )

options(warn=1) # print any warnings as they occur -- only works in sequential mode (1 cpu)

DelayDiff.start.time <- Sys.time()

# Load-balanced parallel computation of all species (requires snowfall package)
results.list <- sfClusterApplyLB(1:N.spp, Delay.Diff.fun)

# how long did it take?
# print(Sys.time()-DelayDiff.start.time)

# stop parallel cluster
sfStop()

names(results.list) <- sp.vec

N.yrs.vec <- parms.df$last.yr - parms.df$first.yr + 1		# need to redefine outside of delay diff fxn
names(N.yrs.vec) <- sp.vec

#### FLAG RUNS WITH OBVIOUS ERRORS ######

# create flag for runs in which final depletion doesn't match current draw
for (i in 1:N.spp)
{
      B.target <- results.list[[i]][ ,paste("B",parms.df[i,"delta.yr"],sep="")]
      B0 <- results.list[[i]][ ,paste("B",parms.df[i,"first.yr"],sep="")]
      Depletion <- 1 - results.list[[i]][ ,"Delta.vec"]
      Flag.Miss <- as.numeric( abs( B.target/B0 - Depletion ) > 0.005) # depletion must be within 0.5% of target
      # if biomass in target year is NA, call it a miss (usually due to oscillations)
      Flag.Miss[is.na(Flag.Miss)] <- 1
      results.list[[i]] <- cbind(results.list[[i]], Flag.MissTarget=Flag.Miss)
      rm(B.target, B0, Depletion, Flag.Miss)
}

# create flag that identifies any runs that hit the bounds for B0
# DEFINITION: "hitting" the bound means coming within 1 mt of upper OR 1% of lower bound
for (i in 1:N.spp)
{
   B0.low.tol <- 0.01*parms.df[i,"B0.low"] # B0 "hits" lower bound if it is within 1% of avg. catch
   B0.low.logical  <- abs(parms.df[i,"B0.low"]-results.list[[i]][,N.yrs.vec[i]+10]) < B0.low.tol
   # upper bound is "hit" if B0 is within 1 ton
   B0.high.logical <- abs(parms.df[i,"B0.high"]-results.list[[i]][,N.yrs.vec[i]+10])<1
   test.mat <- cbind(B0.low.logical, B0.high.logical)
   results.list[[i]] <- cbind(results.list[[i]], Flag.HitBounds=as.numeric(apply(test.mat, 1, any)))
}
rm(B0.low.logical, B0.high.logical, test.mat)

# create flag that identifies any runs with an error from optimize(), but no other errors,
# and also flag "good runs" for easy extraction into "results.good"
results.good <- as.list(1:N.spp)
for (i in 1:N.spp)
{
   # identify flagged runs from optimize()
   Optim.flagged  <- results.list[[i]][,"Flag.OptimWarn"]
   # identify runs with no other errors
   No.other.flags <- apply(results.list[[i]][,c("Flag.NegBiomass","Flag.MissTarget","Flag.NAs","Flag.HitBounds")],
                           MARGIN=1, sum)<1
   Good.run.vec <- as.numeric(No.other.flags)
   ## replace NAs with 1
   Good.run.vec[is.na(Good.run.vec)] <- 1
   results.list[[i]] <- cbind(results.list[[i]],
                              Good.Run=Good.run.vec,
                              Flag.OptimWarnOnly=as.numeric(Optim.flagged>0 & No.other.flags>0))
   results.good[[i]] <- results.list[[i]][Good.run.vec>0,]
}
names(results.good) <- sp.vec
rm(Optim.flagged, No.other.flags, Good.run.vec)

# RECORD RESULTS OF ALL ITERATIONS FROM PTFS MODEL AS .CSV FILES
# dir.create(paste(getwd(),"/model_output",sep=""))
# write separate .csv files for each species' PTFS results
# for (i in 1:N.spp)
# {
   # write.table(results.list[[i]], paste(getwd(),"/model_output/", sp.vec[i], "_all_results.csv", sep=""),
               # sep=",", row.names=FALSE)
# }

### save all trajectories that don't go negative, hit bounds, miss the target depletion level, or have NAs in timeseries
# dir.create(paste(getwd(),"/retained_model_output",sep=""))
# for (i in 1:N.spp)
# {
   # write.table(results.good[[i]], paste(getwd(),"/retained_model_output/", sp.vec[i], "_good_results.csv", sep=""),
               # sep=",", row.names=FALSE)
# }

# summarize error messages (negative biomass, hitting bounds of B0, missing target, etc.)
errors.list <- list()
errors.df   <- data.frame(t(rep(NA, 9)))
names(errors.df) <- c("Species","Niter","GoodRuns","Flag.NegBiomass","Flag.MissTarget",
                      "Flag.HitBounds","Flag.NAs","Flag.OptimWarn","Flag.OptimWarnOnly")
for(i in 1:N.spp)
{
   errors.list[[i]] <- results.list[[i]][,c("Good.Run","Flag.NegBiomass","Flag.MissTarget","Flag.HitBounds",
                                            "Flag.NAs","Flag.OptimWarn","Flag.OptimWarnOnly")]
   errors.df[i,1]   <- sp.vec[i]
   errors.df[i,2]   <- Niter
   errors.df[i,3:9] <- apply(errors.list[[i]], MARGIN=2, sum)
}
rm(errors.list)
# errors.df
# write.table(errors.df, paste(getwd(),"/error_summary.csv", sep=""), sep=",", row.names=F, quote=F)

# save parms.df as a comma-delimited text file
# write.table(parms.df, paste(getwd(),"/parms.csv", sep=""), sep=",", row.names=F, quote=F)



# CALC SUMMARY STATS FOR UNFISHED BIOMASS
# create list with B0 draws from each species
B0.yrs <- paste("B",parms.df[,"first.yr"],sep="")
B0.list <- as.list(1:N.spp)
for (i in 1:N.spp)
{
   B0.list[[i]] <- results.good[[sp.vec[i]]][,B0.yrs[i]]
}
names(B0.list) <- sp.vec

x.tmp <- matrix(NA, nrow=N.spp, ncol=6)
for (i in 1:N.spp)
{
   x.tmp[i,1] <- mean(B0.list[[i]])
   x.tmp[i,2:6] <- quantile(B0.list[[i]], probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
}

B0.stats <- cbind.data.frame(Common.Name = parms.df[,"common.name"], x.tmp)
names(B0.stats ) <- c("Common.Name", "B0.Mean", "2.5%", "25%", "50%", "75%", "97.5%")
# print(B0.stats)
# write.table(B0.stats, file=paste(getwd(),"/DB-SRA_B0_summary_stats.csv", sep=""), sep=",", row.names=F, quote=F)
rm(x.tmp)



# CALC SUMMARY STATS FOR OFL IN "DBSRA.OFL.yr" (doesn't have to be the same as "delta.yr")
# create list with OFL timeseries
OFL.list <- as.list(1:N.spp)
for (i in 1:N.spp)
{
   # use "grep" to identify column names that begin (\\b) with "OFL" (period after OFL is wildcard)
   OFL.list[[i]] <- results.good[[sp.vec[i]]][, grep("\\bOFL.", names(results.good[[sp.vec[i]]]))]
}
names(OFL.list) <- sp.vec



OFLStats<- as.data.frame(matrix(NA,nrow=length(yr),ncol=5))

BioStats<- as.data.frame(matrix(NA,nrow=length(yr),ncol=4))

bStats<- as.data.frame(matrix(NA,nrow=length(yr),ncol=4))

FlatResults<- as.data.frame(matrix(NA,nrow=0,ncol=14))

#Go through and squeeze results into usable format

Results<- results.good[[sp.vec]]

for (y in 1:length(yr))
{
	
	WhereYear<- grepl(yr[y],colnames(Results))
	
	
	YearData<- Results[,WhereYear]
		
	b<- YearData[,grepl('B',colnames(YearData))]/Results$BMSY
	
	YearData<- cbind(YearData,b)

	Quants<- as.data.frame(apply(YearData,2,quantile,probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))
	
	Means<- as.data.frame(apply(YearData,2,mean))

	SD<- as.data.frame(apply(YearData,2,sd))

	
	TempResults<- data.frame(yr[y],Results[,1:8],YearData,lands.df$catch.mt[y])
	
	colnames(TempResults)<- NA
	
	FlatResults<- rbind(FlatResults,TempResults)
	
	OFLMean<- grepl('OFL',rownames(Means))
	OFLStats[y,]<- data.frame(yr[y],Means[OFLMean,1],Quants[1,OFLMean],Quants[5,OFLMean],SD[OFLMean,1],check.names=F)
	
	BioMean<- grepl('B',rownames(Means))
	BioStats[y,]<- data.frame(yr[y],Means[BioMean,1],Quants[1, BioMean],Quants[5, BioMean],check.names=F)
	
	bMean<- grepl('b',rownames(Means))
	bStats[y,]<- data.frame(yr[y],Means[bMean,1],Quants[1, bMean],Quants[5, bMean],check.names=F)
	
	
		
}



colnames(FlatResults)<- c( 'Year',colnames(Results)[1:8],'OFL','B','P','BvBmsy','Catch')
colnames(OFLStats)<- c( 'Year','OFL','LowerQuantile','UpperQuantile','SD')
colnames(BioStats)<- c( 'Year','Biomass','LowerQuantile','UpperQuantile')
colnames(bStats)<- c( 'Year','BvBmsy','LowerQuantile','UpperQuantile')


pdf(file=paste(FigureFolder,'DBSRA Biomass Boxplots.pdf'))
boxplot(B~Year,data= FlatResults,outline=F)
title(xlab='Year',ylab='Biomass (mt)')
dev.off()

pdf(file=paste(FigureFolder,'DBSRA BvBmsy Boxplots.pdf'))
boxplot(BvBmsy~Year,data= FlatResults, outline=F)
title(xlab='Year',ylab='B/Bmsy')
dev.off()

pdf(file=paste(FigureFolder,'DBSRA OFL Boxplots.pdf'))
boxplot(OFL~Year,data= FlatResults, outline=F)
title(xlab='Year',ylab='OFL (mt)')
dev.off()

pdf(file=paste(FigureFolder,'DBSRA CatchvOFL Boxplots.pdf'))
boxplot((Catch/OFL)~Year,data= FlatResults, outline=F)
title(xlab='Year',ylab='Catch/OFL')
dev.off()

pdf(file=paste(FigureFolder,'DBSRA Parameter Posteriors.pdf'))
par(mfrow=c(4,2))
hist(FlatResults$M.vec,xlab='M',main=NULL)
hist(FlatResults$FMSYtoM.vec,xlab='FMSYtoM',main=NULL)
hist(FlatResults$BMSYtoB0.vec,xlab='BMSYtoB0',main=NULL)
hist(FlatResults$Delta.vec,xlab='Depletion',main=NULL)
hist(FlatResults$FMSY,xlab='Fmsy',main=NULL)
hist(FlatResults$EMSY,xlab='Emsy',main=NULL)
hist(FlatResults$BMSY,xlab='Bmsy',main=NULL)
hist(FlatResults$MSY,xlab='MSY',main=NULL)
dev.off()

PlotColors<- heat.colors(4,alpha=0.6)


pdf(file=paste(FigureFolder,'DBSRA Line Plots.pdf'))
par(mfrow=c(3,1))
plot(OFLStats$OFL~OFLStats$Year,type='l',bty='n',xlab='Year',ylab='OFL (mt)')
polygon(x=c(OFLStats$Year,rev(OFLStats$Year)),y=c(OFLStats$UpperQuantile,rev(OFLStats$LowerQuantile)),col= PlotColors[1],border=F)

plot(BioStats$Biomass ~ BioStats $Year,type='l',bty='n',xlab='Year',ylab='Biomass (mt)')
polygon(x=c(BioStats $Year,rev(BioStats $Year)),y=c(BioStats $UpperQuantile,rev(BioStats $LowerQuantile)),col= PlotColors[2],border=F)

plot(bStats$BvBmsy~ bStats$Year,type='l',bty='n',xlab='Year',ylab='B/Bmsy',ylim=c(0,2.5))
polygon(x=c(bStats$Year,rev(bStats$Year)),y=c(bStats$UpperQuantile,rev(bStats$LowerQuantile)),col= PlotColors[3],border=F)
abline(h=1,lty=2)
dev.off()

OFL.yrs <- paste("OFL",parms.df[,"DBSRA.OFL.yr"],sep="")
x.tmp <- matrix(NA, nrow=N.spp, ncol=6)
for (i in 1:N.spp)
{
   x.tmp[i,1] <- mean(OFL.list[[i]][,OFL.yrs[i]])
   x.tmp[i,2:6] <- quantile(OFL.list[[i]][,OFL.yrs[i]], probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
}



OFL.stats <- cbind.data.frame(Common.Name = parms.df[,"common.name"],
                      Delta.Year = parms.df[,"delta.yr"],
                      OFL.Year = parms.df[,"DBSRA.OFL.yr"],
                      x.tmp)
names(OFL.stats) <- c("Common.Name", "Delta.Year", "OFL.Year", "OFL.Mean", "2.5%", "25%", "50%", "75%", "97.5%")
# print(OFL.stats)
# write.table(OFL.stats, file=paste(getwd(),"/DB-SRA_OFL_summary_stats.csv", sep=""), sep=",", row.names=F, quote=F)
rm(x.tmp)


Flag<- NA

 Output<- data.frame(yr,'DBSRA', lands.df$catch.mt,OFLStats$OFL, OFLStats$LowerQuantile, OFLStats$UpperQuantile, OFLStats$SD,'OFL',Flag)

colnames(Output)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')

Output$Method<- 'DBSRA'

Output$Metric<- 'OFL'

return(list(Output=Output,Details=list(AllResults=FlatResults,BioSeries=BioStats,BvBmsySeries=bStats,DCAC= DCAC.summary)))
}

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


LBAR<-function(LengthDat,LagLength,Weight,IncludeMPA,ReserveYr,OutsideBoundYr,Iterations,BootStrap,LifeError,Lc=NA)
{
##################
  ###### LBAR ######
  ##################
  #Source: Based off of Ault et al. 1998
  #Summary: Calculate F based on mean length of catch, assumption of natural mortality
  
  ################
  #### Inputs ####
  ################
  
  # LengthData: The matrix of length observations
  # LagLength: The number of years of historic length data to lag when calculating Lc. A lag of 1 calculates an 
  #           independent Lc for each year. A lag that is the same number of years as there is data calculates a 
  #           a single Lc for all years. An intermediate lag functions as a moving average.
  # Weight: The weight assigned to historic length data. 1 weights all years equally.
  # IncludeMPA: 1= Include length data from MPA, 0=Exclude length data from MPA. If no data from MPA is included,
  #           M is calculated using the mean of 3 life history based methods.            
  # ReserveYr: Either year of reserve implementation or NA
  # OutsideBoundYr: Either the year mgmt changed outside reserve (for bounding), or use NA for unbounded Lbar method
  # Iterations: Number of iterations to run the model. 1 makes the model run with just the default data
  # BootStrap: 1 means data are bootstrapped, 0 means not
  # Lc: You can manually set an Lc (for example, a min size limit), or set to NA and it will be calculated for you.It
  #           can be calculated as a separate Lc for each year, a single Lc for all years, or a moving avg.

  
  ################
  #### Organize Data ####
  ################
  
  # LengthDat<- LengthData
  
  # LagLength<- 1
  
  # Weight<- 0.2
  
  # IncludeMPA<- 0
  
  # Iterations<- 1000
  
  # BootStrap<- 1
  
  # LifeError<- 1
  
  # ReserveYr <- 2011
  
  # OutsideBoundYr <- NA
  
  # Lc <- NA
  
  ManualLc<- Lc
  
  Years<- sort(unique(LengthDat$Year))
  
  SampleSize<- NA

  Output<- as.data.frame(matrix(NA,nrow=length(Years),ncol=9))
  
  colnames(Output)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')
  
  MCOutput<- as.data.frame(matrix(NA,nrow=length(Years),ncol=7))
  
  colnames(MCOutput)<- c('Iteration','Year','Method','SampleSize','Value','Metric','Flag')
  
  MCDetails<- as.data.frame(matrix(NA,nrow=length(Years)*Iterations,ncol=11))
  
  colnames(MCDetails)<- c('TotalMortality','FishingMortality','NaturalMortality','Lbar_Inside',
                          'Lbar_Outside','Llam_Inside','Llam_Outside','Lc','Iteration','Year','BeddFmsy')  
  
  
  LengthDat<- LengthDat[is.na(LengthDat$Length)==F,]
  
  if (IncludeMPA==0) #Remove MPA data if needed
  {
    LengthDat<- LengthDat[LengthDat$MPA==0,]
  }
  
  c<- 0
  
  BaseFish<- Fish #Set Fish to base life history parameters
  
  for (i in 1:Iterations) #Loop over monte carlo runs
  {
    
    if (dim(LengthDat)[1]>0) #If any data are available
    {
      
      if (i>1 & LifeError==1) #Apply life history uncertainty 
      {
        
        Fish<- BaseFish
        
        Fish<- ApplyLifeHistoryError()
        
      }
      
      for (y in 1:length(Years)) #Loop over years
      {
        
        c<- c+1
        
        Flag<- 'None'
        
        Lengths<- LengthDat$Length[LengthDat$Year==Years[y]]
        
        SampleSize[y]<- length(Lengths)
        
        dYr<- y #current year
        
        weight<- Weight
        
        TotYrs <- length(Years)
        
        if (dYr == 1){
          tempLaggedYears<- seq(from=1,to=LagLength,by=1) # Pull out lagged years
        } else if (dYr == length(Years)){
          tempLaggedYears<- seq(from=dYr-(LagLength-1),to=dYr,by=1) # Pull out lagged years
        } else {
          start <- floor(LagLength/2)
          end <- ceiling(LagLength/2)
          if((dYr-start)<1) {    #Calculate correction if this is below 1, add it to end.
            CorrNum1 <- 1 - (dYr-start)
            start <- start - CorrNum1
            end <- end + CorrNum1
          }
          if((dYr+(end-1))>length(Years)){
            CorrNum2 <- dYr+(end-1)-length(Years)
            start <- start + CorrNum1
            end <- end - CorrNum1
          } 
          tempLaggedYears<- seq(from=dYr-start,to=dYr+(end-1),by=1)
        }
        
        Lag<-t(as.matrix(tempLaggedYears[tempLaggedYears>0])) 
        
        weights<- as.numeric(weight^t(apply(Lag,1,rev))) #Weight assigned to each year
        
        LaggedYears<- Years[Lag]
        
        TempData<- LengthDat[LengthDat$Year %in% LaggedYears,] #Pull out all the length data in the lagged years
        
        if (i>1 & BootStrap==1) #Bootstrap data if desired
        {
          
          NumPoints<- 1:dim(TempData)[1]  
          
          BootSample<- sample(NumPoints,length(NumPoints),replace=T)
          
          TempData<- TempData[BootSample,]	
        }
        
        sizestore<- rep(NA,length(LaggedYears))

        ##Loop over lagged years to calculate an Lc that will work for all years
        if(is.na(ManualLc)){
        for (d in 1:length(LaggedYears)) #Loop over lagged years
        {
          
          YearlySizes<- sort(TempData$Length[TempData$Year==LaggedYears[d]])
          
          #Assigns Lc to be the mode
          x<-table(YearlySizes) #create histogram of lengths
          maxPos <- which.max(x) #find the mode of the catch distribution
          sizestore[d] <- as.numeric(names(x)[min(maxPos)]) #Store the size corresponding to the first mode
            
          }
        #Calculate weighted Lc past lagged years
        # print(sizestore)
    
        Lc<- round(sum(sizestore*weights)/sum(weights),10) 
        }
          
        # print(Lc)
        
        #Calculate Llam -- used for bounding.
        
        if(IncludeMPA == 0){
          LlamIn <- NA
        } else if (Years[y]-ReserveYr < 2){
          #Need at least 2 years since reserve went in to calculate M
          LlamIn <- NA
        } else {
          #Calculate the Expected size of fish protected since reserve.
          AgeLc <- AgeAtLength(Lc,Fish,Fish$LengthError)
          AgeBound <- AgeLc+(Years[y]-ReserveYr)
          LlamIn = LengthAtAge(AgeBound,Fish, Fish$LengthError)     
        }
        
        if(is.na(OutsideBoundYr)){
          LlamOut <- Fish$Linf
        } else {
          AgeLc <- AgeAtLength(Lc,Fish, Fish$LengthError)
          AgeBound <- AgeLc+(Years[y]-OutsideBoundYr)
          LlamOut = LengthAtAge(AgeBound,Fish, Fish$LengthError)
        }
        
        #Calculate mean weights 
        if (is.na(LlamIn)){
          LbarIn <- NA
        } else {
          MPAdat <- subset(TempData, MPA == 1 & Year == Years[y] & Length >= Lc & Length < LlamIn)
          LbarIn <- mean(MPAdat$Length,na.rm=TRUE)
        }
        
         if (is.na(OutsideBoundYr)){      #Should outside length data be bounded?        
          Outdat <- subset(TempData, MPA == 0 & Year == Years[y] & Length >= Lc)
          # Outdat <- subset(TempData, TempData[TempData$MPA == 0 & TempData$Year == Years[y] & TempData$Length >= Lc,])

          LbarOut <- mean(Outdat$Length,na.rm=TRUE)
         } else {
           Outdat <- subset(TempData, MPA == 0 & Year == Years[y] & Length >= Lc & Length > LlamOut)
           LbarOut <- mean(Outdat$Length,na.rm=TRUE)
         }
       # print(c("In",Lc,LbarIn,LlamIn,"Out",Lc, LbarOut,LlamOut))  
        
        #################
        #### Analyze Data ####
        ################
        
        Bounded.BH <- function(z,Linf,K,Lc,Llam,Lbar) #Bounded Beverton-Holt function (see Ault. et al. 1998)
        {((z*(Lc-Lbar)+K*(Linf-Lbar))/(z*(Llam-Lbar)+K*(Linf-Lbar)))-(((Linf-Llam)/(Linf-Lc))^(z/K))}
        
        #allresults=as.data.frame(matrix(NA,nrow=1,ncol=6)) 
        
        #colnames(allresults)=c('M', 'F','Z','Lbar','FvFmsy','adjust')
        EmpM <- FALSE
        
        if (!is.na(LlamIn)){
        #Test whether the M function has a root  
          test<- NULL #Test whether the Z function has a root
          s<- seq(0.0000001,10,.01)
          for (t in 1:length(s))
          {
            
            test[t]<- Bounded.BH(s[t],Fish$Linf,Fish$vbk,Lc,LlamIn,LbarIn)    
          }
          
          if ((test[1]*test[length(s)])<0){ #multiply first and last to check for opposite signs.
            EmpM <- TRUE
          }
        }
          
          if (EmpM == FALSE | IncludeMPA == 0 | length(LengthDat$MPA == 1) == 0 | Years[y]-ReserveYr < 2)
          {
            #If no MPA data, of no root, or told not to use MPA, or not enough time has passed since 
            #MPA implementation then use sum of LHI methods
            M<- mean(c((1.2*Fish$vbk),(4.22/Fish$MaxAge),(exp(1.71-1.084*log(Fish$MaxAge)))))
            
          } else {
            # Calculate Bounded to estimate M
            BOUND_M <-(uniroot(function(z) Bounded.BH(z,Fish$Linf,Fish$vbk,Lc,LlamIn,LbarIn), c(0.0000001,10)))$root  #M   
            M <- mean(c(BOUND_M,(1.2*Fish$vbk),(4.22/Fish$MaxAge),(exp(1.71-1.084*log(Fish$MaxAge)))))
          }
        
       
        
        test<- NULL #Test whether the Z function has a root
        s<- seq(0.0000001,10,.01)
        for (t in 1:length(s))
        {
          
          test[t]<- Bounded.BH(s[t],Fish$Linf,Fish$vbk,Lc, LlamOut,LbarOut)
         
        }
        
        if ((test[1]*test[length(s)])>=0) #If the function has no root (aka no possibly natural mortality)
        {
          Flag<- 'Uniroot failed, function has no zero root'
          
          MCOutput[c,]<- data.frame(i,Years[y],'LBAR',SampleSize[y],-999,'FvM',Flag)
          
        }
          #If function has a root calcualte Z/M using bounded BH function
        
        if ((test[1]*test[length(s)])<0) 
        {
          
          
          BOUND_Z <-(uniroot(function(z) Bounded.BH(z,Fish$Linf,Fish$vbk,Lc,LlamOut,LbarOut), c(0.0000001,10)))$root  #M	 
          
          #Test Ault's numbert
          # BOUND_Z <-(uniroot(function(z) Bounded.BH(z,939,0.13,400,798,493), c(0.0000001,10)))$root  #M	 
          
          
          if (M > BOUND_Z) #If the result doesn't seem possible, calculate M from life history invariants
          {
            
            Flag<- 'Real M greater than Z, using Fish$M estimate'
            M<- Fish$M
            
          }
          
          
          BOUND_F <-BOUND_Z-M #Fishing mortality
          if (BOUND_F<0)
          {
            Flag<- 'LBAR Failed, Estimated Z too low'
          }
          
          FvFmsy <- BOUND_F/M #Calculate F/Fmsy
          
          ################
          #### Store Results ####
          ################
          
          MCOutput[c,]<- c(i,Years[y],'LBAR',SampleSize[y],FvFmsy,'FvM',Flag)
          
          
          MCDetails$TotalMortality[c]<- BOUND_Z
          
          MCDetails$FishingMortality[c]<- BOUND_F
          
          MCDetails$NaturalMortality[c]<- M
          
          MCDetails$Lbar_Inside[c]<- LbarIn
          
          MCDetails$Lbar_Outside[c]<- LbarOut
          
          MCDetails$Llam_Inside[c]<- LlamIn
          
          MCDetails$Llam_Outside[c]<- LlamOut
          
          MCDetails$Lc[c]<- Lc
          
          MCDetails$Iteration[c]<- i
          
          MCDetails$Year[c]<- Years[y]
          
          BeddFmsy<- (0.6*Fish$vbk)/(0.67-Lc/Fish$Linf)
          
          MCDetails$BeddFmsy[c]<- BeddFmsy

          
        } #Close if there's a root
        
        
      } #Close Year Loop	
      
    }
    if (dim(LengthDat)[1]==0) #If there are no data
      
    {
      Output<- c(Years,'LBAR',NaN,'FvM',NaN,NaN,NaN,'No Usable Data')
    }
    
  }
  
  TrueIteration<- MCOutput$Iteration==1
  
  TrueOutput<- MCOutput[TrueIteration,]
  
  # MCOutput<- MCOutput[TrueIteration==F,]
  
  Output$Year<- Years
  
  Output$Method<- 'LBAR'
  
  Output$Value<- TrueOutput$Value
  
  Output$LowerCI<- NaN
  
  Output$UpperCI<- NaN
  
  Output$SD<- NaN
  
  Output$Metric<- 'FvM'
  
  Output$Flag<- TrueOutput$Flag
  
  Output$SampleSize<- SampleSize
  
  if (Iterations>1) # If there are multiple iterations, process Monte Carlo output
  {
    
    for (y in 1:length(Years))
    {
      Where<- MCOutput$Year==Years[y]
      
      Temp<- MCOutput[Where,]
      
      TempValue<- sort(as.numeric(Temp$Value))
      
      Bottom<- ceiling(.025*length(TempValue))
      
      Top<- ceiling(.975*length(TempValue))
      
      MeanMetric<- mean(as.numeric(Temp$Value),na.rm=T)
      
      LowerCI<- TempValue[Bottom]
      
      UpperCI<- TempValue[Top]
      
      SD<- sd(TempValue[Bottom:Top],na.rm=T)
      
      Output$Year<- Years
      
      Output$Method[y]<- 'LBAR'
      
      Output$Value[y]<- TrueOutput$Value[y]
      
      Output$LowerCI[y]<- LowerCI
      
      Output$UpperCI[y]<- UpperCI
      
      Output$SD[y]<- SD
      
      Output$Metric[y]<- 'FvM'
      
    }
  }
  
  ################
  #### Produce Figures ####
  ################
  
  NumNans<-(is.na(MCDetails)[,1])
  if (sum(NumNans)<length(NumNans))
  {
    
    pdf(file=paste(FigureFolder,' LBAR FvM Boxplots.pdf',sep=''))
    boxplot((MCDetails$FishingMortality/M)~MCDetails$Year,frame=F,xlab='Year',ylab='F/M',notch=T,outline=F,width=SampleSize)
    dev.off()
    
    pdf(file=paste(FigureFolder,' LBAR Fishing Mortality Boxplots.pdf',sep=''))
    boxplot((MCDetails$FishingMortality)~MCDetails$Year,frame=F,xlab='Year',ylab='F',notch=T,outline=F,width=SampleSize)
    dev.off()
    
    pdf(file=paste(FigureFolder,' LBAR Total Mortality Boxplots.pdf',sep=''))
    boxplot((MCDetails$TotalMortality)~MCDetails$Year,frame=F,xlab='Year',ylab='Z',notch=T,outline=F, width=SampleSize)
    dev.off()
    
    pdf(file=paste(FigureFolder,' LBAR Mean Length Boxplots.pdf',sep=''))
    if(IncludeMPA == 1 & sum(is.na(MCDetails$Lbar_Inside))<length(MCDetails$Lbar_Inside)){
    par(mfrow=c(1,2))
    boxplot((MCDetails$Lbar_Inside)~MCDetails$Year,frame=F,xlab='Year',ylab='Mean Length',notch=T,main = "Mean Lengths Inside",outline=F, width=SampleSize)
    boxplot((MCDetails$Lbar_Outside)~MCDetails$Year,frame=F,xlab='Year',ylab='Mean Length',notch=T,main = "Mean Lengths Outside",outline=F, width=SampleSize)
    } else {
      boxplot((MCDetails$Lbar_Outside)~MCDetails$Year,frame=F,xlab='Year',ylab='Mean Length',notch=T,main = "Mean Lengths Outside",outline=F, width=SampleSize)
    }
    dev.off()
  }
  
  Fish<- BaseFish
 
  return(list(Output=Output,Details=MCDetails))	
}

   CatchCurve<- function(LengthDat,CatchCurveWeight,WeightedRegression, ReserveYr,OutsideBoundYr,ManualM,Iterations,BootStrap,LifeError,HistInterval)
  {
  ##################
  ###### CatchCurve ######
  ##################
  #Source:Loosely based on Kay & Wilson 2012, and Sparre & Venema 1998
  #Summary: Estiamte natural and total mortality from slope of age frequency histograms inside and outside of marine reserves
  
  ################
  #### Inputs ####
  ################
  
  
  #LengthDat: Length Data
  #CatchCurveWeight: The weight assigned to the catch curve relative to other methods of calculating mortality.'AgeBased' slides catch curve weight in proportion to the age of the MPA
  #ReserveYr: The year of MPA implementation
  #OutsideBoundYr: A year marker for any selectivity changes outside the MPA
  #ManualM: 1 or 0, 1 sets natural mortality to lit/LHI based value, 0 trys to estimate from MPA
  #Iterations: The number of iterations to run (1 runs only with default data)
  #BootStrap: 1 means bootstrap data in monte carlo runs
  #LifeError: 1 means that that life history has monte carlo based error introduced
  #HistInterval: The bin width of the size histograms
  
     # Fish<- BaseFish
     # Iterations<- 10
     # CatchCurveWeight<- 'AgeBased'
     # LengthDat<- LengthData
     # BootStrap<- 1
     # LifeError<- 1
     # HistInterval<- 1
	 # ManualM<- 1
   	# WeightedRegression<- 1 #marks whether to use a weighted regression
   # ReserveYr<- 2011
    # OutsideBoundYr<- NA
  
	SampleSize<- NA
  
   # LengthDat<- LengthDat[LengthDat$MPA==0,]
  
  ################
  #### Process Data ####
  ################
  
  CCWeight<- CatchCurveWeight
  
  Years<- sort(unique(LengthDat$Year))

  MCOutput<- as.data.frame(matrix(NA,nrow=length(Years)*Iterations,ncol=7))
  
  colnames(MCOutput)<- c('Iteration','Year','Method','SampleSize','Value','Metric','Flag')
  
  MCDetails<- as.data.frame(matrix(NA,nrow=length(Years)*Iterations,ncol=4))
  
  colnames(MCDetails)<- c('Iteration','Year','FishingMortality','NaturalMortality')
  
  Output<- as.data.frame(matrix(NA,nrow=length(Years),ncol=9))
  
  colnames(Output)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')
  
  Details<- as.data.frame(matrix(NA,nrow=length(Years)*Iterations,ncol=9))
  
  colnames(Details)<- c('Year','FishingMortality','LowerCI','UpperCI','SD','NaturalMortality','LowerCI','UpperCI','SD')
  
  BaseFish<- Fish
  
  c<- 0
  for (i in 1:Iterations) #Loop over Monte Carlo run
  {
    
    if (i>1 & LifeError==1) #Apply life history errors
    {
      Fish<- BaseFish
      
      Fish<- ApplyLifeHistoryError()
    }
    
    
    for (y in 1:length(Years)) #Loop over years
    {

      c<- c+1
      
      Flag<- 'None'
      
      TempLengthDat<- LengthDat[LengthDat$Year==Years[y],]
      
      TempLengthDat<- TempLengthDat[is.na(TempLengthDat$Length)==F,]
      
      SampleSize[y]<- dim(TempLengthDat)[1]
      
      if (i>1 & BootStrap==1) #Resample length data inside and outside of MPAs
      {
        
        FishedDat<- TempLengthDat[TempLengthDat$MPA==0,]

        NumPoints<- 1:dim(FishedDat)[1]	
        
        BootSample<- sample(NumPoints,length(NumPoints),replace=T)
        
        TempFishedDat<- FishedDat[BootSample,]
      
       	 MPADat<- TempLengthDat[TempLengthDat$MPA==1,]
        
        NumPoints<- 1:dim(MPADat)[1]	
  		
  		BootSample<- sample(NumPoints,length(NumPoints),replace=T)
        
        TempMPADat<- MPADat[BootSample,]

       TempLengthDat<- rbind(TempFishedDat,TempMPADat)

       
        
      }
      
      # TempLengthDat$Length[TempLengthDat$Length>=Fish$Linf]<- Fish$Linf*.98
      
      TempLengthDat<- TempLengthDat[TempLengthDat$Length<Fish$Linf,] #Only use length data less than Linf (might want to fix this)
      
      TempLengthDat$Age<- pmax(.01,AgeAtLength(TempLengthDat$Length,Fish,Fish$LengthError	)) #Calculate age at length

      BinBreaks<- seq(from=0,to=ceiling(max(TempLengthDat$Age,na.rm=T)+1),by=HistInterval)

      AllHist<- DanHist(TempLengthDat$Age,BinBreaks)
      
      TotalPeak<- which(AllHist$Frequency==max(AllHist$Frequency))[1]
      
      # BinBreaks<- 0:(Fish$MaxAge)

      # BinBreaks<- seq(from=0,to=(Fish$MaxAge),by=HistInterval)

      AgeDistFished<- DanHist(TempLengthDat$Age[TempLengthDat$MPA==0],BinBreaks) #Create histogram of age data outside MPA
      
      AgeDistFished$MPA<- 'Fished Area'
      
      FishedPeak<- which(AgeDistFished$Frequency ==max(AgeDistFished$Frequency))[1] #Identidy the mode of the fished age hist.

      # FishedPeak<- TotalPeak
      
      FishedAllObserved<- (which(AgeDistFished$Frequency>0))
          
      AgeDistMPA <- DanHist(TempLengthDat$Age[TempLengthDat$MPA==1],BinBreaks) #Create histogram of age data inside MPA
      
      AgeDistMPA$MPA<- 'MPA'
      
      MPAPeak<- which(AgeDistMPA$Frequency == max(AgeDistMPA$Frequency))[1] #Identify the mode of the MPA age hist.

      # MPAPeak<- TotalPeak
      
      MPAAllObserved<- (which(AgeDistMPA$Frequency>0))
      

      # if ((length(MPAAllObserved)>1 & length(FishedAllObserved)>1) & (is.na(AgeDistMPA$LogFrequency[MPAPeak])==F & is.na(AgeDistFished$LogFrequency[FishedPeak])==F)) #If you've got any real histogram to work with

      if (length(FishedAllObserved)>1 & is.na(AgeDistFished$LogFrequency[FishedPeak])==F) #If you've got any real histogram to work with

      {
      
      if (sum(TempLengthDat$MPA,na.rm=T)>0)
      {
      MPALastObserved<- MPAAllObserved[length(MPAAllObserved)-1] #Find second to last observed age in the MPA
      
        
        if (WeightedRegression==1)
        {
        	MPALastObserved<- MPALastObserved+1
        }
        
        if (is.na(ReserveYr)==F & Years[y]>ReserveYr) # Bound the last observed age by species young enough to only be impacted by the MPA
        {
        	BoundedLastObserved<- MPAPeak+(Years[y]-ReserveYr)
        	
        	if (BoundedLastObserved<MPALastObserved)
        	{
        		MPALastObserved<- BoundedLastObserved #Adjsut last observed age by MPA age
        	}	
        }
      
      
      if (MPAPeak==MPALastObserved & length(MPAAllObserved)!=1) #Stupid fix
      {
        MPAPeak<- MPAAllObserved[length(MPAAllObserved)-1]
      }
      
      if (length(MPAAllObserved)==1)
      {
      	Flag<- 'Catch Curve cannot run - Only 1 Age group observed'
      }
      
      NumAgeGroups<- length(AgeDistMPA$Frequency)
      
      MPANumPoints<- length(MPAPeak:MPALastObserved)
      
      MPACatchCurve<- lm((AgeDistMPA$LogFrequency[MPAPeak:MPALastObserved]) ~ (AgeDistMPA$Age[MPAPeak:MPALastObserved]),na.action='na.omit') #Fit catch curve between peak and final point
      
      PredictedMPAValues<- predict(MPACatchCurve,data.frame(Ages=seq(from=BinBreaks[MPAPeak],to=BinBreaks[MPALastObserved],length.out=MPANumPoints)))	
      
           if (WeightedRegression==1)
        {
        	
        	    RegWeights<- pmax(0,PredictedMPAValues)/sum(pmax(0,PredictedMPAValues))
        	      
        	      MPACatchCurve<- lm((AgeDistMPA$LogFrequency[MPAPeak:MPALastObserved]) ~ (AgeDistMPA$Age[MPAPeak:MPALastObserved]),na.action='na.omit',weights= RegWeights) #Fit catch curve between peak and final point
      
      PredictedMPAValues<- predict(MPACatchCurve,data.frame(Ages=seq(from=BinBreaks[MPAPeak],to=BinBreaks[MPALastObserved],length.out=MPANumPoints)))	

        }
        
      
      CatchCurveNaturalMortality<- -MPACatchCurve$coefficients[2] #Natural mortality: slope of MPA catch curve

      
      if (CCWeight=='AgeBased') #Allows the weight of the MPA based M to increase the older the MPA is
      {
      	CatchCurveWeight<- Years[y]-ReserveYr
      }
      
      MortalityWeight<- c(CatchCurveWeight,1,1,1)
      
      MortalityWeight<- MortalityWeight/sum(MortalityWeight)
      if(is.na(CatchCurveNaturalMortality))
      {
      	Flag<- 'Catch Curve could not estimate M'	
      }
      NaturalMortality<-  sum(MortalityWeight*c(CatchCurveNaturalMortality,(1.2*Fish$vbk),(3/Fish$MaxAge),(exp(1.71-1.084*log(Fish$MaxAge)))),na.rm=T)/sum(MortalityWeight) #Calcualte weighted natural mortality
      }
  
      if (ManualM==1 | sum(TempLengthDat$MPA,na.rm=T)==0 | (is.na(ReserveYr)==F & (Years[y]-ReserveYr)<2)) #Use lit or LHI based M if set that way, if there is no MPA, or if the MPA is less than 2 years old
      {
      
      MortalityWeight<- c(1,1,1,1)
      
      MortalityWeight<- MortalityWeight/sum(MortalityWeight)

      NaturalMortality<-  sum(MortalityWeight*c(Fish$M,(1.2*Fish$vbk),(3/Fish$MaxAge),(exp(1.71-1.084*log(Fish$MaxAge)))),na.rm=T)/sum(MortalityWeight) #Calcualte weighted natural mortality
      # show(NaturalMortality)

      	Flag<- 'No MPA based M possible - derived from LHI and Lit'
      }
 
      FishedLastObserved<- FishedAllObserved[length(FishedAllObserved)-1]
                
                 if (WeightedRegression==1)
        {
        	
        	  FishedLastObserved <-  FishedLastObserved+1

        }
        
      
      
        if (is.na(OutsideBoundYr)==F & Years[y]>OutsideBoundYr) # Allows bounding if there was a selectivity intervention outside the MPA
        {
        	BoundedLastObserved<- FishedPeak+(Years[y]-OutsideBoundYr)
        	
        	if (BoundedLastObserved<FishedLastObserved)
        	{
        		FishedLastObserved<- BoundedLastObserved
        	}	
        }
      
      FishedNumPoints<- length(FishedPeak:FishedLastObserved)

      if (FishedPeak==FishedLastObserved & length(FishedAllObserved)!=1  )
      {

		Flag<- 'Catch Curve cannot run, no righthand side of age distribution'
        # FishedPeak<- FishedAllObserved[length(FishedAllObserved)-1]
        
      }
      
      if (length(FishedAllObserved)==1)
      {
      	Flag<- 'Catch Curve cannot run - Only 1 Age group observed'
      }
      
      if (FishedPeak<FishedLastObserved)
      {
      FishedNumPoints<- length(FishedPeak:FishedLastObserved)
      
      FishedCatchCurve<- lm(AgeDistFished$LogFrequency[FishedPeak:FishedLastObserved] ~ AgeDistFished$Age [FishedPeak:FishedLastObserved], na.action='na.omit')#Fit catch curve between peak and final point
      
      PredictedFishedValues<- predict(FishedCatchCurve,data.frame(Ages=seq(from=BinBreaks[FishedPeak],to=BinBreaks[FishedLastObserved],length.out=FishedNumPoints)))	
      
                 if (WeightedRegression==1)
        {
        	
        	    RegWeights<- pmax(0, PredictedFishedValues)/sum(pmax(0, PredictedFishedValues))
        	      
      FishedCatchCurve<- lm(AgeDistFished$LogFrequency[FishedPeak:FishedLastObserved] ~ AgeDistFished$Age [FishedPeak:FishedLastObserved], na.action='na.omit',weights=RegWeights)#Fit catch curve between peak and final point
      
      PredictedFishedValues<- predict(FishedCatchCurve,data.frame(Ages=seq(from=BinBreaks[FishedPeak],to=BinBreaks[FishedLastObserved],length.out=FishedNumPoints)))	

        }
      
      
      
      CatchCurveTotalMortality<- -FishedCatchCurve$coefficients[2]
      }
      else
      {
      	CatchCurveTotalMortality<- NA
      }
      # TotalMortalityWeights<- c(CatchCurveWeight,1)
      
      # TotalMortalityWeights<-    TotalMortalityWeights/sum(TotalMortalityWeights)
      
      TotalMortality<- CatchCurveTotalMortality
      
      # TotalMortality<- sum(TotalMortalityWeights* c(CatchCurveTotalMortality,LBAR(TempLengthDat,1,1,0)$Details$TotalMortality))/sum(TotalMortalityWeights)
      
      FishingMortality<- TotalMortality-NaturalMortality
      
      
      if (FishingMortality<0 & is.na(FishingMortality)==F)
      {
        Flag<- 'Catch-Curve not working-Negative Fishing Mortality'
      }
 

      AgeDist<- rbind(AgeDistMPA,AgeDistFished)
      
      if (i==1) #Plot catch curves analysis on first iteration
      {
       
       
       if (sum(TempLengthDat$MPA,na.rm=T)>0)
       {
        pdf(file=paste(FigureFolder,i,Years[y],' Catch Curve Analysis.pdf'))
        LayoutMatrix<- matrix(c(1:2),nrow=2,ncol=1)
        PlotLayout<- layout(mat=LayoutMatrix)		
        plot(LogFrequency ~ Age,data=AgeDistMPA,xlab='Age',ylab='ln Frequency',bty='n',col=3,cex=2,pch=16,main='Inside Reserve',xlim=c(0,max(TempLengthDat$Age+1)+1),ylim=c(0,ceiling(max(AllHist$LogFrequency,na.rm=T))))
        lines(seq(from=BinBreaks[MPAPeak],to=BinBreaks[MPALastObserved],length.out=MPANumPoints),PredictedMPAValues,lty=2,lwd=2)
        text(ceiling(max(TempLengthDat$Age+1))*.65,.95*max(PredictedMPAValues),labels=paste('M=',round(NaturalMortality,2) ))
        plot(LogFrequency ~ Age,data=AgeDistFished,xlab='Age',ylab='ln Frequency',bty='n',col=2,cex=2,pch=16,main='Outside Reserve',xlim=c(0,max(TempLengthDat$Age+1)+1),ylim=c(0,ceiling(max(AllHist$LogFrequency,na.rm=T))))
        lines(seq(from=BinBreaks[FishedPeak],to=BinBreaks[FishedLastObserved],length.out=FishedNumPoints),PredictedFishedValues,lty=2,lwd=2)
        text(ceiling(max(TempLengthDat$Age+1))*.65,max(PredictedFishedValues)*.95,labels=paste('Z=',round(TotalMortality,2), '; F=', round(FishingMortality,2) ))
        dev.off()
        }
        else
        {
         pdf(file=paste(FigureFolder,i,Years[y],' Catch Curve Analysis.pdf'))
        plot(LogFrequency ~ Age,data=AgeDistFished,xlab='Age',ylab='ln Frequency',bty='n',col=2,cex=2,pch=16,main='Outside Reserve',xlim=c(0,max(TempLengthDat$Age+1)+1),ylim=c(0,ceiling(max(AllHist$LogFrequency,na.rm=T))))
        lines(seq(from=BinBreaks[FishedPeak],to=BinBreaks[FishedLastObserved],length.out=FishedNumPoints),PredictedFishedValues,lty=2,lwd=2)
        text(ceiling(max(TempLengthDat$Age+1))*.65,max(PredictedFishedValues)*.95,labels=paste('Z=',round(TotalMortality,2), '; F=', round(FishingMortality,2) ))
        dev.off()
    	
        } 
        
      }
      
      
      # MCOutput[c,]<- c(i,Years[y],'CatchCurve',FishingMortality/NaturalMortality,'FvFmsy',Flag)
      
      # MCDetails[c,]<- c(i,Years[y],FishingMortality,NaturalMortality)
      
      MCOutput$Iteration[c]<- i
      MCOutput$Year[c]<- Years[y]
      MCOutput$Method[c]<- 'CatchCurve'
      MCOutput$Value[c]<- FishingMortality/NaturalMortality
      MCOutput$Metric[c]<- 'FvM'
      MCOutput$Flag[c]<- Flag
	  MCOutput$SampleSize[c]<- SampleSize[y]

      MCDetails$Iteration[c]<- i
      MCDetails$Year[c]<- Years[y]
      MCDetails$FishingMortality[c]<- FishingMortality
      MCDetails$NaturalMortality[c]<- NaturalMortality

      
    }
    
    
    if  (length(MPAAllObserved)==1 | length(FishedAllObserved)==1 )
    {
     Flag<- 'Catch Curve cannot run - Only 1 Age group observed'
     MCOutput$Flag[c]<- Flag
           MCOutput$Iteration[c]<- i
      MCOutput$Year[c]<- Years[y]
      MCOutput$Method[c]<- 'CatchCurve'
      MCOutput$Value[c]<- NA
      MCOutput$Metric[c]<- 'FvM'
      MCOutput$Flag[c]<- Flag


      MCDetails$Iteration[c]<- i
      MCDetails$Year[c]<- Years[y]
      MCDetails$FishingMortality[c]<- NA      
      MCDetails$NaturalMortality[c]<- NA
  
    }
    MCOutput$Iteration[c]<- i

    }
    
  } #Close monte carlo loop
  
  
  #####################################
  ##### Proces Monte Carlo Data #######
  #####################################
  
  TrueIteration<- MCOutput$Iteration==1
  
  TrueOutput<- MCOutput[TrueIteration,]
  
  TrueDetails<- MCDetails[TrueIteration,]
  
  MCOutput<- MCOutput[TrueIteration==F,]
  
  Output$Year<- Years
  
  Output$Method<- 'CatchCurve'
  
  Output$Value<- TrueOutput$Value
  
  Output$LowerCI<- NA
  
  Output$UpperCI<- NA
  
  Output$SD<- NA
  
  Output$Metric<- 'FvM'
  
  Output$Flag<-TrueOutput$Flag 
  
  Output$SampleSize<- SampleSize
  
  if (Iterations>1)
  {
    
    for (y in 1:length(Years))
    {
      Where<- MCOutput$Year==Years[y]
      
      Temp<- MCOutput[Where,]
      
      if(sum(Temp$Value,na.rm=T)>0)
      {
      TempValue<- sort(as.numeric(Temp$Value))
      
      Bottom<- ceiling(.025*length(TempValue))
      
      Top<- ceiling(.975*length(TempValue))
      
      MeanMetric<- mean(as.numeric(Temp$Value),na.rm=T)
      
      LowerCI<- TempValue[Bottom]
      
      UpperCI<- TempValue[Top]
      
      SD<- sd(TempValue[Bottom:Top],na.rm=T)
      
      # Output[y,]<- c(Years[y],'CatchCurve',TrueOutput$Value[y],LowerCI,UpperCI,SD,'FvFmsy',TrueOutput$Flag[y])
      
      Output$Year[y]<- Years[y]
  
  Output$Method[y]<- 'CatchCurve'
  
  Output$Value[y]<- TrueOutput$Value[y]
  
  Output$LowerCI[y]<- LowerCI
  
  Output$UpperCI[y]<- UpperCI
  
  Output$SD[y]<- SD
  
  Output$Metric[y]<- 'FvM'
  
  Output$Flag[y]<-TrueOutput$Flag[y]
      }
      
    }
  }
  
  #####################################
  ##### Plot outputs #######
  #####################################
  
  pdf(file=paste(FigureFolder,' Catch Curve FvM Boxplots.pdf',sep=''))
  boxplot((MCDetails$FishingMortality/MCDetails$NaturalMortality)~MCDetails$Year,frame=F,xlab='Year',ylab='F/M',notch=T,outline=F,na.rm=T,width=SampleSize)
  dev.off()
  
  pdf(file=paste(FigureFolder,' Catch Curve Fishing Mortality Boxplots.pdf',sep=''))
  boxplot((MCDetails$FishingMortality)~MCDetails$Year,frame=F,xlab='Year',ylab='F',notch=T,outline=F, width=SampleSize)
  dev.off()
  
  pdf(file=paste(FigureFolder,' Catch Curve Natural Mortality Boxplots.pdf',sep=''))
  boxplot((MCDetails$NaturalMortality)~MCDetails$Year,frame=F,xlab='Year',ylab='M',notch=T,outline=F, width=SampleSize)
  dev.off()
  
  Fish<- BaseFish
  
  
    return(list(Output=Output,Details=MCDetails,MonteCarlo=MCOutput))		
  }



DensityRatio<- function(DenDat,LagLength,Weight,Form,Iterations,BootStrap)
{
  ##################
  ###### DensityRatio ######
  ##################
  #Source:Loosely based on methods described in McGilliard et al. 2011 and Babcock & MacCall 2011
  #Summary: Estimate N/K from density inside and outside of marine reserves
  
  # DenDat: The raw density data
  # LagLength: The number of years of lagged data to use
  # Weight: The weight assigned to historic data
  
  

  ############################
  ### Process Density Data ###	
  ############################
  
  ### Figure out timeline of data you want ###
  
  Years<- sort(unique(DenDat$Year))
  
  lag<- LagLength #lagged years
  
  weight<- Weight
  
  Output<- as.data.frame(matrix(NA,nrow=length(Years),ncol=8))
  
  colnames(Output)<- c('Year','Method','Value','LowerCI','UpperCI','SD','Metric','Flag')
  
  MCOutput<- as.data.frame(matrix(NA,nrow=length(Years)*Iterations,ncol=6))
  
  colnames(MCOutput)<- c('Iteration','Year','Method','Value','Metric','Flag')
  
  MCDetails<- as.data.frame(matrix(NA,nrow=length(Years)*Iterations,ncol=4))
  
  colnames(MCDetails)<- c('Iteration','Year','FishedDensity','MPADensity')
  
  Flag<- 'None'
  c<- 0
  
  BaseFish<- Fish
  
  for (i in 1:Iterations)
  {
    
    if (i>1)
    {
      Fish<- BaseFish
      
      Fish<- ApplyLifeHistoryError()
    }
    
    for (y in 1:length(Years))
    {
      
      c<- c+1
      
      Flag<- 'None'
      
      dYr<- y
      
      tempLaggedYears<- seq(from=dYr-lag,to=dYr,by=1)
      
      Lags<- t(as.matrix(tempLaggedYears[tempLaggedYears>0])) #Index of CPUEs that you want to grab for weighted calculation
      
      LaggedYears<- Years[Lags]
      
      weights<- weight^t(apply(Lags,1,rev)) #weight assigned to each year
      
      TempDenDat<- DenDat
      
      
      
       if (i>1 & BootStrap==1) #Resample length data inside and outside of MPAs
      {
        
       TempDenStorage<- as.data.frame(matrix(NA,nrow=0,ncol=dim(TempDenDat)[2]))
       
       colnames(TempDenStorage)<- colnames(TempDenDat)
       
        cc<- 0
         
        for (yy in 1:length(LaggedYears))
        {
        	
        
        MPADat<- TempDenDat[TempDenDat$MPA==1 & TempDenDat$Year==LaggedYears[yy],]
        
        FishedDat<- TempDenDat[TempDenDat$MPA==0 & TempDenDat$Year==LaggedYears[yy],]
        
        NumPoints<- 1:dim(MPADat)[1]	
        
        BootSample<- sample(NumPoints,length(NumPoints),replace=T)
        
        TempMPADat<- MPADat[BootSample,]
        
        NumPoints<- 1:dim(FishedDat)[1]	
        
        BootSample<- sample(NumPoints,length(NumPoints),replace=T)
        
        TempFishedDat<- FishedDat[BootSample,]
        
        TempDat<- rbind(TempFishedDat,TempMPADat)
        
        Size<- dim(TempDat)[1]
        
        TempDenStorage[cc+1:Size,]<- TempDat
        
        cc<- cc+Size
        
        }
        
        TempDenDat<- TempDenStorage
      }

      
      WeightedDensity<- CalculateDensity(TempDenDat,LaggedYears,weights,Form)
      
      inside<- WeightedDensity$MPADensity
      
      outside<- WeightedDensity$FishedDensity
      
      WeightedRatio<- outside/inside
      
      WeightedIn<- inside
      
      WeightedOut<- outside
      
      
      MCOutput$Iteration[c]<- i
      
      MCOutput$Year[c]<- Years[y]
      
      MCOutput$Method[c]<- 'DensityRatio'
      
      MCOutput$Value[c]<- WeightedRatio
      
      MCOutput$Metric[c]<- 'N/k'
      
      MCOutput$Flag[c]<- Form
      
      MCDetails$Iteration[c]<- i
      
      MCDetails$Year[c]<- Years[y]
      
      MCDetails$FishedDensity[c]<- WeightedOut
      
      MCDetails$MPADensity[c]<- WeightedIn
      
      
    } #Close year loop		
    
    
  } #Close iteration loop
  
  
  TrueIteration<- MCOutput$Iteration==1
  
  TrueOutput<- MCOutput[TrueIteration,]
  
  TrueDetails<- MCDetails[TrueIteration,]
  
  MCOutput<- MCOutput[TrueIteration==F,]
  
  Output$Year<- Years
  
  Output$Method<- 'DensityRatio'
  
  Output$Value<- TrueOutput$Value
  
  Output$LowerCI<- NA
  
  Output$UpperCI<- NA
  
  Output$SD<- NA
  
  Output$Metric<- 'N/K'
  
  Output$Flag<-TrueOutput$Flag 
  
  if (Iterations>1)
  {
    
    for (y in 1:length(Years))
    {
      Where<- MCOutput$Year==Years[y]
      
      Temp<- MCOutput[Where,]
      
      TempValue<- sort(as.numeric(Temp$Value))
      
      # pdf(file=paste(FigureFolder,Years[y],' Density Ratio Histogram.pdf',sep=''))
      # hist(TempValue,xlab='N/K',main=NA)
      # dev.off()
      
      Bottom<- ceiling(.025*length(TempValue))
      
      Top<- ceiling(.975*length(TempValue))
      
      MeanMetric<- mean(as.numeric(Temp$Value),na.rm=T)
      
      LowerCI<- TempValue[Bottom]
      
      UpperCI<- TempValue[Top]
      
      SD<- sd(TempValue[Bottom:Top],na.rm=T)
      
      Output[y,]<- c(Years[y],'DensityRatio',TrueOutput$Value[y],LowerCI,UpperCI,SD,'N/K',TrueOutput$Flag[y])
      
    }
    
  }
  
  MCDetails$MPADensity[MCDetails$MPADensity==0]<- NA
  
  pdf(file=paste(FigureFolder,' Density Ratio NvK Boxplots.pdf',sep=''))
  boxplot((MCDetails$FishedDensity/MCDetails$MPADensity)~MCDetails$Year,frame=F,xlab='Year',ylab='Fished Density/Unfished Density',notch=T,outline=F,ylim=c(0,5))
  dev.off()
  
  pdf(file=paste(FigureFolder,' Density Ratio Fished Density Boxplots.pdf',sep=''))
  boxplot((MCDetails$FishedDensity)~MCDetails$Year,frame=F,xlab='Year',ylab='Density',notch=T,outline=F)
  dev.off()
  
  pdf(file=paste(FigureFolder,' Density Ratio MPA Density Boxplots.pdf',sep=''))
  boxplot((MCDetails$MPADensity)~MCDetails$Year,frame=F,xlab='Year',ylab='Density',notch=T,outline=F)
  dev.off()
  
  Fish<- BaseFish
  
  return(list(Output=Output,Details=MCDetails))
  
}

LBSPR<-function(LengthDat,EstimateM,Iterations,BootStrap,LifeError,LengthBins)
{
  ######################
  ###### LBSPR #########
  ######################
  #Source: Based on the length based SPR methods developed by Prince, Valencia, Adrian
  #Summary: Estimate SPR by estimating selectivity, F using observed and predicted length frequency data
  
   # setwd("/Users/danovando/Desktop/Bren/SFG Work/DPSA")
  
# #   LengthDat<- LengthData
  
  # EstimateM<- 0   
  
  # Iterations<-  10 
  
  # BootStrap<- 1
  
  # LifeError<- 1
  
  # LengthBins<- 1
  
  ######################
  ###### Inputs #########
  ######################
  
  #   LengthDat: LengthData
  #   EstimateM: 1 estiamtes M using catch curve 
  #   Iterations: Number of iterations to run (1 runs only on the default data)
  #   BootStrap: 1 bootstaps the length data in the monte carlo run
  
  
  ######################
  ###### Process Data #########
  ######################
  LengthDat<- LengthDat[is.na(LengthDat$Length)==F,]
  
  # LengthDat$Length[LengthDat$Length>Fish$Linf]<- Fish$Linf * 0.98
  
  # LengthDat<- LengthDat[LengthDat$Length<Fish$Linf,]
  
  SampleSize<- NA
  
  Years<- sort(unique(LengthDat$Year))
  
  Output<- as.data.frame(matrix(NA,nrow=length(Years),ncol=9))
  
  colnames(Output)<- c('Year','Method','SampleSize','Value','LowerCI','UpperCI','SD','Metric','Flag')
  
  MCOutput<- as.data.frame(matrix(NA,nrow=length(Years)*Iterations,ncol=9))
  
  colnames(MCOutput)<- c('Iteration','Year','Method','Value','LowerCI','UpperCI','SD','Metric','Flag')
  
  Flag<- 'None'
  
  Details<- as.data.frame(matrix(NA,nrow=length(Years),ncol=4))
  
  colnames(Details)<- c('Year','FvM','SelL50','SelL95')
  
  MCDetails<- as.data.frame(matrix(NA,nrow=length(Years)*Iterations,ncol=5))
  
  colnames(MCDetails)<- c('Iteration','Year','FvM','SelL50','SelL95')
  
  ###################################
  # Source LBSPR function #    
  #################################### 
  CurrentDir<- getwd()
  AssessDir <- paste(CurrentDir,"LBSPR",sep="/")
  setwd(AssessDir)
  # source("LBSPRAssessmentFun_May21.R")
  ###################################
  # Source all other functions needed#    
  #################################### 
  
  Call.functions.Fun    <- function(dropboxDir) {
    cd <- dropboxDir
    File.Path  <- paste(cd,"All_Functions", sep="/")
    for (nm in list.files(File.Path)) {
      source(file.path(File.Path, nm))
    }
  }
  Call.functions.Fun(AssessDir)	
  
  ###################################
  # load R2ADMB and compile tpl file #    
  ###################################
  #install.packages("R2admb")  #You will need to install the R2admb package the first time, and make sure 
  #the path is set so that R and admb can communicate.
  library("R2admb")
  compile_admb("LBSPR_AssessFun",verbose=FALSE)
  setwd(CurrentDir)  #changes directory back to main folder.
  
  ##Set Length bins
  # LengthBins <- 5 #1cm
  
  BaseFish<- Fish
  
  c<- 0
  
  CatchCurveResiduals<- as.data.frame(matrix(NA,nrow=0,ncol=2))

  AgeResiduals<- as.data.frame(matrix(NA,nrow=0,ncol=2))

  
  for (i in 1:Iterations)
  {
    
    if (i>1 & LifeError==1) #Apply life history error
    {
      Fish<- BaseFish
      Fish<- ApplyLifeHistoryError()
    }
    
    for (y in 1:length(Years)) #loop over years
    {
      #Read in the Size data and Species parameter file
      c<- c+1
      
      WhereFished<- LengthDat$Year==Years[y] & LengthDat$MPA==0
      
      WhereAll<- LengthDat$Year==Years[y] 
      
      AnyMPA<- sum(LengthDat$MPA[WhereAll])>0
      
      ## Convert lenth data into appropriate format-- vector of raw size data
      CatchatLength <- LengthDat$Length[WhereFished]
      
      SampleSize[y]<- length(CatchatLength)
      # CatchatLength <- LengthDat$Length[WhereAll]
      
      if (i>1 & BootStrap==1) #Bootstrap length data 
      {
        
        NumPoints<- 1:length(CatchatLength)
        
        BootSample<- sample(NumPoints,length(NumPoints),replace=T)
        
        CatchatLength<- CatchatLength[BootSample]	
      }
      
      
      EstimatedM<-'None'
      
      if (EstimateM==1 & AnyMPA==T) #Estimate M using catch curve if you want
      {
        AllCatchAtLength<- LengthDat[WhereAll,]
        
        EstimatedM<- CatchCurve(AllCatchAtLength,1,1,1,0,0,1)$Details$NaturalMortality
        
      }
      
      # LengthDat <- read.csv("ExampleData.csv")
      # CatchatLength <- LengthDat$Redtail.F
      
      # Call the Assessment Function
      #print(paste("Pre-Funct LengthBins",LengthBins))
      
     
      Estimates <- LBSPR_SingleSpeciesAssessmentfun(CatchatLength,AssessDir,CurrentDir,LengthBins,Years[y],EstimatedM,Fish) #Run LBSPR in ADMB (see SubFunctions)
      
      
      CatchCurveResiduals<- rbind(CatchCurveResiduals,Estimates$CatchCurveResiduals)


      AgeResiduals<- rbind(AgeResiduals,Estimates$AgeResiduals)
      
      MCOutput[c,]<- c(i,Estimates$Output)
      
      MCDetails[c,]<- c(i,Estimates$Details)
      # Sel50, Sel95, F/M, and SPR stored in Estimates.
    } #Close year loop	
  } #Close iteration loop
  
  
  ########################################
  ###### Process Monte Carlo Data #########
  #########################################
  
  TrueIteration<- MCOutput$Iteration==1
  
  TrueOutput<- MCOutput[TrueIteration,]
  
  TrueDetails<- MCDetails[TrueIteration,]
  
  # MCOutput<- MCOutput[TrueIteration==F,]
  
  Output$Year<- Years
  
  Output$Method<- 'LBSPR'
  
  Output$Value<- TrueOutput$Value
  
  Output$LowerCI<- NA
  
  Output$UpperCI<- NA
  
  Output$SD<- NA
  
  Output$Metric<- 'SPR'
  
  Output$Flag<-TrueOutput$Flag 
  
  Output$SampleSize<- SampleSize
  
  if (Iterations>1)
  {
    
    for (y in 1:length(Years))
    {
      Where<- MCOutput$Year==Years[y]
      
      Temp<- MCOutput[Where,]
      
      TempValue<- sort(as.numeric(Temp$Value))
      
      Bottom<- ceiling(.025*length(TempValue))
      
      Top<- ceiling(.975*length(TempValue))
      
      MeanMetric<- mean(as.numeric(Temp$Value),na.rm=T)
      
      LowerCI<- TempValue[Bottom]
      
      UpperCI<- TempValue[Top]
      
      SD<- sd(TempValue[Bottom:Top],na.rm=T)
      
      Output[y,]<-c(Years[y],'LBSPR',SampleSize[y],TrueOutput$Value[y],LowerCI,UpperCI,SD,'SPR',TrueOutput$Flag[y])
      
    }
    
  }
  
  ######################
  ###### Make Plots #########
  ######################
  
  pdf(file=paste(FigureFolder,'Age Residuals Boxplots.pdf',sep=''))
  boxplot(AgeResiduals$Residuals ~ AgeResiduals$Age,frame=F,xlab='Age',ylab='Residuals',notch=T,outline=F)
  abline(h=0)
  dev.off()
  

  pdf(file=paste(FigureFolder,'Cohort Residuals Boxplots.pdf',sep=''))
  boxplot(CatchCurveResiduals$Residuals ~ CatchCurveResiduals$Cohort,frame=F,xlab='Cohort',ylab='Residuals',notch=T,outline=F)
  abline(h=0)
  dev.off()
  
  
  pdf(file=paste(FigureFolder,' LBSPR SPR Boxplots.pdf',sep=''))
  boxplot(MCOutput$Value ~ MCOutput$Year,frame=F,xlab='Year',ylab='SPR',notch=T,outline=F,width=SampleSize)
  dev.off()
  
  pdf(file=paste(FigureFolder,' LBSPR FvM Boxplots.pdf',sep=''))
  boxplot((MCDetails$FvM)~MCDetails$Year,frame=F,xlab='Year',ylab='F/M',notch=T,outline=F, width=SampleSize)
  dev.off()
  
  
  Fish<- BaseFish
  
  return(list(Output=Output,Details=MCDetails))
}

UnishedDensity<- function(DenDat,LagLength,Weight,Form,Iterations,BootStrap)
{
  ##################
  ###### UnfishedDensity ######
  ##################
  #Source:Loosely based on methods McClanahan et al. 2011
  #Summary: Estimate 
  
  # DenDat: The raw density data
  # LagLength: The number of years of lagged data to use
  # Weight: The weight assigned to historic data
  
  

  ############################
  ### Process Density Data ###	
  ############################
  
  ### Figure out timeline of data you want ###
  
  Years<- sort(unique(DenDat$Year))
  
  lag<- LagLength #lagged years
  
  weight<- Weight
  
  Output<- as.data.frame(matrix(NA,nrow=length(Years),ncol=8))
  
  colnames(Output)<- c('Year','Method','Value','LowerCI','UpperCI','SD','Metric','Flag')
  
  MCOutput<- as.data.frame(matrix(NA,nrow=length(Years)*Iterations,ncol=6))
  
  colnames(MCOutput)<- c('Iteration','Year','Method','Value','Metric','Flag')
  
  MCDetails<- as.data.frame(matrix(NA,nrow=length(Years)*Iterations,ncol=4))
  
  colnames(MCDetails)<- c('Iteration','Year','FishedDensity','MPADensity')
  
  Flag<- 'None'
  c<- 0
  
  BaseFish<- Fish
  
  for (i in 1:Iterations)
  {
    
    if (i>1)
    {
      Fish<- BaseFish
    }
    
    for (y in 1:length(Years))
    {
      
      c<- c+1
      
      Flag<- 'None'
      
      dYr<- y
      
      tempLaggedYears<- seq(from=dYr-lag,to=dYr,by=1)
      
      Lags<- t(as.matrix(tempLaggedYears[tempLaggedYears>0])) #Index of CPUEs that you want to grab for weighted calculation
      
      LaggedYears<- Years[Lags]
      
      weights<- weight^t(apply(Lags,1,rev)) #weight assigned to each year
      
      TempDenDat<- DenDat
      
      
      
       if (i>1 & BootStrap==1) #Resample length data inside and outside of MPAs
      {
        
       TempDenStorage<- as.data.frame(matrix(NA,nrow=0,ncol=dim(TempDenDat)[2]))
       
       colnames(TempDenStorage)<- colnames(TempDenDat)
       
        cc<- 0
         
        for (yy in 1:length(LaggedYears))
        {
        	
        
        MPADat<- TempDenDat[TempDenDat$MPA==1 & TempDenDat$Year==LaggedYears[yy],]
        
        FishedDat<- TempDenDat[TempDenDat$MPA==0 & TempDenDat$Year==LaggedYears[yy],]
        
        NumPoints<- 1:dim(MPADat)[1]	
        
        BootSample<- sample(NumPoints,length(NumPoints),replace=T)
        
        TempMPADat<- MPADat[BootSample,]
        
        NumPoints<- 1:dim(FishedDat)[1]	
        
        BootSample<- sample(NumPoints,length(NumPoints),replace=T)
        
        TempFishedDat<- FishedDat[BootSample,]
        
        TempDat<- rbind(TempFishedDat,TempMPADat)
        
        Size<- dim(TempDat)[1]
        
        TempDenStorage[cc+1:Size,]<- TempDat
        
        cc<- cc+Size
        
        }
        
        TempDenDat<- TempDenStorage
      }

      
      WeightedDensity<- CalculateDensity(TempDenDat,LaggedYears,weights,Form)
      
      inside<- WeightedDensity$MPADensity
      
      outside<- WeightedDensity$FishedDensity
      
      WeightedRatio<- outside/inside
      
      WeightedIn<- inside*10
      
      WeightedOut<- outside*10
      
      
      MCOutput$Iteration[c]<- i
      
      MCOutput$Year[c]<- Years[y]
      
      MCOutput$Method[c]<- 'UnfishedDensity'
      
      MCOutput$Value[c]<- WeightedIn
      
      MCOutput$Metric[c]<- 'kg/ha'
      
      MCOutput$Flag[c]<- Form
      
      MCDetails$Iteration[c]<- i
      
      MCDetails$Year[c]<- Years[y]
      
      MCDetails$FishedDensity[c]<- WeightedOut
      
      MCDetails$MPADensity[c]<- WeightedIn
      
      
    } #Close year loop		
    
    
  } #Close iteration loop
  
  
  TrueIteration<- MCOutput$Iteration==1
  
  TrueOutput<- MCOutput[TrueIteration,]
  
  TrueDetails<- MCDetails[TrueIteration,]
  
  MCOutput<- MCOutput[TrueIteration==F,]
  
  Output$Year<- Years
  
  Output$Method<- 'UnfishedDensity'
  
  Output$Value<- TrueOutput$Value
  
  Output$LowerCI<- NA
  
  Output$UpperCI<- NA
  
  Output$SD<- NA
  
  Output$Metric<- 'kg/ha'
  
  Output$Flag<-TrueOutput$Flag 
  
  if (Iterations>1)
  {
    
    for (y in 1:length(Years))
    {
      Where<- MCOutput$Year==Years[y]
      
      Temp<- MCOutput[Where,]
      
      TempValue<- sort(as.numeric(Temp$Value))
      
      
      Bottom<- ceiling(.025*length(TempValue))
      
      Top<- ceiling(.975*length(TempValue))
      
      MeanMetric<- mean(as.numeric(Temp$Value),na.rm=T)
      
      LowerCI<- TempValue[Bottom]
      
      UpperCI<- TempValue[Top]
      
      SD<- sd(TempValue[Bottom:Top],na.rm=T)
      
      Output[y,]<- c(Years[y],'UnfishedDensity',TrueOutput$Value[y],LowerCI,UpperCI,SD,'kg/ha',TrueOutput$Flag[y])
      
    }
    
  }
  
  MCDetails$MPADensity[MCDetails$MPADensity==0]<- NA
  
  pdf(file=paste(FigureFolder,'Unfished Density Boxplots.pdf',sep=''))
  boxplot((MCDetails$MPADensity)~MCDetails$Year,frame=F,xlab='Year',ylab='Unfished Density [kg/ha]',notch=T,outline=F)
  dev.off()
   
  Fish<- BaseFish
  
  return(list(Output=Output,Details=MCDetails))
  
}


DTREE<- function(CPUEDat,LengthDat,Options)
{
  
}