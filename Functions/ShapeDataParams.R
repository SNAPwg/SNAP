ShapeDatParams<-function(datafile)
{
DataParams <- NULL
DataParams$histCatchFD    <-Samp[grep('histCatchFD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]              # Collect historical Catch Data?
DataParams$histEffortFD    <-Samp[grep('histEffortFD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]           	  # Collect historical CPUE Data?
DataParams$histStartYr	 	<-Samp[grep('histStartYr',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	            # time step in which historical data was first collected
DataParams$Aggregate   	 	<-Samp[grep('Aggregate',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	# Aggregate data over space? 1 for yes, 2 for no.
DataParams$SampStartYr	 	<-Samp[grep('SampStartYr',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]            	# time step in which current sampling program begins
DataParams$SampFreq		    <-Samp[grep('SampFreq',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]       	# Sampling Frequency
DataParams$sigHistCatch  	<-Samp[grep('sigHistCatch',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]         		  	# Standard deviation of observation error on total catches. If no error, sigHistCatch=0.
DataParams$catchBias	  	<-Samp[grep('catchBias',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	              # Bias in catch reporting. Scalar between (-1, 1). 0 equals no bias.
DataParams$Sel50FI       	<-Samp[grep('Sel50FI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	      # Survey selectivity at age
DataParams$Sel95FI	     	<-Samp[grep('Sel95FI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]            	# Survey selectivity at age
DataParams$Sel50VB		    <-Samp[grep('Sel50VB',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]       	#  Selectivity for von Bertlanffy surveys (usually uniform above a specific age)
DataParams$Sel95VB    	 	<-Samp[grep('Sel95VB',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]         	# Selectivity for von Bertlanffy surveys (usually uniform above a specific age)
DataParams$SurveyF	    	<-Samp[grep('SurveyF',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	          # Fishing rate associated with survey effort
DataParams$nAgeFD         <-Samp[grep('nAgeFD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]   # Number of fishery dependent age samples collected from each patch
DataParams$nAgeFI         <-Samp[grep('nAgeFI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]              # Number of fishery independent age samples collected from each patch
DataParams$nSizeFD      	<-Samp[grep('nSizeFD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]           	  # Number of fishery dependent size samples collected from each patch
DataParams$nSizeFI      	<-Samp[grep('nSizeFI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	            # Number of fishery independent size samples collected from each patch
DataParams$numVB      	 	<-Samp[grep('numVB',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	# Number of VB age-size samples collected from each patch
DataParams$VBbins	      	<-Samp[grep('VBbins',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]            	#  Size bins associated with von Bertalanffy survey
DataParams$lengthObsSD_FI <-Samp[grep('lengthObsSD_FI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]              # Observation Error on fishery independent lengths
DataParams$lengthObsSD_FD <-Samp[grep('lengthObsSD_FD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]           	  # Observation Error on fishery dependent lengths
DataParams$AgeObsSD_FI	 	<-Samp[grep('AgeObsSD_FI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	            # Observation Error on fishery independent ages
DataParams$AgeObsSD_FD   	<-Samp[grep('AgeObsSD_FD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	# Observation Error on fishery dependent ages
DataParams$FISurvType   	<-Samp[grep('FISurvType',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]            	# Spatial pattern for fishery independent surveys (0 = input a .csv file specifying custom survey distribution, 1 = sample all patches, 2=sample no patches, 3 = sample random 20% of patches, 4= sample near and far site, 5=sample MPAs)
DataParams$FDSurvType		  <-Samp[grep('FDSurvType',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]       	# Spatial pattern for fishery dependent surveys (0 = input a .csv file specifying custom survey distribution, 1 = sample all fished patches, 2=sample no patches, 3 = sample random 20% of patches, 4= sample near and far site)
DataParams$Catch_FD       <-Samp[grep('Catch_FD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	            # Collect fishery dependent catch data? 1 = yes, 2 = no.
DataParams$Catch_FI   	 	<-Samp[grep('Catch_FI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	# Collect fishery independent catch data? 1 = yes, 2 = no.
DataParams$Effort_FI	  	<-Samp[grep('Effort_FI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]            	# Collect fishery independent abundance index data? 1 = yes, 2 = no.
DataParams$Effort_FD      <-Samp[grep('Effort_FD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]              # Collect fishery dependent CPUE data? 1 = yes, 2 = no.
DataParams$Ages_FD   	    <-Samp[grep('Ages_FD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]           	  # Collect fishery dependent age composition data? 1 = yes, 2 = no.
DataParams$Ages_FI	    	<-Samp[grep('Ages_FI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	            # Collect fishery independent age composition data? 1 = yes, 2 = no.
DataParams$Sizes_FD   	 	<-Samp[grep('Sizes_FD',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]	# Collect fishery dependent size composition data? 1 = yes, 2 = no.
DataParams$Sizes_FI     	<-Samp[grep('Sizes_FI',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]            	# Collect fishery independent size composition data? 1 = yes, 2 = no.
DataParams$VBSurvey		    <-Samp[grep('VBSurvey',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]       	# Collect age and size data to estimate VB params? 1 = yes, 2 = no.
DataParams$Survey_q       <-Samp[grep('Survey_q',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]            	# Survey catchability. Might be different depending on mode of survey (such as dive surveys)
DataParams$sigSurvey		    <-Samp[grep('sigSurvey',Samp[,ncol(Samp)]),seq(1,ncol(Samp)-1)]       	# Standard deviation of observation error on total catches. If no error, sigHistCatch=0.
return(DataParams)
}