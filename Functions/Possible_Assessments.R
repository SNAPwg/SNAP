Find_Possible_Assessments <- function(Fishery, Data, time_steps)
{

  Files <- names(Data)

  assessments <- c('cpue_trend', 'catch_trend','length_look')

  possible_requirements <- c('assessment','n_lengths','y_catches','y_cpues')

  assessment_requirements <- as.data.frame(matrix(NA,nrow = length(assessments), ncol = length(possible_requirements)))

  colnames(assessment_requirements) <- possible_requirements

  assessment_requirements$assessment <- assessments

  assessment_requirements$n_lengths[assessment_requirements$assessment == 'length_look'] <- 100

  assessment_requirements$y_catches[assessment_requirements$assessment == 'catch_trend'] <- 2

  assessment_requirements$y_cpues[assessment_requirements$assessment == 'cpue_trend'] <- 2

  if ( any(grepl('_LengthData',Files)) )
  {
    LengthData <- Data[[which(grepl('_LengthData',Files))]]

  }
  if ( any(grepl('_CatchData',Files)) )
  {
    CatchData <-  Data[[which(grepl('_CatchData',Files))]]
  }
  if ( any(grepl('_EffortData',Files)) )
  {
    EffortData <- Data[[which(grepl('_EffortData',Files))]]
  }
  if ( any(grepl('_CPUEData',Files)) )
  {
    CPUEData <- Data[[which(grepl('_CPUEData',Files))]]
  }
  if ( any(grepl('_DensityData',Files)) )
  {
    DensityData <- Data[[which(grepl('_DensityData',Files))]]
  }
  if ( any(grepl('_MapData',Files)) )
  {
    Locations <- Data[[which(grepl('_MapData',Files))]]
  }

  length_samples <-  LengthData %>%
    filter(TimeStep %in% time_steps) %>%
    summarize(n_lengths = sum(is.na(Length)==F))

  catches <-  CatchData %>%
    summarize(y_catches = sum(is.na(Catch)==F))

  cpues <-  CPUEData %>%
    summarize(y_cpues = sum(is.na(Index)==F))

  what_we_have <- data.frame(n_lengths = length_samples,y_catches = catches, y_cpues = cpues)


  whats_missing <- data.frame(assessment = assessments, what_we_have[rep(1,length(assessments)),]
                              >= assessment_requirements[,2:dim(assessment_requirements)[2]], stringsAsFactors = F)

  total_requirements <- assessment_requirements %>%
    gather('requirement','value', 2:dim(assessment_requirements)[2]) %>%
    group_by(assessment) %>%
    summarise(number_of_requirements = sum(is.na(value)==F))

  requirements_available <- whats_missing %>%
    gather('requirement','sufficient', 2:dim(whats_missing)[2]) %>%
    group_by(assessment) %>%
    summarise(number_of_requirements_met = sum(sufficient==T & is.na(sufficient) ==F))

  what_can_we_do <- join(total_requirements,requirements_available, by = 'assessment')

  possible_assessments <- with(what_can_we_do, assessment[number_of_requirements == number_of_requirements_met])

  return(list(possible_assessments = possible_assessments, whats_missing= whats_missing ))
}