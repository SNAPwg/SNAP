RepMat<- function(Data,Reps)
{
  # RepMat: A function to create repeated entries of a dataframe or matrix ---------------
  # Data is the row(s) or column(s) that you want to replice. 
  # Reps is the number of times you want it replicated
  # Direction indicates whether you want to rep rows or columns
  
  #   Data<- ((TempStock[TempStock$Year==MaxYear,]))
  
  RepDex<- rep(1:dim(Data)[1],Reps)
  
  RepData<- Data[RepDex,]
  
  return(RepData)
}


