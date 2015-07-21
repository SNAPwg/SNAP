#' This function Predict the value of kslope parameter from the LM and Pars.  
#' @name PredictKSlope  
#' @title Predict the value of kslope parameter from the LM and Pars.
#' @param Pars An object of class \code{list} that contains all parameters
#'   required to run GTG model. 
#' @author Adrian Hordyk 
#' @export

PredictKSlope <- function(Pars) {
  data(KSlopeMod)
  KSlope <- exp(predict(KSlopeMod, newdata=list(LinfMK=log(Pars$Linf)-log(Pars$MK), CVLinf=Pars$CVLinf, RelL50=(Pars$L50/Pars$Linf), MatDelta=Pars$L95 - Pars$L50,  Mpow=Pars$Mpow)))  
  return(as.numeric(KSlope))
}
