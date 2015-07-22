LengthAtAge<- function(Ages,Fish,Error)
{
  LenSD<- Error*(1+Fish$VBErrorSlope*Ages/Fish$MaxAge)

  RawLengths<- Fish$Linf*(1-exp(-Fish$vbk*(Ages-Fish$t0)))

  LengthWithError<- RawLengths*rlnorm(length(Ages),meanlog = 0, sdlog = LenSD)

  return(LengthWithError)
}