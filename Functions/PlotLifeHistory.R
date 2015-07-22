PlotLifeHistory<- function(Fish)
{

  Lengths<- LengthAtAge(1:Fish$MaxAge,Fish,0)

  Maturity<- MaturityAtAge(Lengths,Fish)

  LifeMat<- data.frame(1:Fish$MaxAge,Lengths,Maturity)

  colnames(LifeMat)<- c('Age','Length','PercentMature')


  LifePlot<- (ggplot(LifeMat,aes(x=Age,y=Length,color=PercentMature))+geom_line(size=2) +
              ylab('Length (mm)') +
              scale_colour_gradient(limits=c(0, 1), low="steelblue2",high='red'))

  return(LifePlot)
  dev.off()
}
