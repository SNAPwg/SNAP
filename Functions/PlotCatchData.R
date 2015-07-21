PlotCatchData<- function(CatchDat,FigureFolder)
{
  
  library(zoo,quietly = T)
  # CatchDat=CatchData
  
  InterpCatch<- na.approx(CatchDat$Catch)
  
  InterpMarker=as.numeric(is.na(CatchDat$Catch))
  
  PointStyle=21*InterpMarker
  PointStyle[PointStyle==0]<- 16
  
  pdf(file=paste(FigureFolder,' Catch History.pdf',sep=''),width=7,height=4) 
  par(mai=c(1,1,1,1.5))
  Month<- CatchDat$Month
  Month[Month==-999]<- 'Total'
  TimeName=paste(CatchDat$Year,Month,sep='-')
  plot(InterpCatch,xaxt='n',pch=PointStyle,col=InterpMarker+1,ylab=CatchDat$Units[1],xlab=NA,cex=1.2,pty='m',bty='n')
  axis(1,at=1:length(TimeName),labels=TimeName,las=0)
  legend('right',pch=c(16,21),col=c(1,2),legend=c('Real','Interpolated'),xpd=T,bty='n',inset=-.3)
  dev.off()
  
}