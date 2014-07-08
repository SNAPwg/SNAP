InitialPop<-function(R0,M,kmax,SpaceR,SpaceC,MoveProb,coords,Movement,Graphs=T,burn)
{
#==set up initial population============
SpaceNumAtLenT[1,,,1]<-R0/(SpaceR*SpaceC)
for(h in 2:kmax)
 SpaceNumAtLenT[1,,,h]<-SpaceNumAtLenT[1,,,h-1]*exp(-M)

#==plus group
SpaceNumAtLenT[1,,,kmax]<-SpaceNumAtLenT[1,,,kmax-1]*exp(-M)/(1-exp(-M))

#==allow movement to equilibrate==========
if(Graphs==T)
{
dev.new()
filled.contour(SpaceNumAtLenT[1,,,1])
}
for(time in 2:burn)
{
 for(age in 1:kmax)
  {
   tempAge<-array(dim=c(SpaceR,SpaceC,SpaceR*SpaceC))
   Movers<-MoveProb[age]*SpaceNumAtLenT[time -1,,,age]
   NonMovers<-(1-MoveProb[age])*SpaceNumAtLenT[time -1,,,age]

   for(h in 1:nrow(coords))
    tempAge[,,h]<-Movers[coords[h,1],coords[h,2]]*Movement[,,h]

   Moved<-apply(tempAge,c(1,2),sum)
   SpaceNumAtLenT[time ,,,age]<-Moved+NonMovers
   }
if(Graphs==T)
{
dev.new()
filled.contour(SpaceNumAtLenT[time ,,,1])
}
}
return(SpaceNumAtLenT[1:burn,,,])
}
