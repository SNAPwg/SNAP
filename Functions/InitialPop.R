InitialPop<-function(R0,M,kmax,SpaceR,SpaceC,MoveProb,coords,Movement,Graphs=T,burn)
{
SpaceNumAtAgeT<-array(dim=c(burn,SpaceR,SpaceC,kmax))
#==set up initial population============
SpaceNumAtAgeT[1,,,1]<-R0/(SpaceR*SpaceC)
for(h in 2:kmax)
 SpaceNumAtAgeT[1,,,h]<-SpaceNumAtAgeT[1,,,h-1]*exp(-M)
#SpaceNumAtAgeT[1,,,]
#==plus group
SpaceNumAtAgeT[1,,,kmax]<-SpaceNumAtAgeT[1,,,kmax-1]*exp(-M)/(1-exp(-M))

#==allow movement to equilibrate==========
if(Graphs==T)
{
filled.contour(SpaceNumAtAgeT[1,,,1])
}
for(time in 2:burn)
{
 for(age in 1:kmax)
  {
   tempAge<-array(dim=c(SpaceR,SpaceC,SpaceR*SpaceC))
   Movers<-MoveProb[age]*SpaceNumAtAgeT[time -1,,,age]
   NonMovers<-(1-MoveProb[age])*SpaceNumAtAgeT[time -1,,,age]

   for(h in 1:nrow(coords))
    tempAge[,,h]<-Movers[coords[h,1],coords[h,2]]*Movement[,,h]

   Moved<-apply(tempAge,c(1,2),sum)
   SpaceNumAtAgeT[time ,,,age]<-Moved+NonMovers
   }
if(Graphs==T)
{
filled.contour(SpaceNumAtAgeT[time ,,,1])
}

}
return(SpaceNumAtAgeT[1:burn,,,])
}
