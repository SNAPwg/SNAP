VisualizeMovement<-function()
{
#==visualize movement
tempSpace<-array(dim=c(SpaceR,SpaceC,SpaceC*SpaceR))
 coords<-NULL
 for(x in 1:SpaceR)
  {
   temp<-cbind(rep(x,SpaceC),seq(1,SpaceC))
   coords<-rbind(coords,temp)
  }

for(x in 1:10)
{
 dev.new()
 filled.contour(SpaceNum)
for(h in 1:nrow(coords))
 {
 tempSpace[,,h]<-Movement[,,h]*SpaceNum[coords[h,1],coords[h,2]]
 }
 SpaceNum<-apply(tempSpace,c(1,2),sum)
}
}