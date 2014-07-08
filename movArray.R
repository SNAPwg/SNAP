
movArray<-function(SpaceR,SpaceC,sdx,sdy)
{
 SpaceIn<-array(dim=c(SpaceR,SpaceC,SpaceC*SpaceR))
 # bivariate normal diffusion
 # list the coordinates 
 coords<-NULL
 for(x in 1:SpaceR)
  {
   temp<-cbind(rep(x,SpaceC),seq(1,SpaceC))
   coords<-rbind(coords,temp)
  }

 # population the probability array of moving 
 for(h in 1:nrow(coords))
  for(j in 1:SpaceR)
   for(k in 1:SpaceC)
     SpaceIn[j,k,h]<-exp(-(((coords[h,1]-j))^2/(2*sdx^2)+((coords[h,2]-k))^2/(2*sdy^2)))

 # normalize so that it sums to 1 (effectively wraps around)
 for(h in 1:nrow(coords))
  SpaceIn[,,h]<-SpaceIn[,,h]/sum(SpaceIn[,,h])

return(SpaceIn)
}

#==directional (add a rho?)   
   


