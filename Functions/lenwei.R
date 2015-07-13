lenwei<-function(Linf,K,t0,a,b,kmax)
{
 #Length weight at age
 l<-rep(0,kmax)
 w<-rep(0,kmax)
 for(k in 1:kmax)
  {
   l[k]<-Linf*(1-exp(-K*k-t0))
   w[k]<-a*l[k]^b
  }
 list(l=l,w=w)
}