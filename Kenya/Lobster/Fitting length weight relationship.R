setwd("~/Dropbox/SNAP/SNAP Working Group Repository/SNAPwg/Kenya/Lobster")
Lobs <- read.csv("Length-weight.csv")

str(Lobs)
attach(Lobs)
plot(Length,Weight)

Mod <- nls(Weight ~ A*(Length^B),start = list(A= 1, B=3),main="P. ornatus Length-Weight")
coef(Mod)

lens <- seq(1,150,by=1)

LenWei <- coef(Mod)[1]*(lens^coef(Mod)[2])
lines(lens,LenWei)
legend("topleft",legend=c(paste("A = ",round(coef(Mod)[1],3),sep=""),paste("B = ",round(coef(Mod)[2],3),sep="")))


#Len at min weight
#66.5mm

coef(Mod)[1]*(67^coef(Mod)[2])

hist(Length,breaks=seq(0,160,1))$counts
table(Length)
hist(Length,breaks=seq(0,160,1))$mids


# Selectivity at length: sel50 = 63.5, sel95 = 79.5 
#VB growth
Lt2Age<-function(Linf,k,t0,lt)
{
  age.lt<-(t0-(log(1-(lt/Linf))/k))
  return(age.lt)
}

Sel50<-Lt2Age(150,0.57,0,63.5)
Sel95<-Lt2Age(150,0.57,0,79.5)
MinSize <- Lt2Age(150,0.57,0,66.5)

#First reproduction at 75mm, McFarlane and Moore 1986
Lt2Age(150,0.57,0,75)
