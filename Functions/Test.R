Test<- function(x)
{
  Blarg<- x^2
  Blarg2<- NULL
  Blarg2$Test<- 1
  b<- x^2
  d<- 'test'
  return(list(b=b,d=d))
}