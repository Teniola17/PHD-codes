BinaryExactCI <- function(
  x,
  n,
  alpha=0.05)
{

alpha2 <- alpha/2

Low = alpha/2
High = 1-alpha/2

pLow = qbeta( Low, x+(x==0), n-x+1)
pHigh = qbeta( High, x+1, n-x+((n-x)==0))

nms = c( paste(round(100*Low,1),"%",sep=""),paste(round(100*High,1),"%",sep="") )

 CI = cbind(pLow,pHigh)
 colnames(CI) = nms


return( CI )
}


