## better to use rlm or lqa
bicor.test <- function(x, y) {
  b = bicor(x,y)
  n = min(length(na.omit(x)), length(na.omit(y)))
  stat = abs(b * sqrt( (n-2)/(1-(b^2))))
  p.value = 1 - pt( stat, df=n-2 )
  pr = cor.test(x, y)
  r = pr$estimate
  res = c(bicor=b, b.t=stat, 
          r, pr$statistic, 
          b.p.value=p.value, r.p.value=pr$p.value)
  return(res)
}
