fast.cor.FS = function(x, y, min_su){
  
  require(FCBF)
  
  dis = discretize_exprs(t(x))
  
  fcbf = fcbf(dis, y, verbose = T, minimum_su = min_su)
  res = x[,fcbf$index]
  
  return(res)
}
