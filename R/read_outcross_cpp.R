read_outcross_cpp <- function(infile)
{
  X <- READ_OUTCROSS(infile)
  m <- matrix(X$geno,ncol = X$n.mar)
  colnames(m) <- X$marker
  l <- list(geno = m, n.ind = X$n.ind, n.mar = X$n.mar,
            segr.type = X$segr.type, segr.type.num = X$segr.type.num,
            n.phe = 0, pheno = numeric(0), input = infile)
  attr(l,"class") <- "outcross"
  return(l)
}
