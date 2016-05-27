#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: record.R                                                      ##
## Contains: record                                                    ##
##                                                                     ##
## Written by Marcelo Mollinari                                        ##
## copyright (c) 2007-9, Marcelo Mollinari                             ##
##                                                                     ##
## First version: 11/29/2009                                           ##
## Last update: 12/2015                                                ##
## The detailed description of the algorithm was removed by Augusto    ##
## Garcia, on 2015/07/25, since it was not compiling in new R versions ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################



##' Recombination Counting and Ordering
##'
##' Implements the marker ordering algorithm \emph{Recombination Counting and
##' Ordering} (\cite{Van Os et al., 2005}).
##'
##' \emph{Recombination Counting and Ordering} (\emph{RECORD}) is an algorithm
##' for marker ordering in linkage groups. It is not an exhaustive search
##' method and, therefore, is not computationally intensive. However, it does
##' not guarantee that the best order is always found. The only requirement is
##' a matrix with recombination fractions between markers.
##'
##' After determining the order with \emph{RECORD}, the final map is
##' constructed using the multipoint approach (function
##' \code{\link[onemap]{map}}).
##'
##' @param input.seq an object of class \code{sequence}.
##' @param times integer. Number of replicates of the RECORD procedure.
##' @param LOD minimum LOD-Score threshold used when constructing the pairwise
##' recombination fraction matrix.
##' @param max.rf maximum recombination fraction threshold used as the LOD
##' value above.
##' @param tol tolerance for the C routine, i.e., the value used to evaluate
##' convergence.
##' @return An object of class \code{sequence}, which is a list containing the
##' following components: \item{seq.num}{a \code{vector} containing the
##' (ordered) indices of markers in the sequence, according to the input file.}
##' \item{seq.phases}{a \code{vector} with the linkage phases between markers
##' in the sequence, in corresponding positions. \code{-1} means that there are
##' no defined linkage phases.} \item{seq.rf}{a \code{vector} with the
##' recombination frequencies between markers in the sequence. \code{-1} means
##' that there are no estimated recombination frequencies.}
##' \item{seq.like}{log-likelihood of the corresponding linkage map.}
##' \item{data.name}{name of the object of class \code{outcross} with the raw
##' data.} \item{twopt}{name of the object of class \code{rf.2pts} with the
##' 2-point analyses.}
##' @author Marcelo Mollinari, \email{mmollina@@usp.br}
##' @seealso \code{\link[onemap]{make.seq}} and \code{\link[onemap]{map}}
##' @references Mollinari, M., Margarido, G. R. A., Vencovsky, R. and Garcia,
##' A. A. F. (2009) Evaluation of algorithms used to order markers on genetics
##' maps. \emph{Heredity} 103: 494-502.
##'
##' Van Os, H., Stam, P., Visser, R.G.F. and Van Eck, H.J. (2005) RECORD: a
##' novel method for ordering loci on a genetic linkage map. \emph{Theoretical
##' and Applied Genetics} 112: 30-40.
##' @keywords utilities
##' @examples
##'
##' \dontrun{
##'   ##outcross example
##'   data(example.out)
##'   twopt <- rf.2pts(example.out)
##'   all.mark <- make.seq(twopt,"all")
##'   groups <- group(all.mark)
##'   LG1 <- make.seq(groups,1)
##'   LG1.rec <- record(LG1)
##'
##'   ##F2 example
##'   data(fake.f2.onemap)
##'   twopt <- rf.2pts(fake.f2.onemap)
##'   all.mark <- make.seq(twopt,"all")
##'   groups <- group(all.mark)
##'   LG1 <- make.seq(groups,1)
##'   LG1.rec <- record(LG1)
##'   LG1.rec
##' }
##'
split_map_batches <- function(x, cores){
  start <- 1
  interval <- ceiling(length(x) / cores)
  while(interval < 80){
    cores = cores - 1
    interval <- ceiling(length(x) / cores)
    if(cores <= 0) stop("Mapping group < 80. No need to apply parallel record")
  }
  end <- interval
  res <- list()
  for(f in 1:cores){
    if(f == cores){
      if(length(start:end) < 80){
        res[[f - 1]] <- c(res[[f - 1]], x[(start + 1):end])
      } else {
        res[[f]] <- x[start:end]
      }
    } else {
      res[[f]] <- x[start:end]
      start <- end
      end <- end + interval
      if(end > length(x)){ end <- length(x)}
    }
  }
  return(res)
}

combine_map_batches <- function(x){
  final <- NULL
  for(f in 1:length(x)){
    if(f == 1){
      final$seq.num <- x[[f]]$seq.num
      final$seq.phases <- x[[f]]$seq.phases
      final$seq.rf <- x[[f]]$seq.rf
      final$seq.like<- x[[f]]$seq.like
    } else {
      final$seq.num <- c(final$seq.num, x[[f]]$seq.num[-1])
      final$seq.phases <- c(final$seq.phases, x[[f]]$seq.phases)
      final$seq.rf <- c(final$seq.rf, x[[f]]$seq.rf)
      final$seq.like<- final$seq.like + x[[f]]$seq.like
    }
    final$seq.like <- final$seq.like / length(x)
    final$data.name <- x[[1]]$data.name
    final$twopt <- x[[1]]$twopt
    class(final)<-"sequence"
  }
  return(final)
}

record.parallel.rfc<-function(input.seq, times=10, cores=10, LOD=0, max.rf=0.5,
                              tol=10E-5, useC = TRUE, domap = TRUE){
  require(parallel)
  require(RcppArmadillo)
  ## checking for correct object
  if(!any(class(input.seq)=="sequence")) stop(deparse(substitute(input.seq))," is
                                              not an object of class 'sequence'")
  n.mrk <- length(input.seq$seq.num)

  ## create reconmbination fraction matrix

  if(class(get(input.seq$twopt))[2]=="outcross")
    r<-get_mat_rf_out(input.seq, LOD=FALSE, max.rf=max.rf, min.LOD=LOD, useC = TRUE)
  else
    r<-get_mat_rf_in(input.seq, LOD=FALSE, max.rf=max.rf, min.LOD=LOD)
  r[is.na(r)]<-0.5
  diag(r)<-0

  ##RECORD algorithm
  X<-r*get(input.seq$data.name, pos=1)$n.ind ## Obtaining X multiplying the MLE of the recombination
  ## fraction by the number of individuals

  if(n.mrk==2)
    return(map(make.seq(get(input.seq$twopt),input.seq$seq.num[1:2],twopt=input.seq$twopt), tol=10E-5))

  ## For three markers (calculation of 3 possible orders: 3!/2)
  else if(n.mrk==3) {
    all.perm<-perm.pars(1:3)
    m.old<-Inf
    for(k in 1:nrow(all.perm)){
      m.new<-COUNT(X, all.perm[k,])
      if(m.new < m.old)
        result.new <- all.perm[k,]; m.old<-m.new
    }
  }

  ## For more than three markers (RECORD algorithm itself)
  else{
    marks<-c(1:n.mrk)
    result.new<-sample(marks)## randomize markers
    results.list <- mclapply(1:times,mc.cores=cores, function(l){ ## loop for replicates the RECORD procedure
      result<-sample(marks, 2)
      for(i in 2:(n.mrk-2)){
        next.mark<-sample(c(1:n.mrk)[-result],1)
        partial<-rep(NA,2*length(result)+1)
        partial[seq(2,(length(partial)-1), by = 2)]<-result
        temp<-10e1000
        for(j in seq(1,(length(partial)), by = 2)){
          partial.temp<-partial
          partial.temp[j]<-next.mark
          temp.new<-CCOUNT(X, partial.temp[is.na(partial.temp)==FALSE] - 1)
          if(temp.new<temp) {
            result<-partial.temp[is.na(partial.temp)==FALSE]
            temp<-temp.new
          }
        }
      }
      next.mark<-c(1:n.mrk)[-result]
      partial<-rep(NA,2*length(result)+1)
      partial[seq(2,(length(partial)-1), by = 2)]<-result
      temp<-10e1000
      for(j in seq(1,(length(partial)), by = 2)){
        partial.temp<-partial
        partial.temp[j]<-next.mark
        temp.new<-CCOUNT(X, partial.temp[ is.na(partial.temp)==FALSE] - 1)
        if(temp.new<temp) {
          result<-partial.temp[ is.na(partial.temp)==FALSE]
          temp<-temp.new
        }
      }
      for(i in 1:(n.mrk-1)){
        for(j in 1:(n.mrk-i)){
          y<-result[j:(j+i)]
          if(j==1){
            if(i!=n.mrk-1){
              perm.order<-c(rev(result[1:(1+i)]),result[(i+2):length(result)])
              COUNT.temp<-CCOUNT(X,perm.order - 1)
              if(COUNT.temp < temp) {
                result<-perm.order
                temp<-COUNT.temp
              }
            }
          }
          else{
            if(j!=(n.mrk-i)){
              perm.order<-c(result[1:(j-1)],rev(result[j:(j+i)]),result[(j+i+1):length(result)])
              COUNT.temp<-CCOUNT(X,perm.order - 1)
              if(COUNT.temp < temp){
                result<-perm.order
                temp<-COUNT.temp
              }
            }
            else{
              perm.order<-c(result[1:(j-1)],rev(result[j:length(result)]))
              COUNT.temp<-CCOUNT(X,perm.order - 1)
              if(COUNT.temp < temp){
                result<-perm.order
                temp<-COUNT.temp
              }
            }
          }
        }
      }
      return(result)
    })
    count.scores <- mclapply(results.list, mc.cores = cores, function(f){
      CCOUNT(X,f - 1)
    })
    min_count <- which.min(unlist(count.scores))
    result.new <- results.list[[min_count]]
  }
  ## end of RECORD algorithm
  if(domap){
    cat("\norder obtained using RECORD algorithm:\n\n", input.seq$seq.num[result.new],".\n\n",
        "now calling map()\n\n")
    batch <- split_map_batches(result.new, cores)
    m <- mclapply(batch, mc.cores = length(batch), function(x){
      s <- make.seq(get(input.seq$twopt),
                    input.seq$seq.num[x],
                    twopt=input.seq$twopt)
      map(s)
    })
    final <- combine_map_batches(m)
    return(final)
  } else {
    cat("\norder obtained using RECORD algorithm:\n\n", input.seq$seq.num[result.new],".\n\n",
        "NOT calling map()\n\n")
    return(result.new)
  }
}

##end of file
