#######################################################################
#                                                                     #
# Package: BatchMap                                                     #
#                                                                     #
# File: cpp.utils.R                                                   #
# Contains: get.bins                                                  #
# These functions are for internal use only                           #
#                                                                     #
# Written Marcelo Mollinari                                           #
#                                                                     #
# First version: 09/2015                                              #
# Last update: 01/2016                                                #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

# This function calls C++ routine to find markers with redundant information
get.bins <- function(geno, exact=TRUE)
{
  bins<-.Call("get_bins",
              geno,
              as.numeric(exact),
              options()$width-6,
              PACKAGE = "BatchMap" )
  return(bins)
}

# This function calls C++ routine for two-point analysis (outcross)
est_rf_out<-function(geno, mrk=0, seg_type=NULL, nind, verbose=TRUE)
{
  r<-.Call("est_rf_out_wrap",
           geno,
           mrk=mrk-1,
           as.numeric(seg_type),
           as.numeric(nind),
           as.numeric(verbose),
           PACKAGE = "BatchMap" )

  if(mrk <= 0)
  {
      names(r)<-c("CC", "CR", "RC", "RR")
      for(i in 1:4) dimnames(r[[i]])<-list(colnames(geno), colnames(geno))
      return(r)
  }
  else
  {
      rownames(r[[1]])<-c("rCC", "rCR", "rRC", "rRR")
      colnames(r[[1]])<-colnames(geno)
      rownames(r[[2]])<-c("LODCC", "LODCR", "LODRC", "LODRR")
      colnames(r[[2]])<-colnames(geno)
      return(r)
  }
}

##' C++ routine for multipoint analysis in outcrossing populations
##'
##' It calls C++ routine that implements the methodology of Hidden
##' Markov Models (HMM) to construct multipoint linkage maps in
##' outcrossing species
##'
##' @param geno matrix of genotypes. Rows represent marker and columns
##'     represent individuals.
##'
##' @param type a vector indicating the type of marker. For more
##'     information see \code{\link[onemap]{read.onemap}}
##'
##' @param phase a vector indicating the linkage phases between
##'     markers. For more information see
##'     \code{\link[onemap]{make.seq}}
##'
##' @param rf.vec a vector containing the recombination fraction
##'     initial values
##'
##' @param verbose If \code{TRUE}, print tracing information.
##'
##' @param tol tolerance for the C routine, i.e., the value used to
##'     evaluate convergence.
##'
##' @return a list containing the re-estimated vector of recombination
##'     fractions and the logarithm of the likelihood
##'
##' @keywords internal
##'
##' @export
est_map_hmm_out<-function(geno, type,  phase, rf.vec=NULL, verbose=TRUE, tol=1e-6)
{
  if(is.null(rf.vec))
    rf.vec<-rep(0.1, (nrow(geno)-1))
  r<-.Call("est_hmm_out",
           geno,
           as.numeric(type),
           as.numeric(phase),
           as.numeric(rf.vec),
           as.numeric(verbose),
           as.numeric(tol),
           PACKAGE = "BatchMap")
  names(r)<-c("rf", "loglike")
  return(r)
}
#end of the file

