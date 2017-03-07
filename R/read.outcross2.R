##' Read data from a full-sib progeny (outcrossing populations)
##'
##' This version implements the \code{read.outcross} function in a faster
##' way. Everything else is essentially the same.
##'
##' @param infile the name of the input file which contains the data to be read.
##' @return An object of class \code{outcross}, i.e., a list with the following
##' components: \item{geno}{a matrix with integers indicating the genotypes
##' read for each marker. Each column contains data for a marker and each row
##' represents an individual.} \item{n.ind}{number of individuals.}
##' \item{n.mar}{number of markers.} \item{segr.type}{a vector with the
##' segregation type of each marker, as \code{strings}.} \item{segr.type.num}{a
##' vector with the segregation type of each marker, represented in a
##' simplified manner as integers, i.e. 1 corresponds to markers of type
##' \code{"A"}; 2 corresponds to markers of type \code{"B1.5"}; 3 corresponds
##' to markers of type \code{"B2.6"}; 4 corresponds to markers of type
##' \code{"B3.7"}; 5 corresponds to markers of type \code{"C.8"}; 6 corresponds
##' to markers of type \code{"D1"} and 7 corresponds to markers of type
##' \code{"D2"}} \item{n.phe}{the number of traits included in the file}
##' \item{pheno}{the name of the phenoytpes} \item{input}{the name of the input file.}
##' @author Adapted from Karl Broman (package \pkg{qtl}) by Gabriel R A
##' Margarido, \email{gramarga@@gmail.com}, later with additions from Luciano C Silva
##' @seealso \code{example} directory in the package source.
##' @references Broman, K. W., Wu, H., Churchill, G., Sen, S., Yandell, B.
##' (2008) \emph{qtl: Tools for analyzing QTL experiments} R package version
##' 1.09-43
##'
##' Lincoln, S. E., Daly, M. J. and Lander, E. S. (1993) Constructing genetic
##' linkage maps with MAPMAKER/EXP Version 3.0: a tutorial and reference
##' manual. \emph{A Whitehead Institute for Biomedical Research Technical
##' Report}.
##'
##' Wu, R., Ma, C.-X., Painter, I. and Zeng, Z.-B. (2002) Simultaneous maximum
##' likelihood estimation of linkage and linkage phases in outcrossing species.
##' \emph{Theoretical Population Biology} 61: 349-363.
##' @keywords IO
##' @examples
##' \dontrun{
##'     outcr_data <-
##' read.outcross2("data_file.txt")
##'   }
read.outcross2 <- function(infile)
{
  X <- READ_OUTCROSS(path.expand(infile))
  m <- matrix(X$geno,ncol = X$n.mar)
  colnames(m) <- X$marker
  l <- list(geno = m, n.ind = X$n.ind, n.mar = X$n.mar,
            segr.type = X$segr.type, segr.type.num = X$segr.type.num,
            n.phe = 0, pheno = numeric(0), input = infile)
  attr(l,"class") <- "outcross"
  return(l)
}
