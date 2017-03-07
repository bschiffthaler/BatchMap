#######################################################################
##                                                                     ##
## Package: BatchMap                                                     ##
##                                                                     ##
## File: pseudo.testcross.split.R                                      ##
## Contains: pseudo.testcross.split                                    ##
##                                                                     ##
## Written by Bastian Schiffthaler                                     ##
## copyright (c) 2017 Bastian Schiffthaler                             ##
##                                                                     ##
##                                                                     ##
## First version: 07/03/2017                                           ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#######################################################################

##' Split a dataset into parent-specific subsets
##'
##' In order to create a map using the pseudo-testcross approach, the
##' data needs to be split into parent specific subsets. This function
##' uses the segregation type to create subset given an input object
##' of type \code{group}.
##'
##' @param input.obj an object of class \code{group}
##'
##' @return A list with two marker sequences for each linkage group. One for
##' sequences with type "D2.15", one for those with type "D1.10".
##'
##' @export
pseudo.testcross.split <- function(input.obj)
{
  ls <- lapply(seq(1,input.obj$n.groups), function(f){
    gr <- which(input.obj$groups == f)
    out <- get(input.obj$data.name)
    st <- out$segr.type[gr]
    p1 <- make.seq(get(input.obj$twopt),gr[st != "D2.15"],
                   twopt = input.obj$twopt)
    p2 <- make.seq(get(input.obj$twopt),gr[st != "D1.10"],
                   twopt = input.obj$twopt)
    list("d1.10" = p1, "d2.15" = p2)
  })
  names(ls) <- paste0("LG",seq(1,input.obj$n.groups))
  return(unlist(ls,recursive = FALSE))
}
