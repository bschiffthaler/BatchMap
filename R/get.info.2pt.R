#######################################################################
#                                                                     #
# Package: BatchMap                                                     #
#                                                                     #
# File: get_rf.R                                                      #
# Contains: get_mat_rf_in, get_vec_rf_in, get_mat_rf_out,             #
#           get_vec_rf_out                                            #
#                                                                     #
# Written by Marcelo Mollinari                                        #
# copyright (c) 2015, MarceloMollinari                                #
#                                                                     #
# First version: 12/2015                                              #
# Last update: 12/2015                                                #
# License: GNU General Public License version 2 (June, 1991) or later #
#                                                                     #
#######################################################################

#For a guiven sequence, this function gets the recombination
#fraction/LOD matrix for outcrossing
get_mat_rf_out<- function(input.seq, LOD=FALSE, max.rf=0.5, min.LOD=0, useC = FALSE) {
     if(!any(class(input.seq)=="sequence")) stop(deparse(substitute(input.seq))," is not an object of class 'sequnece'")
     if(length(input.seq$seq.num) < 2) stop("The sequence must have at least 2 markers")
     n.mrk<-length(input.seq$seq.num)
     mrk.names <- colnames(get(input.seq$data.name, pos=1)$geno)[input.seq$seq.num]
     ## create reconmbination fraction matrix
     if(LOD){
       r <- matrix(NA,n.mrk,n.mrk)
       dimnames(r)<-list(mrk.names, mrk.names)
       for(i in 1:(n.mrk-1)) {
         for(j in (i+1):n.mrk) {
           k<-sort(c(input.seq$seq.num[i], input.seq$seq.num[j]))
           rfs<-sapply(get(input.seq$twopt)$analysis, function(x,i,j) x[i,j], k[2], k[1])
           LODs<-sapply(get(input.seq$twopt)$analysis, function(x,i,j) x[i,j], k[1], k[2])
           ## check if any assignment meets the criteria
           phases <- which((LODs >= min.LOD) & rfs <= max.rf)
           if(length(phases) == 0)
           {
             r[i,j] <- NA
             r[j,i] <- NA
           }
           else
           {
             r.temp<-rfs[phases[which.max(LODs[phases])]]
             if(r.temp > 0.5) r.temp<-0.5
             r[i,j]<-r.temp
             r[j,i]<-max(LODs[phases])
           }
         }
       }
     }
     else
     {
       if(useC){
         r <- GET_RF_MAT_NO_LOD(input.seq$seq.num, n.mrk, get(input.seq$twopt)$analysis$CC,
                                get(input.seq$twopt)$analysis$CR, get(input.seq$twopt)$analysis$RC,
                                get(input.seq$twopt)$analysis$RR, min.LOD, max.rf)
         dimnames(r)<-list(mrk.names, mrk.names)
       } else {
         r <- matrix(NA,n.mrk,n.mrk)
         dimnames(r)<-list(mrk.names, mrk.names)
         for(i in 1:(n.mrk-1)) {
           for(j in (i+1):n.mrk) {
             k<-sort(c(input.seq$seq.num[i], input.seq$seq.num[j]))
             rfs<-sapply(get(input.seq$twopt)$analysis, function(x,i,j) x[i,j], k[2], k[1])
             LODs<-sapply(get(input.seq$twopt)$analysis, function(x,i,j) x[i,j], k[1], k[2])
             ## check if any assignment meets the criteria
             phases <- which((LODs >= min.LOD) & rfs <= max.rf)
             if(length(phases) == 0)
             {
               r[j,i] <- r[i,j] <- NA
             }
             else
             {
               r.temp<-rfs[phases[which.max(LODs[phases])]]
               if(r.temp > 0.5) r.temp<-0.5
               r[j,i]<-r[i,j]<-r.temp
             }
           }
         }
       }
     }
     return(r)
 }

#For a guiven sequence, this function gets the recombination
#fraction/LOD matrix for outcrossing
get_vec_rf_out<- function(input.seq, LOD=FALSE, max.rf=0.5, min.LOD=0, acum=TRUE)
{
    if(!any(class(input.seq)=="sequence"))
        stop(deparse(substitute(input.seq))," is not an object of class 'sequnece'")
    if(length(input.seq$seq.num) < 2) stop("The sequence must have at least 2 markers")
    mat<-get_mat_rf_out(input.seq=input.seq, LOD=LOD, max.rf=max.rf, min.LOD=min.LOD)
    mrk.names<-colnames(mat)
    r<-numeric(length(input.seq$seq.num)-1)
    pair.names<-character(ncol(mat)-1)
    if(LOD) mat<-t(mat)
    for(i in 1:(length(input.seq$seq.num)-1))
    {
        r[i]<-mat[i, i+1]
        pair.names[i]<-paste(mrk.names[i], mrk.names[i+1], sep="-")
    }
    if(acum)
    {
        r<-c(0,cumsum(get(get(".map.fun", envir=.onemapEnv))(r)))
        names(r)<-mrk.names
    }
    else
        names(r)<-pair.names
    return(r)
}


# end of file
