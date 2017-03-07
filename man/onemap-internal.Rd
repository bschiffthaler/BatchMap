\name{onemap-internal}
\alias{acum}
\alias{adjust.rf.ril}
\alias{avoid.reverse}
\alias{CCOUNT}
\alias{check.linkage}
\alias{codif.data}
\alias{comb}
\alias{comb.ger}
\alias{compare.outcross}
\alias{Cindex}
\alias{diplo}
\alias{est_rf_out}
\alias{flip_phases}
\alias{generate_all}
\alias{generate_one}
\alias{generate.overlapping.batches}
\alias{generate_rand}
\alias{get.bins}
\alias{get_mat_rf_in}
\alias{get_vec_rf_in}
\alias{get_mat_rf_out}
\alias{GET_RF_MAT_NO_LOD}
\alias{get_twopt_info}
\alias{get_vec_rf_out}
\alias{perm.pars}
\alias{perm.tot}
\alias{phases}
\alias{READ_OUTCROSS}
\alias{rem.amb.ph}
\alias{return.geno}
\alias{ser.ord}
\alias{split_map_batches}
\alias{try.seq.outcross}
\alias{predict_time}
\alias{print.compare}
\alias{print.onemap}
\alias{print.onemap.bin}
\alias{print.order}
\alias{print.sequence}
\title{Internal onemap functions}
\description{
  Functions from \pkg{onemap} not to be directly called by the user.
}
\usage{
acum(w)
avoid.reverse(input)
adjust.rf.ril(r,type,expand)
check.linkage(i, s, cl, geno, st, max.rf, LOD)
codif.data(geno.in,segr.type.in,cross = c("outcross"))
comb(x,y)
comb.ger(f)
compare.outcross(input.seq, n.best, tol, verbose)
Cindex(order,r)
diplo(w, seq.num, seq.phases)
est_rf_out(geno, mrk, seg_type, nind, verbose)
generate_all(input.seq, p, ws)
generate_one(input.seq, p, ws, no_reverse)
generate.overlapping.batches(input.seq, size = 50, overlap = 15, silent = FALSE)
generate_rand(input.seq, p, ws, n, pref)
get.bins(geno, exact)
get_twopt_info(twopt, small, big)
get_mat_rf_out(input.seq, LOD = FALSE, max.rf = 0.5, min.LOD = 0, useC = FALSE)
get_vec_rf_out(input.seq, LOD, max.rf, min.LOD, acum)
perm.pars(v)
perm.tot(v)
phases(input.seq, LOD=0, max.rf=0.50)
predict_time(batches, ws, timings)
rem.amb.ph(M,w,seq.num)
return.geno(segr.type, link.phases)
ser.ord(r,i)
try.seq.outcross(input.seq, mrk, tol, pos, verbose)
\method{print}{compare}(x, \dots)
\method{print}{onemap}(x, \dots)
\method{print}{onemap.bin}(x, \dots)
\method{print}{order}(x, \dots)
\method{print}{order}(x, \dots)
\method{print}{outcross}(x, \dots)
\method{print}{sequence}(x, \dots)
}
\author{Gabriel Rodrigues Alves Margarido and Marcelo Mollinari}
\keyword{internal}
