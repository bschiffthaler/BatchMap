\name{onemap-internal}
\alias{acum}
\alias{adjust.rf.ril}
\alias{avoid.reverse}
\alias{check.linkage}
\alias{codif.data}
\alias{comb}
\alias{comb.ger}
\alias{compare_inbred}
\alias{compare_outcross}
\alias{combine_map_batches}
\alias{Cindex}
\alias{diplo}
\alias{draw.order}
\alias{draw.rf.inter}
\alias{draw.try}
\alias{est_map_hmm_f2}
\alias{est_map_hmm_bc}
\alias{est_rf_bc}
\alias{est_rf_f2}
\alias{est_rf_out}
\alias{generate_all}
\alias{generate_one}
\alias{generate_overlapping_batches}
\alias{generate_rand}
\alias{get.bins}
\alias{get_mat_rf_in}
\alias{get_vec_rf_in}
\alias{get_mat_rf_out}
\alias{get_twopt_info}
\alias{get_vec_rf_out}
\alias{perm.pars}
\alias{perm.tot}
\alias{phases}
\alias{plotFunction.out}
\alias{rem.amb.ph}
\alias{return.geno}
\alias{ripple_seq_inbred}
\alias{ripple_seq_outcross}
\alias{ser.ord}
\alias{split_map_batches}
\alias{try_seq_inbred}
\alias{try_seq_outcross}
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
codif.data(geno.in,segr.type.in)
comb(x,y)
comb.ger(f)
combine_map_batches(x)
compare_inbred(input.seq, n.best, tol, verbose = FALSE)
compare_outcross(input.seq, n.best, tol, verbose)
Cindex(order,r)
diplo(w, seq.num, seq.phases)
draw.order(map.input)
draw.rf.inter(info.graph, scale, cex)
draw.try(base.input, try.input, pos=NULL)
est_rf_out(geno, mrk, seg_type, nind, verbose)
est_rf_f2(geno, mrk, seg_type, nind, verbose)
est_rf_bc(geno, mrk,  nind, type, verbose)
est_map_hmm_f2(geno, rf.vec, verbose, tol)
est_map_hmm_bc(geno, rf.vec, verbose, tol)
generate_all(input.seq, p, ws)
generate_one(input.seq, p, ws, no_reverse)
generate_overlapping_batches(input.seq, size = 50, overlap = 15, silent = FALSE)
generate_rand(input.seq, p, ws, n, pref)
get.bins(geno, exact)
get_mat_rf_in(input.seq, LOD, max.rf, min.LOD)
get_twopt_info(twopt, small, big)
get_vec_rf_in(input.seq, LOD, acum)
get_mat_rf_out(input.seq, LOD, max.rf, min.LOD)
get_vec_rf_out(input.seq, LOD, max.rf, min.LOD, acum)
perm.pars(v)
perm.tot(v)
phases(input.seq, LOD=0, max.rf=0.50)
plotFunction.out(info.graph, cex, main)
predict_time(batches, ws, timings)
rem.amb.ph(M,w,seq.num)
return.geno(segr.type, link.phases)
ripple_seq_inbred(input.seq, ws, ext.w, LOD, tol)
ripple_seq_outcross(input.seq, ws, LOD, tol)
ser.ord(r,i)
split_map_batches(x, cores)
try_seq_inbred(input.seq, mrk, tol, draw.try, pos, verbose)
try_seq_outcross(input.seq, mrk, tol, draw.try, pos, verbose)
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
