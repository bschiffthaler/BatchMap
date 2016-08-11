/*
  OneMap: software for genetic mapping in outcrossing species
  Copyright (C) 2007-2015 Gabriel R A Margarido and Marcelo Mollinari

  This file is part of OneMap.

  OneMap is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  OneMap is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
*/

/*
  File: twopts_bc.cpp

  Description: Function to compute the recombination fraction in
  Backcross experimental populations.

  Written by Marcelo Mollinari

  Escola Superior de Agricultura "Luiz de Queiroz"
  Departamento de Genética - São Paulo, Brazil
  Contact: mmollina@usp.br
  First version: 09/2015
  Last update: 11/2015
*/

#include <Rcpp.h>
#include "twopts_bc.h"
#include "progressBar.h"
#define TOL 0.00001
#define rf_TOL_min 1e-50
#define rf_TOL_max 0.5-1e-50
#define LN_5 log(1.0/2.0)
using namespace Rcpp;
using namespace std;


RcppExport SEXP est_rf_bc_wrap(SEXP geno_R, SEXP mrk_R, SEXP n_ind_R, SEXP type_R, SEXP verbose_R)
{
  Rcpp::NumericVector geno = Rcpp::as<Rcpp::NumericVector>(geno_R);
  long n_ind = Rcpp::as<long>(n_ind_R);
  long type = Rcpp::as<long>(type_R);
  bool verbose = Rcpp::as<bool>(verbose_R);
  long mrk = Rcpp::as<long>(mrk_R);
  NumericMatrix z = est_rf_bc(geno, mrk, n_ind, type, verbose);
  return(wrap(z));
}


Rcpp::NumericMatrix est_rf_bc(Rcpp::NumericVector geno, long mrk,
			      long n_ind, long type, bool verbose)
{
  long n_mar=((long)geno.size()/n_ind);
  double rtemp, l, l0, mis=0, nr=0;
  Rcpp::NumericMatrix r(n_mar, n_mar);
  NumericMatrix rm(2, n_mar);
  long chunk=((n_mar*n_mar)-n_mar)/20, ct1=0, ct2=1, a1, a2, a3;
  if(mrk < 0)
    {
      if(verbose && n_mar > 100)
	Rcpp::Rcout << "Computing " << ((n_mar*n_mar)-n_mar)/2 << " recombination fractions:\n\n" << "\t0%\t.";
      else if(verbose==1)
	Rcpp::Rcout << "Computing " << ((n_mar*n_mar)-n_mar)/2 << " recombination fractions ... \n";
      a1=0;
      a2=n_mar-1;
    }
  else
    {
      a1=mrk;
      a2=mrk+1;
    }
  for(long i=a1; i < a2; i++)
    {
      if(mrk < 0){
	if(verbose && n_mar > 100)
	  {
	    if(ct1%(chunk/(chunk/10))==0)
	      {
	      progressBar(Rcpp::Rcout, ct1, ((n_mar*n_mar)-n_mar)/2, 80);
	      }
	  }
      }
      if(mrk < 0) a3=(i+1);
      else a3=0;
      R_CheckUserInterrupt(); /* check for ^C */
      for(long j=a3; j  < n_mar; j++)
	{
	  ct1++;
	  nr=mis=0;
	  std::vector<long> k_sub(&geno[i*n_ind],&geno[i*n_ind+n_ind]);
	  std::vector<long> k1_sub(&geno[j*n_ind],&geno[j*n_ind+n_ind]);
	  for(long k=0; k < n_ind; k++)
	    {
	      if(k_sub[k]==0 || k1_sub[k]==0) mis++;
	      else if((k_sub[k] != k1_sub[k])) nr++;
	    }
	  rtemp=nr/(n_ind-mis);
	  if(rtemp < rf_TOL_min) rtemp=rf_TOL_min;
	  if(rtemp > rf_TOL_max)
	    {
	      if(mrk < 0)
		{
		  r(j,i)=rf_TOL_max;
		  r(i,j)=0.0;
		}
	      else
		{
		  rm(0,j)=rf_TOL_max;
		  rm(1,j)=0.0;
		}
	    }
	  else{
	    l=(n_ind-(mis+nr))*log(1-rtemp)+nr*log(rtemp);
	    l0=LN_5*(n_ind-mis);
	    if(mrk < 0){
	      if(type==0)  // backcross
		r(j,i)=rtemp;
	      else if(type==1) // ril self mating
		r(j,i)=rtemp*2/(1+2*rtemp);
	      else if(type==2) // ril sib mating
		r(j,i)=rtemp*4/(1+6*rtemp);
	    r(i,j)=(l-l0)/log(10.0);
	    }
	    else
	      {
		if(type==0)  // backcross
		  rm(0,j)=rtemp;
		else if(type==1) // ril self mating
		  rm(0,j)=rtemp*2/(1+2*rtemp);
		else if(type==2) // ril sib mating
		  rm(0,j)=rtemp*4/(1+6*rtemp);
		rm(1,j)=(l-l0)/log(10.0);
	      }
	  }
	}
    }
  if(mrk < 0){
    if(verbose==1 && n_mar > 100) Rcpp::Rcout << "\t100%\n";
    return(r);
  }
  else
    {
      return(rm);
    }
}
