//////////////////////////////////////////////////////////////////////////
// MCMCmetropSO.cc is based on MCMCmetrop1R and samples from partition 
// space to generate a sample of a user-defined objective function over
// districting plans. 
//
// Santiago Olivella
// Dept. of Political Science
// Washington University in St. Louis
// olivella@wustl.edu
// 
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
// file for more information.
//

// Copyright for MCMCmetrop1R (C) 2003-2007 Andrew D. Martin and Kevin 
// M. Quinn
// Copyright for MCMCmetrop1R (C) 2007-present Andrew D. Martin, Kevin 
// M. Quinn, and Jong Hee Park
//////////////////////////////////////////////////////////////////////////

#ifndef MCMCMETROPSO_CC
#define MCMCMETROPSO_CC

#include "matrix.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "MCMCrng.h"
#include "MCMCfcds.h"

#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

#include <Rdefines.h>
#include <Rinternals.h>


using namespace std;
using namespace scythe;

//function to get candidate locations without replacement

double* gen_cands(SEXP candFun, SEXP currCuts, SEXP sampenv){
  SEXP R_fcall;
  int nCuts = length(currCuts);
  PROTECT(R_fcall = lang2(candFun,R_NilValue));
  SETCADR(R_fcall,currCuts);
  SEXP funval;
  PROTECT(funval = eval(R_fcall,sampenv));
  if (!isReal(funval)) error("`genFun' must return a double");
  double* candidates = new double [nCuts];
  for (unsigned int i = 0; i < nCuts; ++i) {
    candidates[i] = REAL(funval)[i];
  }
  UNPROTECT(2);
  return candidates;
}


// to evaluate the districting function
double user_fun_dist_eval(SEXP fun, SEXP theta, SEXP myframe) {
  SEXP R_fcall;
  if(!isFunction(fun)) error("`fun' must be a function");
  if(!isEnvironment(myframe)) error("myframe must be an environment");
  PROTECT(R_fcall = lang2(fun, R_NilValue));
  SETCADR(R_fcall, theta);
  SEXP funval;
  PROTECT(funval = eval(R_fcall, myframe));
  if (!isReal(funval)) error("`fun' must return a double");
  double fv = REAL(funval)[0];
  if (fv == R_PosInf) error("`fun' returned +Inf");
  if (R_IsNaN(fv) || R_IsNA(fv)) error("`fun' returned NaN or NA");
  UNPROTECT(2);
  return fv;
}

template <typename RNGTYPE>
void MCMCmetropSO_impl (rng<RNGTYPE>& stream, SEXP& fun, 
                        SEXP& theta, SEXP& myframe,
						SEXP& candFun, 
						SEXP& sampenv,
                        unsigned int burnin,
                        unsigned int mcmc, unsigned int thin,
                        unsigned int verbose,
                        SEXP& sample_SEXP)
{
  // define constants
  const unsigned int npar = length(theta);
  const unsigned int tot_iter = burnin + mcmc;
  const unsigned int nsamp = mcmc / thin;
  
  // initialize matrix to hold the sample
  Matrix<> sample(nsamp, npar+1, false);

  // put theta into a Scythe Matrix 
  if (!isReal(theta)) error("`theta' must be a double");
  double* theta_data = REAL(theta);
  const int theta_nr = length(theta);
  const int theta_nc = 1;
  Matrix <> theta_M (theta_nc, theta_nr, theta_data);
  theta_M = t(theta_M);

  // evaluate userfun at starting value
  double userfun_cur =  user_fun_dist_eval(fun, theta, myframe);
  //if (! logfun) 
    userfun_cur = std::log(userfun_cur);
  
  
  // THE METROPOLIS SAMPLING
  unsigned int count = 0;
  unsigned int accepts = 0;
  for (unsigned int iter = 0; iter < tot_iter; ++iter) {

    // generate candidate value of theta 
    double* theta_canSO = gen_cands(candFun,theta,sampenv);
    Matrix <> theta_can_M (theta_nc, theta_nr, theta_canSO);
    theta_can_M = t(theta_can_M);
    
    // put theta_can_M into a SEXP
    SEXP theta_can;
    PROTECT(theta_can = allocVector(REALSXP, npar));
    for (unsigned int i = 0; i < npar; ++i) {
      //if (!isReal(theta_can_M(i))) error("`theta_can_M' must be a double");
      REAL(theta_can)[i] = theta_can_M(i);
    }
    
    // evaluate user function fun at candidate theta
    double userfun_can = user_fun_dist_eval(fun, theta_can, myframe);
    //if (! logfun)
      userfun_can = std::log(userfun_can);
    const double ratio = std::exp(userfun_can - userfun_cur);
    
    if (stream() < ratio) {
      for (unsigned int i = 0; i < npar; ++i) {
	REAL(theta)[i] = theta_can_M(i);
      }      
      //theta = theta_can;
      theta_M = theta_can_M;
      userfun_cur = userfun_can;
      ++accepts;
    }
    UNPROTECT(1);      

    // store values in matrices
    if ((iter%thin) == 0 && iter >= burnin) {
      for (unsigned int j = 0; j <= npar; j++){
	if (j == npar){
	  sample(count, j) = userfun_cur;
	}
	else{
	  //if (!isReal(theta[j])) error("`theta j' must be a double");
	  sample(count, j) = REAL(theta)[j];
	}
      }
      ++count;
    }
          
    if (verbose && iter % verbose == 0) {
      Rprintf("MCMCmetropSO iteration %i of %i \n", (iter+1), tot_iter);
      Rprintf("function value = %10.5f\n", userfun_cur);
      Rprintf("theta = ");
      for (unsigned int i = 0; i < npar; ++i)
        Rprintf("%4.0f", REAL(theta)[i]);
      Rprintf("\nMetropolis acceptance rate = %3.5f\n\n", 
          static_cast<double>(accepts) / static_cast<double>(iter+1));	
    } 

    
    R_CheckUserInterrupt(); // allow user interrupts
  }

  // put the sample into a SEXP and return it   
  //sample_SEXP = PROTECT(allocMatrix(REALSXP, nsamp, npar));
  for (unsigned int i = 0; i < nsamp; ++i) {
    for (unsigned int j = 0; j <= npar; ++j) {
      REAL(sample_SEXP)[i + nsamp * j] = sample(i,j);
    }
  }
  //UNPROTECT(1);


  // print the the acceptance rate to the console in a way that 
  // everyone (even Windows users) can see
  Rprintf("\n\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  Rprintf("The Metropolis acceptance rate was %3.5f", 
    static_cast<double>(accepts) / static_cast<double>(tot_iter));
  Rprintf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
}


extern "C" {

  // the function that actually does the sampling and returns a value to R
  SEXP MCMCmetropSO_cc(SEXP fun, SEXP theta, SEXP myframe, 
		       SEXP candFun, 
		       SEXP sampenv,
		       SEXP burnin_R,
		       SEXP mcmc_R, SEXP thin_R, 
		       SEXP verbose, SEXP lecuyer_R, SEXP seedarray_R,
		       SEXP lecuyerstream_R)
  {
    // put rng stuff together
    int seedarray[6];
    for(int i=0; i<6; ++i) seedarray[i] = INTEGER(seedarray_R)[i];
    int uselecuyer_cc = INTEGER(lecuyer_R)[0];
    int lecuyerstream_cc = INTEGER(lecuyerstream_R)[0];
    int* uselecuyer = &uselecuyer_cc;
    int* lecuyerstream = &lecuyerstream_cc;

    const unsigned int npar = length(theta);
    const unsigned int nsamp = INTEGER(mcmc_R)[0] / INTEGER(thin_R)[0];
    SEXP sample_SEXP;
    PROTECT(sample_SEXP = allocMatrix(REALSXP, nsamp, npar+1));
    MCMCPACK_PASSRNG2MODEL(MCMCmetropSO_impl, fun, theta, myframe,
			   candFun,sampenv,
			   INTEGER(burnin_R)[0], INTEGER(mcmc_R)[0], 
			   INTEGER(thin_R)[0],
			   INTEGER(verbose)[0], sample_SEXP);
    UNPROTECT(1);
      
    // return the sample
    return sample_SEXP;
  }
}

#endif
