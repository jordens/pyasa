/*
 *   pyasa - python bindings for Adaptive Simulated Annealing
 */

#include "asa_usr_asa.h"
#include "asa_rand.h"

#if ASA_LIB
static LONG_INT *asa_rand_seed;
#endif


/* Here is a good random number generator */

#define MULT ((LONG_INT) 25173)
#define MOD ((LONG_INT) 65536)
#define INCR ((LONG_INT) 13849)
#define FMOD ((double) 65536.0)

#if ASA_LIB
/***********************************************************************
* LONG_INT asa_seed - returns initial random seed
***********************************************************************/

#if HAVE_ANSI
LONG_INT
asa_seed (LONG_INT seed)
#else
LONG_INT
asa_seed (seed)
     LONG_INT seed;
#endif
{
  static LONG_INT rand_seed;

  if (fabs ((double) seed) > 0) {
    asa_rand_seed = &rand_seed;
    rand_seed = seed;
  }

  return (rand_seed);
}
#endif /* ASA_LIB */

/***********************************************************************
* double myrand - returns random number between 0 and 1
*	This routine returns the random number generator between 0 and 1
***********************************************************************/

#if HAVE_ANSI
double
myrand (LONG_INT * rand_seed)
#else
double
myrand (rand_seed)
     LONG_INT *rand_seed;
#endif
  /* returns random number in {0,1} */
{
#if TRUE                        /* (change to FALSE for alternative RNG) */
  *rand_seed = (LONG_INT) ((MULT * (*rand_seed) + INCR) % MOD);
  return ((double) (*rand_seed) / FMOD);
#else
  /* See "Random Number Generators: Good Ones Are Hard To Find,"
     Park & Miller, CACM 31 (10) (October 1988) pp. 1192-1201.
     ***********************************************************
     THIS IMPLEMENTATION REQUIRES AT LEAST 32 BIT INTEGERS
     *********************************************************** */
#define _A_MULTIPLIER  16807L
#define _M_MODULUS     2147483647L      /* (2**31)-1 */
#define _Q_QUOTIENT    127773L  /* 2147483647 / 16807 */
#define _R_REMAINDER   2836L    /* 2147483647 % 16807 */
  long lo;
  long hi;
  long test;

  hi = *rand_seed / _Q_QUOTIENT;
  lo = *rand_seed % _Q_QUOTIENT;
  test = _A_MULTIPLIER * lo - _R_REMAINDER * hi;
  if (test > 0) {
    *rand_seed = test;
  } else {
    *rand_seed = test + _M_MODULUS;
  }
  return ((double) *rand_seed / _M_MODULUS);
#endif /* alternative RNG */
}

/***********************************************************************
* double randflt
***********************************************************************/

#if HAVE_ANSI
double
randflt (LONG_INT * rand_seed)
#else
double
randflt (rand_seed)
     LONG_INT *rand_seed;
#endif
{
  return (resettable_randflt (rand_seed, 0));
}

/***********************************************************************
* double resettable_randflt
***********************************************************************/

#if HAVE_ANSI
double
resettable_randflt (LONG_INT * rand_seed, int reset)
#else
double
resettable_randflt (rand_seed, reset)
     LONG_INT *rand_seed;
     int reset;
#endif
  /* shuffles random numbers in random_array[SHUFFLE] array */
{

  /* This RNG is a modified algorithm of that presented in
   * %A K. Binder
   * %A D. Stauffer
   * %T A simple introduction to Monte Carlo simulations and some
   *    specialized topics
   * %B Applications of the Monte Carlo Method in statistical physics
   * %E K. Binder
   * %I Springer-Verlag
   * %C Berlin
   * %D 1985
   * %P 1-36
   * where it is stated that such algorithms have been found to be
   * quite satisfactory in many statistical physics applications. */

  double rranf;
  unsigned kranf;
  int n;
  static int initial_flag = 0;
  LONG_INT initial_seed;
#if ASA_SAVE
  /* random_array[] local to all of asa_usr.c set at top of file */
#else
  static double random_array[SHUFFLE];  /* random variables */
#endif

  if (*rand_seed < 0)
    *rand_seed = -*rand_seed;

  if ((initial_flag == 0) || reset) {
    initial_seed = *rand_seed;

    for (n = 0; n < SHUFFLE; ++n)
      random_array[n] = myrand (&initial_seed);

    initial_flag = 1;

    for (n = 0; n < 1000; ++n)  /* warm up random generator */
      rranf = randflt (&initial_seed);

    rranf = randflt (rand_seed);

    return (rranf);
  }

  kranf = (unsigned) (myrand (rand_seed) * SHUFFLE) % SHUFFLE;
  rranf = *(random_array + kranf);
  *(random_array + kranf) = myrand (rand_seed);

  return (rranf);
}
