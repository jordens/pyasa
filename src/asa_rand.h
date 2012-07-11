/*
 *   pyasa - python bindings for Adaptive Simulated Annealing
 */

#if ASA_LIB
  LONG_INT asa_seed (LONG_INT seed);
#endif
  double myrand (LONG_INT * rand_seed);
  double randflt (LONG_INT * rand_seed);
  double resettable_randflt (LONG_INT * rand_seed, int reset);

#define SHUFFLE 256             /* size of random array */
