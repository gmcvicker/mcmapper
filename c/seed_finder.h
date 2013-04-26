#ifndef __SEED_FINDER_H__
#define __SEED_FINDER_H__

#include "seed_table.h"


typedef struct {
  SeedTable *seed_tab;

  int read_len;
  int seed_len;
  int n_seed;

  /* seed matches for each possible seed offset in read 
   * read_len
   */
  SeedMatch *match;

  /* arrays (n_seed x read_len) used for dynamic programming algorithm
   * used to find best seeds.
   */
  long **lowest_match;
  int **seed_start;

  /*
   * arrays (length: n_seed) of best non-overlapping seeds
   * and their offsets into the read
   */
  SeedMatch **best_seeds;
  int *best_read_offsets;

} SeedFinder;


void seed_finder_free(SeedFinder *sf);

SeedFinder *seed_finder_new(SeedTable *seed_tab, unsigned int read_len,
			    unsigned int seed_len, unsigned int n_seed);

void seed_finder_best_seeds(SeedFinder *sf, unsigned char *read_nucs);


#endif
