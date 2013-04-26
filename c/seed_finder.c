
#include "memutil.h"
#include "util.h"
#include "seed_finder.h"
#include "seed_table.h"



SeedFinder *seed_finder_new(SeedTable *seed_tab, unsigned int read_len,
			    unsigned int seed_len, unsigned int n_seed) {
  SeedFinder *sf;
  int i, j, min_read_len;

  if(n_seed < 1) {
    my_err("%s:%d: invalid number of seeds (%d)\n", 
	   __FILE__, __LINE__, n_seed);
  }
  
  sf = my_new(SeedFinder, 1);

  sf->seed_tab = seed_tab;
  sf->read_len = read_len;
  sf->seed_len = seed_len;
  sf->n_seed = n_seed;
  sf->match = my_new(SeedMatch, read_len);


  min_read_len = seed_len * n_seed;
  if(read_len < min_read_len) {
    my_err("%s:%d min read len is %d for %d non-overlapping "
	   "seeds of length %u\n", __FILE__, __LINE__, min_read_len, 
	   n_seed, seed_len);
  }

  /* make arrays with dimensions n_seed x read_len */
  sf->lowest_match = my_new(long *, n_seed);
  sf->seed_start = my_new(int *, n_seed);
  for(i = 0; i < n_seed; i++) {
    sf->lowest_match[i] = my_new(long, read_len);
    sf->seed_start[i] = my_new(int, read_len);

    for(j = 0; j < read_len; j++) {
      sf->lowest_match[i][j] = -1;
      sf->seed_start[i][j] = -1;
    }

  }

  sf->best_seeds = my_new(SeedMatch *, n_seed);
  sf->best_read_offsets = my_new0(int, n_seed);

  return sf;
}



void seed_finder_free(SeedFinder *sf) {
  int i;

  for(i = 0; i < sf->n_seed; i++) {
    my_free(sf->lowest_match[i]);
    my_free(sf->seed_start[i]);
  }  
  my_free(sf->lowest_match);
  my_free(sf->seed_start);

  my_free(sf->seed_start);

  my_free(sf->best_seeds);
  my_free(sf->best_read_offsets);

  my_free(sf);
}




void seed_finder_best_seeds(SeedFinder *sf, unsigned char *read_nucs) {
  int i, j, first, left_bound, right_bound, offset;
  long total_match;

  /* lookup number of matches for each possible seed from this read */
  for(i = 0; i < sf->read_len - sf->seed_len + 1; i++) {
    seed_table_lookup(sf->seed_tab, &read_nucs[i], &sf->match[i]);
  }

  /* populate two arrays, with rows corresponding to each seed
   *  lowest_match - fewest total seed matches (for this seed
   *    and seeds to right) if seed i starts at this position or greater
   *  seed_start - seed start position for seed i that yields
   *      fewest matches
   */                        
  for(i = 0; i < sf->n_seed; i++) {

    first = TRUE;

    /* left/right bounds for where seed i can start */
    right_bound = sf->read_len - (sf->seed_len * (i + 1)) - 1;
    left_bound = sf->seed_len * (sf->n_seed - i - 1);
    
    for(j = right_bound; j >= left_bound; j--) {
      if(i == 0) {
	/* this is right-most seed */
	total_match = 0;
      } else {
	/* if we start seed here, what are the fewest total matches 
	 * of the seeds starting to the right of this one?
	 */
	total_match = sf->lowest_match[i-1][j + sf->seed_len];
      }
      
      if(first || (sf->match[j].n_match + total_match) < sf->lowest_match[i][j+1]) { 
	/* starting seed at this position gives fewest total matches */
	sf->lowest_match[i][j] = sf->match[j].n_match + total_match;
	sf->seed_start[i][j] = j;
	first = FALSE;
      } else {
	/* starting seed to right of this position is better */
	sf->lowest_match[i][j] = sf->lowest_match[i][j+1];
	sf->seed_start[i][j] = sf->seed_start[i][j+1];
      }
    }
  }
  
  /** traceback and set best seeds, as well as their offsets into read **/
  offset = sf->seed_start[sf->n_seed-1][0];

  sf->best_read_offsets[0] = offset;
  sf->best_seeds[0] = &sf->match[offset];

  for(i = 1; i < sf->n_seed; i++) {
    offset = sf->seed_start[sf->n_seed-i-1][offset + sf->seed_len];
    sf->best_read_offsets[i] = offset;
    sf->best_seeds[i] = &sf->match[offset];
  }

  /* TODO: could re-order so that seeds with fewest matches are
   * tried first
   */

  /* for debugging: */
  /* fprintf(stderr, "    all seeds:\n    "); */
  /* for(i = 0; i < sf->read_len - sf->seed_len + 1; i++) { */
  /*   fprintf(stderr, " %u", sf->match[i].n_match); */
  /* } */
  /* fprintf(stderr, "\n    lowest matches:\n"); */
  /* for(i = 0; i < sf->n_seed; i++) { */
  /*   fprintf(stderr, "     "); */
  /*   for(j = 0; j < sf->read_len; j++) { */
  /*     fprintf(stderr, " %ld", sf->lowest_match[i][j]); */
  /*   } */
  /*   fprintf(stderr, "\n"); */
  /* } */
  /* fprintf(stderr, "\n    best non-overlapping seeds:\n"); */
  /* for(i = 0; i < sf->n_seed; i++) { */
  /*   fprintf(stderr, "      seed%d: n_match:%u read_offset:%u\n", */
  /* 	    i, sf->best_seeds[i]->n_match, sf->best_read_offsets[i]); */
  /* } */
  
}
