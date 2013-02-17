
#include <stdio.h>

#include "util.h"
#include "kmer.h"
#include "memutil.h"
#include "seed_table.h"

/**
 * Creates a new SeedTable data structure and returns it
 */
SeedTable *seed_table_new(int seed_len, size_t buf_sz) {
  SeedTable *seed_tab;

  seed_tab = my_new(SeedTable, 1);
  seed_tab->seed_len = seed_len;
  seed_tab->n_seed = kmer_num_kmers(seed_len);

  seed_tab->n_seed_match = my_new0(unsigned int, seed_tab->n_seed);

  /* linked lists of seed matches, one for each possible seed */
  seed_tab->seed_match = my_new0(SeedMatch *, seed_tab->n_seed);

  /* allocate memory for seed matches */
  fprintf(stderr, "allocating buffer for seed matches: %lu x %lu"
	  " = %lu bytes\n", sizeof(SeedMatch), buf_sz,  
	  buf_sz*sizeof(SeedMatch));
  seed_tab->total_match = 0;
  seed_tab->seed_match_buf = my_new(SeedMatch, buf_sz);
  seed_tab->seed_match_buf_sz = buf_sz;
  
  return seed_tab;
}


void seed_table_free(SeedTable *seed_tab) {
  unsigned int i;
  SeedMatch *match, *prev_match;

  /* free seed matches */
  for(i = 0; i < seed_tab->n_seed; i++) {
    match = seed_tab->seed_match[i];

    while(match) {
      prev_match = match->prev;
      my_free(match);
      match = prev_match;
    }
  }

  my_free(seed_tab->n_seed_match);
  my_free(seed_tab->seed_match);

  my_free(seed_tab->seed_match_buf);

  my_free(seed_tab);
}




/**
 * Converts provided array of nucleotides to a kmer_id, and adds a 
 * seed match to the provided seed table.
 */
void seed_table_add_match(SeedTable *seed_tab, unsigned int offset,
			  unsigned char *nucs) {
  unsigned int kmer_id;
  SeedMatch *match;

  kmer_id = kmer_nucs_to_id(nucs, seed_tab->seed_len);

  /* create and init new match data structure,
   * using memory already allocated in seed match buffer
   */
  if(seed_tab->total_match >= seed_tab->seed_match_buf_sz) {
    my_err("out of buffer space for seed matches\n");
  }
  match = &seed_tab->seed_match_buf[seed_tab->total_match];
  seed_tab->total_match += 1;
  match->offset = offset;

  /* TODO: may not be necessary to save identical adjacent matches
   * e.g. if seed is AAAAAAAAA, may be no need to record that it
   * matches to multiple adjacent positions.
   */

  /* maintain linked list of matches to same seed */
  match->prev = seed_tab->seed_match[kmer_id];
  seed_tab->seed_match[kmer_id] = match;
  seed_tab->n_seed_match[kmer_id] += 1;
}


/**
 * Returns the number of matches to the provided seed
 */
unsigned int seed_table_n_match(SeedTable *seed_tab, unsigned char *nucs) {
  unsigned int kmer_id;
  kmer_id = kmer_nucs_to_id(nucs, seed_tab->seed_len);
  return seed_tab->n_seed_match[kmer_id];
}


/**
 * writes a seed table to a gzipped file in binary format
 */
void seed_table_write(SeedTable *seed_table, gzFile f) {
  unsigned int i, j, magic;
  SeedMatch *match;
  
  magic = SEED_TABLE_MAGIC;
  util_gzwrite_one(f, magic);
  util_gzwrite_one(f, seed_table->seed_len);

  for(i = 0; i < seed_table->n_seed; i++) {
    /* write number of matches to this seed */
    util_gzwrite_one(f, seed_table->n_seed_match[i]);

    /* write offsets for each seed match */
    match = seed_table->seed_match[i];
    for(j = 0; j < seed_table->n_seed_match[i]; j++) {
      if(match == NULL) {
	my_err("%s:%d: expected %d seed matches, but only got %d",
	       __FILE__, __LINE__, seed_table->n_seed_match[i], j);
      }

      util_gzwrite_one(f, match->offset);
      
      match = match->prev;
    }
  }
}


/**
 * Creates a new SeedTable by reading a binary file
 */
SeedTable *seed_table_read(const char *filename) {
  unsigned int i, j;
  unsigned int magic;
  int seed_len;
  unsigned int total_match;
  SeedTable *seed_tab;
  SeedMatch *match;
  gzFile f;
  
  f = util_must_gzopen(filename, "rb");

  util_gzread_one(f, magic);
  
  if(magic != SEED_TABLE_MAGIC) {
    my_err("%s:%d: seed table magic number does not match. "
	   "incorrect byte order?", __FILE__, __LINE__);
  }

  util_gzread_one(f, seed_len);
  util_gzread_one(f, total_match);

  seed_tab = seed_table_new(seed_len, total_match);

  for(i = 0; i < seed_tab->n_seed; i++) {
    /* read how many matches to this seed there are */
    util_gzread_one(f, seed_tab->n_seed_match[i]);
    
    for(j = 0; j < seed_tab->n_seed_match[i]; j++) {
      /* init new match data structure, using memory already allocated
       * in seed match buffer
       */
      if(seed_tab->total_match >= seed_tab->seed_match_buf_sz) {
	my_err("out of buffer space for seed matches\n");
      }
      match = &seed_tab->seed_match_buf[seed_tab->total_match];
      seed_tab->total_match += 1;

      util_gzread_one(f, match->offset);
      match->prev = seed_tab->seed_match[i];
      seed_tab->seed_match[i] = match;
    }
  }

  gzclose(f);

  return seed_tab;
}
