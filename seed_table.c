
#include <stdio.h>
#include <limits.h>

#include "util.h"
#include "kmer.h"
#include "memutil.h"
#include "seed_table.h"

/**
 * Creates a new SeedTable data structure and returns it
 */
SeedTable *seed_table_new(int seed_len) {
  SeedTable *seed_tab;

  seed_tab = my_new(SeedTable, 1);
  seed_tab->seed_len = seed_len;
  seed_tab->n_seed = kmer_num_kmers(seed_len);

  /* number of matches for each seed (init to 0) */
  seed_tab->n_match = my_new0(unsigned int, seed_tab->n_seed);

  /* arrays of seed matches (init to NULL) */
  seed_tab->match = my_new0(unsigned int *, seed_tab->n_seed);

  /* number of matches that have actually been added */
  seed_tab->cur = my_new0(unsigned int, seed_tab->n_seed);

  seed_tab->total_match = 0;
  seed_tab->match_buf = NULL;
  
  return seed_tab;
}


void seed_table_free(SeedTable *seed_tab) {
  /* free seed matches */
  my_free(seed_tab->match_buf);
  my_free(seed_tab->match);
  my_free(seed_tab->n_match);
  my_free(seed_tab->cur);

  my_free(seed_tab);
}



/**
 * Counts a seed match, but does not actually add its location
 * to the seed table
 */
void seed_table_count_match(SeedTable *seed_tab, unsigned char *nucs) {
  unsigned int kmer_id;
  kmer_id = kmer_nucs_to_id(nucs, seed_tab->seed_len);
  /* increment number of matches to this kmer */

  if(seed_tab->n_match[kmer_id] == UINT_MAX) {
    my_err("%s:%d maximum number of seed matches (%u) exceeded for kmer %u",
	   __FILE__, __LINE__, UINT_MAX, kmer_id);
  }
  seed_tab->n_match[kmer_id] += 1;
  seed_tab->total_match += 1;
}


/**
 * initializes memory for saving genomic locations of matches
 */
void seed_tab_init_match_mem(SeedTable *seed_tab) {
  unsigned int i;
  long idx;


  /* allocate enough memory to hold all matches */
  seed_tab->match_buf = my_new(unsigned int, seed_tab->total_match);

  /* each match array points to a different location in the buffer */
  idx = 0;
  for(i = 0; i < seed_tab->n_seed; i++) {
    if(seed_tab->n_match[i]) {
      seed_tab->match[i] = &seed_tab->match_buf[idx];
      /* advance ptr by number of matches */
      idx += seed_tab->n_match[i];
    } else {
      seed_tab->match[i] = NULL;
    }
  }
}


/**
 * Adds location for a seed match to the provided seed table.
 */
void seed_table_add_match(SeedTable *seed_tab, unsigned int offset,
			  unsigned char *nucs) {
  unsigned int kmer_id, i;

  if(seed_tab->match_buf == NULL) {
    seed_tab_init_match_mem(seed_tab);
  }

  kmer_id = kmer_nucs_to_id(nucs, seed_tab->seed_len);
  /* cur is number of matches already added to array */
  i = seed_tab->cur[kmer_id];
  if(i >= seed_tab->n_match[kmer_id]) {
    my_err("%s:%d: more matches than expected to kmer", __FILE__, __LINE__);
  }

  /* add genomic position (offset) to match array */
  seed_tab->match[kmer_id][i] = offset;

  /* update cur to point to next element of match array */
  seed_tab->cur[kmer_id] += 1;
}


/**
 * Returns the number of matches to the provided seed and sets provided
 * pointer to the array of genomic offsets (if the ptr is not NULL).
 */
unsigned int seed_table_get_matches(SeedTable *seed_tab, 
				    unsigned char *nucs, 
				    unsigned int **match_offsets) {
  unsigned int kmer_id;
  kmer_id = kmer_nucs_to_id(nucs, seed_tab->seed_len);

  if(match_offsets) {
    *match_offsets = seed_tab->match[kmer_id];
  }
  
  return seed_tab->n_match[kmer_id];
}




/**
 * writes a seed table to a gzipped file in binary format
 */
void seed_table_write(SeedTable *seed_table, gzFile f) {
  unsigned int i, j, magic;
  
  magic = SEED_TABLE_MAGIC;
  util_gzwrite_one(f, magic);
  util_gzwrite_one(f, seed_table->seed_len);
  util_gzwrite_one(f, seed_table->total_match);

  /* write number of matches to each seed */
  for(i = 0; i < seed_table->n_seed; i++) {
    util_gzwrite_one(f, seed_table->n_match[i]);
  }

  /* write offsets for each seed match */
  for(i = 0; i < seed_table->n_seed; i++) {
    for(j = 0; j < seed_table->n_match[i]; j++) {
      util_gzwrite_one(f, seed_table->match[i][j]);
    }
  }
}


/**
 * Creates a new SeedTable by reading a binary file
 */
SeedTable *seed_table_read(const char *filename) {
  unsigned int i;
  unsigned int magic;
  int seed_len;
  unsigned int total_match;
  SeedTable *seed_tab;
  gzFile f;
  
  f = util_must_gzopen(filename, "rb");

  util_gzread_one(f, magic);
  
  if(magic != SEED_TABLE_MAGIC) {
    my_err("%s:%d: seed table magic number does not match. "
	   "incorrect byte order?", __FILE__, __LINE__);
  }

  util_gzread_one(f, seed_len);
  util_gzread_one(f, total_match);

  fprintf(stderr, "seed_len: %u, total_match: %u\n",
	  seed_len, total_match);
  seed_tab = seed_table_new(seed_len);
  seed_tab->total_match = total_match;

  fprintf(stderr, "reading seed match table\n");
  for(i = 0; i < seed_tab->n_seed; i++) {
    /* read how many matches to this seed there are */
    util_gzread_one(f, seed_tab->n_match[i]);
  }

  /* allocate memory to hold matches */
  seed_tab_init_match_mem(seed_tab);

  /* read matches */
  for(i = 0; i < seed_tab->n_seed; i++) {
    if(seed_tab->n_match[i] > 0) {
      util_must_gzread(f, seed_tab->match[i], 
		       seed_tab->n_match[i] * sizeof(unsigned int));
    }

    if((i % 1000000) == 0) {
      fprintf(stderr, ".");
    }
  }

  fprintf(stderr, "\n");

  gzclose(f);

  return seed_tab;
}
