
#include <stdio.h>
#include <limits.h>
#include <string.h>

#include "nuc.h"
#include "ambi.h"
#include "util.h"
#include "kmer.h"
#include "memutil.h"
#include "seed_table.h"


/**
 * Creates a new SeedTable data structure and returns it
 */
SeedTable *seed_table_new(int seed_len) {
  SeedTable *seed_tab;
  int i;

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

  /* buffer to hold unambiguous nucleotide arrays */
  seed_tab->unambig_nucs = my_new(unsigned char *, 
				  SEED_TABLE_MAX_UNAMBIG);
  for(i = 0; i < SEED_TABLE_MAX_UNAMBIG; i++) {
    seed_tab->unambig_nucs[i] = my_new(unsigned char, 
				       seed_tab->seed_len);
  }
  
  return seed_tab;
}

/**
 * Frees memory allocated for seed table
 */
void seed_table_free(SeedTable *seed_tab) {
  int i;
  
  /* free seed matches */
  my_free(seed_tab->match_buf);
  my_free(seed_tab->match);
  my_free(seed_tab->n_match);
  my_free(seed_tab->cur);
  
  for(i = 0; i < SEED_TABLE_MAX_UNAMBIG; i++) {
    my_free(seed_tab->unambig_nucs[i]);
  }
  my_free(seed_tab->unambig_nucs);

  my_free(seed_tab);
}





/**
 * Counts a seed match, but does not actually add its location
 * to the seed table
 */
void seed_table_count_match(SeedTable *seed_tab, unsigned char *nucs) {
  unsigned int kmer_id;
  unsigned char **unambig;
  int n_unambig, i;

  /* convert seeds with ambiguity codes to all possible 
   * non-ambiguous seqs 
   */
  if(ambi_has_ambi(nucs, seed_tab->seed_len)) {
      n_unambig = ambi_resolve(nucs, seed_tab->seed_len,
			       seed_tab->unambig_nucs, 
			       SEED_TABLE_MAX_UNAMBIG);

      if(n_unambig == 0) {
	my_warn("seed contains too many ambiguous nucleotides");
	return;
      }
      unambig = seed_tab->unambig_nucs;
  } else {
    /* no ambiguous nucleotides, just use original seed */
    unambig = &nucs;
    n_unambig = 1;
  }
  for(i = 0; i < n_unambig; i++) {
    kmer_id = kmer_nucs_to_id(unambig[i], seed_tab->seed_len);
    /* increment number of matches to this kmer */

    if(seed_tab->n_match[kmer_id] == UINT_MAX) {
      my_err("%s:%d maximum number of seed matches (%u) "
	     "exceeded for kmer %u", __FILE__, __LINE__, 
	     UINT_MAX, kmer_id);
    }
    seed_tab->n_match[kmer_id] += 1;
    seed_tab->total_match += 1;
  }
}




/**
 * initializes memory for saving genomic locations of matches
 */
void seed_tab_init_match_mem(SeedTable *seed_tab) {
  unsigned int i;
  long idx;
  char *mem_str;

  /* allocate enough memory to hold all matches */
  mem_str = util_long_to_comma_str(seed_tab->total_match * 
				   sizeof(unsigned int));
  fprintf(stderr, "allocating %s bytes memory to hold seed matches\n", 
	  mem_str);
  my_free(mem_str);

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
  unsigned int kmer_id, i, j;
  int n_unambig;
  unsigned char **unambig;

  if(seed_tab->match_buf == NULL) {
    seed_tab_init_match_mem(seed_tab);
  }

  /* convert seeds with ambiguity codes to all possible 
   * non-ambiguous seqs 
   */
  if(ambi_has_ambi(nucs, seed_tab->seed_len)) {
      n_unambig = ambi_resolve(nucs, seed_tab->seed_len,
			       seed_tab->unambig_nucs, 
			       SEED_TABLE_MAX_UNAMBIG);

      if(n_unambig == 0) {
	my_warn("seed contains too many ambiguous nucleotides");
	return;
      }
      unambig = seed_tab->unambig_nucs;
  } else {
    /* no ambiguous nucleotides, just use original seed */
    unambig = &nucs;
    n_unambig = 1;
  }
  for(i = 0; i < n_unambig; i++) {
    kmer_id = kmer_nucs_to_id(unambig[i], seed_tab->seed_len);

    /* cur is number of matches already added to array */
    j = seed_tab->cur[kmer_id];
    if(j >= seed_tab->n_match[kmer_id]) {
      my_err("%s:%d: more matches than expected to kmer", 
	     __FILE__, __LINE__);
    }

    /* add genomic position (offset) to match array */
    seed_tab->match[kmer_id][j] = offset;

    /* update cur to point to next element of match array */
    seed_tab->cur[kmer_id] += 1;
  }
}


/**
 * sets attributes of the provided SeedMatch data structure to give
 * the total number of matches and kmer id(s) for the provided
 * nucleotide array. The kmer ids can be used to index matches in
 * the seed match table. There can be multiple kmer ids, because the
 * nucleotides are allowed to contain ambiguity codes.
 */
void seed_table_lookup(SeedTable *seed_tab, unsigned char *nucs, 
		       SeedMatch *seed_match) {
  unsigned int kmer_id, i;
  unsigned char **unambig;

  if(ambi_has_ambi(nucs, seed_tab->seed_len)) {
    /* this seed contains ambiguous nucleotides. convert
     * to all possible seeds containing non-ambiguous nucleotides.
     */
    seed_match->n_kmer = ambi_resolve(nucs, seed_tab->seed_len,
				      seed_tab->unambig_nucs,
				      SEED_TABLE_MAX_UNAMBIG);
    
    if(seed_match->n_kmer == 0) {
      my_warn("seed contains too many ambiguous nucleotides");
    }
    unambig = seed_tab->unambig_nucs;

    /* count total number of matches, set kmer ids */
    seed_match->n_match = 0;
    for(i = 0; i < seed_match->n_kmer; i++) {
      kmer_id = kmer_nucs_to_id(unambig[i], seed_tab->seed_len);
      seed_match->n_match += seed_tab->n_match[kmer_id];
      seed_match->kmer_ids[i] = kmer_id;
    }
  } else {
    /* there were no ambiguous nucleotides, just use original seed */
    kmer_id = kmer_nucs_to_id(nucs, seed_tab->seed_len);
    seed_match->n_match = seed_tab->n_match[kmer_id];
    seed_match->n_kmer = 1;
    seed_match->kmer_ids[0] = kmer_id;
  }
}



/**
 * Returns the total number of genomic matches for the provided nucleotide
 * array
 */
unsigned int seed_table_n_match(SeedTable *seed_tab, unsigned char *nucs) {
  SeedMatch match;

  seed_table_lookup(seed_tab, nucs, &match);
  return match.n_match;
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

  fprintf(stderr, "reading seed matches\n");
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
