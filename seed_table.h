#ifndef __SEED_TABLE_H__
#define __SEED_TABLE_H__

#include <zlib.h>

/* magic number, used to verify byte order and formatting of binary seed
 * table is correct.
 */
#define SEED_TABLE_MAGIC 457401295

typedef struct SeedMatch_t SeedMatch;

struct SeedMatch_t {
  /* offset that can be converted to genomic coordinate */
  unsigned int offset; 
  SeedMatch *prev; /* ptr to prev match or NULL if first match */
};


typedef struct {
  int seed_len;        /* seed length */
  unsigned int n_seed; /* number of seeds */

  unsigned int *n_seed_match; /* number of matches for each seed */
  SeedMatch **seed_match;     /* linked lists of matches for each seed */

  size_t total_match;       /* total number of matches in table so far */
  SeedMatch *seed_match_buf; /* buffer of allocated memory for seed matches */
  size_t seed_match_buf_sz; /* size of seed match buffer */
} SeedTable;



void seed_table_add_match(SeedTable *seed_tab, unsigned int offset,
			  unsigned char *nucs);

void seed_table_write(SeedTable *seed_table, gzFile f);

SeedTable *seed_table_new(int seed_len, size_t buf_sz);
SeedTable *seed_table_read(const char *filename);
void seed_table_free(SeedTable *seed_table);

unsigned int seed_table_n_match(SeedTable *seed_tab, unsigned char *nucs);


#endif
