#ifndef __SEED_TABLE_H__
#define __SEED_TABLE_H__

#include <zlib.h>

/* magic number, used to verify byte order and formatting of binary seed
 * table is correct.
 */
#define SEED_TABLE_MAGIC 457401295
#define SEED_TABLE_MAX_UNAMBIG 128


typedef struct {
  int seed_len;        /* seed length */
  unsigned int n_seed; /* number of seeds */

  unsigned int total_match;  /* total number of matches in table so far */

  unsigned int *n_match; /* number of genomic matches for each seed */
  unsigned int **match;  /* array of genomic coords (offsets) seed  */
  unsigned int *cur;

  /* memory buffer for storing seed matches */
  unsigned int *match_buf;

  /* memory buffer for storing unambiguous seeds (parsed
   * from seeds containing ambiguity codes)
   */
  unsigned char **unambig_nucs;

} SeedTable;



typedef struct {
  unsigned int n_match;  /* total number of genomic matches */
  unsigned int n_kmer;  /* number of unambiguous kmers from seed */

  /* ids of kmers, used to lookup genomic locations of matches in seed table */
  unsigned int kmer_ids[SEED_TABLE_MAX_UNAMBIG];
} SeedMatch;



void seed_table_count_match(SeedTable *seed_tab, unsigned char *nucs);

void seed_table_add_match(SeedTable *seed_tab, unsigned int offset,
			  unsigned char *nucs);

void seed_table_lookup(SeedTable *seed_tab, unsigned char *nucs,
		       SeedMatch *seed_match);

unsigned int seed_table_n_match(SeedTable *seed_tab, unsigned char *nucs);

void seed_table_write(SeedTable *seed_table, gzFile f);

SeedTable *seed_table_new(int seed_len);
SeedTable *seed_table_read(const char *filename);
void seed_table_free(SeedTable *seed_table);


#endif
