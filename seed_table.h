#ifndef __SEED_TABLE_H__
#define __SEED_TABLE_H__

#include <zlib.h>

/* magic number, used to verify byte order and formatting of binary seed
 * table is correct.
 */
#define SEED_TABLE_MAGIC 457401295


typedef struct {
  int seed_len;        /* seed length */
  unsigned int n_seed; /* number of seeds */



  unsigned int total_match;  /* total number of matches in table so far */

  unsigned int *n_match; /* number of genomic matches for each seed */
  unsigned int **match;  /* array of genomic coords (offsets) seed  */
  unsigned int *cur;

  /* memory buffer for storing seed matches */
  unsigned int *match_buf;

} SeedTable;



void seed_table_count_match(SeedTable *seed_tab, unsigned char *nucs);

void seed_table_add_match(SeedTable *seed_tab, unsigned int offset,
			  unsigned char *nucs);

unsigned int seed_table_get_matches(SeedTable *seed_tab, 
				    unsigned char *nucs,
				    unsigned int **match_offsets);

void seed_table_write(SeedTable *seed_table, gzFile f);

SeedTable *seed_table_new(int seed_len);
SeedTable *seed_table_read(const char *filename);
void seed_table_free(SeedTable *seed_table);


#endif
