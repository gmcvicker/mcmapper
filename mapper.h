
#ifndef __MAPPER_H__
#define __MAPPER_H__

#include "chr_table.h"
#include "seed_table.h"

#define MAP_CODE_NONE 0   /* does not map */
#define MAP_CODE_UNIQUE 1 /* maps uniquely */
#define MAP_CODE_MULTI 2  /* maps multiple times */


typedef struct {
  unsigned int len;
  unsigned int read_idx;
  unsigned int n_match;
  unsigned int *matches;
} MapSeed;


typedef struct {
  /* fwd nucleotides of read */
  unsigned char *fwd_nucs;

  /* rev nucleotides of read */
  unsigned char *rev_nucs;

  /* length of read */
  unsigned int len;

  MapSeed fwd_seed;
  MapSeed rev_seed;

  /* offset and strand giving genomic location that read maps to */
  unsigned int map_offset;
  char map_strand;
  
  /* does this read map uniquely, multiple times, or not at all? */
  int map_code;
} MapRead;



void mapper_map_one_read(SeedTable *seed_tab, unsigned char *genome_nucs,
			 long genome_len, MapRead *read);

void mapper_write_read(ChrTable *chr_tab, MapRead *read);


unsigned char *mapper_read_seqs(ChrTable *chr_tab, char **fasta_files, 
				int n_fasta_files);




#endif 
