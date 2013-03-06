
#ifndef __MAPPER_H__
#define __MAPPER_H__

#include "chr_table.h"
#include "seed_table.h"

#define MAP_CODE_NONE 0   /* does not map */
#define MAP_CODE_UNIQUE 1 /* maps uniquely */
#define MAP_CODE_MULTI 2  /* maps multiple times */



typedef struct {
  SeedTable *seed_tab;
  ChrTable *chr_tab;
  unsigned char *genome_nucs;
  long genome_len;
  int allow_mismatch;
} Mapper;


typedef struct {
  unsigned int len;            /* length of seed */
  unsigned int read_idx;       /* index into read */
  SeedMatch match;
} MapSeed;



typedef struct {
  /* fwd nucleotides of read */
  unsigned char *fwd_nucs;

  /* rev nucleotides of read */
  unsigned char *rev_nucs;

  /* length of read */
  unsigned int len;

  /* offset and strand giving genomic location that read maps to */
  unsigned int map_offset;
  char map_strand;
  
  /* does this read map uniquely, multiple times, or not at all? */
  int map_code;

  /* how many mismatches? */
  int n_mismatch;

  /* set to TRUE if read contains an N */
  int has_n;
} MapRead;



void mapper_map_one_read(Mapper *mapper, MapRead *read);

Mapper *mapper_init(SeedTable *seed_tab, ChrTable *chr_tab,
		    char **fasta_files, int n_fasta_files,
		    int max_mismatch);

void mapper_free(Mapper *mapper);

unsigned char *mapper_read_seqs(ChrTable *chr_tab, char **fasta_files, 
				int n_fasta_files);




#endif 
