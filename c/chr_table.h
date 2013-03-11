
#ifndef __CHR_TABLE_H__
#define __CHR_TABLE_H__

#include "seq.h"
#include "seqcoord.h"
#include "chr.h"

typedef struct {
  Chromosome *chr_array;
  int n_chr;
  unsigned int total_chr_len;

  /* These are offsets for the start of each chromosome. We pretend
   * that the entire genome sequence is concatenated together in one
   * large array. This makes it easy to represent a genomic position
   * with a single long integer.
   */
  unsigned int *offset;
} ChrTable;


ChrTable *chr_table_read(const char *filename);
void chr_table_write(FILE *f, ChrTable *chr_tab);
void chr_table_free(ChrTable *chr_tab);

int chr_table_lookup(ChrTable *chr_tab, Seq *seq);

int chr_table_offset_to_coord(ChrTable *chr_tab, unsigned int offset, 
			      SeqCoord *c);

#endif
