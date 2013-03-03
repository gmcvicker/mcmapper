#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "seq.h"
#include "seqcoord.h"
#include "chr.h"
#include "util.h"
#include "memutil.h"

#include "chr_table.h"

/**
 * Reads a UCSC chromInfo.txt file, and creates a ChrTable data
 * structure from it. The table contains offsets for both forward and
 * reverse strands of each chromosome so that genomic coordinates can
 * be represented with a single unsigned 32bit integer.
 */
ChrTable *chr_table_read(const char *filename) {
  ChrTable *chr_tab = my_new(ChrTable, 1);
  int i;
  unsigned int len, prev_len;
  
  /* read chromosome names and lengths from file */
  chr_tab->chr_array = chr_read_file(filename, &chr_tab->n_chr);
  
  /* set offsets, total length */
  prev_len = chr_tab->total_chr_len = 0;
  chr_tab->offset = my_new(unsigned int, chr_tab->n_chr);
  for(i = 0; i < chr_tab->n_chr; i++) {
    len = chr_tab->chr_array[i].len;
    chr_tab->offset[i] = chr_tab->total_chr_len;
    chr_tab->total_chr_len += len;

    if(chr_tab->total_chr_len < prev_len) {
      my_err("genome length exceeds max %ld", UINT_MAX);
    }
    prev_len = chr_tab->total_chr_len;
  }
  
  return chr_tab;
}


/**
 * Prints information from ChrTable datastructure to provided file
 *
 */
void chr_table_write(FILE *f, ChrTable *chr_tab) {
  int i;

  for(i = 0; i < chr_tab->n_chr; i++) {
    fprintf(f, "%s %ld %u\n", chr_tab->chr_array[i].name,
	    chr_tab->chr_array[i].len, chr_tab->offset[i]);
  }
}


/**
 * Frees memory allocated for chr_table datastructure.
 */
void chr_table_free(ChrTable *chr_tab) {
  chr_array_free(chr_tab->chr_array, chr_tab->n_chr);
  my_free(chr_tab->offset);
  my_free(chr_tab);
}


/**
 * Returns index of provided sequence in chromosome table, by
 * comparing the name of the sequence against those in the table. Also
 * checks that length of chromosome is correct.
 */
int chr_table_lookup(ChrTable *chr_tab, Seq *seq) {
  int i;

  for(i = 0; i < chr_tab->n_chr; i++) {
    if(strcmp(chr_tab->chr_array[i].name, seq->name) == 0) {
      /* found entry with matching name */
      if(chr_tab->chr_array[i].len != seq->len) {
	my_err("%s:%d: length of %s chromosome table (%d) "
	       "does not match length of sequence (%d)\n",
	       __FILE__, __LINE__, seq->name, chr_tab->chr_array[i].len, 
	       seq->len);
      }

      return i;
    }
  }

  my_err("%s:%d: could not find sequence %s in chromosome table",
	 __FILE__, __LINE__, seq->name);

  return -1;
}


/**
 * Converts an integer offset to a sequence coordinate. Populates
 * provided SeqCoord data structure (does not create a new one),
 * with start, end and chromosome. Returns index of chromosome
 * in chr_tab->chr_array.
 */
int chr_table_offset_to_coord(ChrTable *chr_tab, unsigned int offset, 
			      SeqCoord *c) {
  int i, found_offset, chr_idx;
  
  found_offset = FALSE;

  for(i = 0; i < chr_tab->n_chr; i++) {
    if(offset >= chr_tab->offset[i]) {
      c->start = offset - chr_tab->offset[i] + 1;
      c->end = c->start;
      c->strand = STRAND_FWD;
      c->chr = &chr_tab->chr_array[i];
      c->seqname = NULL;
      found_offset = TRUE;
      chr_idx = i;
    }
  }

  if(!found_offset) {
    my_err("%s:%d: invalid offset %d", __FILE__, __LINE__, offset);
  }

  if(c->end > c->chr->len) {
    my_err("%s:%d: invalid offset %d implies position past end of chromosome",
	   __FILE__, __LINE__, offset);
  }

  return chr_idx;
}

