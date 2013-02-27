/* TODO: 
 * - Allow for polymorphisms (ambiguity codes).
 * - Reduce memory footprint. Only allow up to MAX_SEED_MATCH matches?
 */

#include <stdio.h>
#include <stdlib.h>

#include "memutil.h"
#include "util.h"
#include "err.h"
#include "seqcoord.h"
#include "seq.h"
#include "nuc.h"

#include "ambi.h"
#include "chr_table.h"
#include "seed_table.h"


int has_n(unsigned char *nucs, int len) {
  int i;
  for(i = 0; i < len; i++) {
    if(nucs[i] == NUC_N) {
      return TRUE;
    }
  }
  return FALSE;
}



void count_matches(ChrTable *chr_tab, SeedTable *seed_tab, Seq *seq) {
  unsigned int i;
  
  /* just call this to verify that chromosome is in table */
  chr_table_lookup(chr_tab, seq);
  
  /* loop over fwd strand of chromosome */
  fprintf(stderr, "counting seeds\n");
  for(i = 0; i < seq->len - seed_tab->seed_len + 1; i++) {
    if(has_n(&seq->sym[i], seed_tab->seed_len)) {
      continue;
    }

    if((i % 1000000) == 0) {
      fprintf(stderr, ".");
    }
    /* count occurance of this seed */
    seed_table_count_match(seed_tab, &seq->sym[i]);
  }
  fprintf(stderr, "\n");
}




void add_matches(ChrTable *chr_tab, SeedTable *seed_tab, Seq *seq) {
  int chr_tab_idx;
  unsigned int i, offset;
  
  chr_tab_idx = chr_table_lookup(chr_tab, seq);

  offset = chr_tab->offset[chr_tab_idx];
  
  /* loop over fwd strand of chromosome */
  fprintf(stderr, "loading seeds from fwd strand\n");
  for(i = 0; i < seq->len - seed_tab->seed_len + 1; i++) {
    if(has_n(&seq->sym[i], seed_tab->seed_len)) {
      continue;
    }
    if((i % 1000000) == 0) {
      fprintf(stderr, ".");
    }
    seed_table_add_match(seed_tab, offset + i, &seq->sym[i]);
  }
  fprintf(stderr, "\n");
}



int main(int argc, char **argv) {
  char **fasta_files;
  int seed_len, i, n_fasta_files;
  ChrTable *chr_tab;
  SeedTable *seed_tab;
  Seq *seq;
  gzFile gzf, out_gzf;
  char *out_filename;
  
  if(argc < 4) {
    fprintf(stderr, "usage: %s <seed_len> <chromInfo.txt> "
	    "<output_seed_index.gz> [chr1.fa.gz [chr2.fa.gz [...]]]\n", 
	    argv[0]);
    exit(2);
  }
  
  
  seed_len = util_parse_long(argv[1]);

  fasta_files = &argv[4];
  n_fasta_files = argc - 4;
  out_filename = argv[3];

  /* read chromosomes and make table containing offsets for both
   * forward and reverse strands (so we can represent genomic
   * coordinates with a single long integer).
   */
  chr_tab = chr_table_read(argv[2]);
  fprintf(stderr, "there are %d chromosomes, total length: %u\n",
	  chr_tab->n_chr, chr_tab->total_chr_len);

  /* create a table to hold seed matches */
  fprintf(stderr, "initializing seed table\n");
  seed_tab = seed_table_new(seed_len);

  /* open an output file to write seed table to */
  if(util_file_exists(out_filename)) {
    my_err("output file %s already exists\n", out_filename);
    exit(2);
  }
  out_gzf = util_must_gzopen(out_filename, "wb");

  /*
   * first pass: count number of matches to each seed 
   */
  seq = seq_new();
  for(i = 0; i < n_fasta_files; i++) {
    fprintf(stderr, "reading sequence from file %s\n", fasta_files[i]);
    gzf = util_must_gzopen(fasta_files[i], "rb");
    while(seq_read_fasta_record(seq, gzf)) {
      fprintf(stderr, "%s %ld\n", seq->name, seq->len);
      fprintf(stderr, "counting seed matches\n");
      count_matches(chr_tab, seed_tab, seq);
    }
    gzclose(gzf);
  }

  /*
   * second pass: store location of each match
   */
  for(i = 0; i < n_fasta_files; i++) {
    fprintf(stderr, "reading sequence from file %s\n", fasta_files[i]);
    gzf = util_must_gzopen(fasta_files[i], "rb");
    while(seq_read_fasta_record(seq, gzf)) {
      fprintf(stderr, "%s %ld\n", seq->name, seq->len);
      fprintf(stderr, "recording seed match positions\n");
      add_matches(chr_tab, seed_tab, seq);
    }
    gzclose(gzf);
  }

  seq_free(seq);

  /* write table to file in binary format */
  fprintf(stderr, "writing seed table to file %s\n", out_filename);
  seed_table_write(seed_tab, out_gzf);
  gzclose(out_gzf);

  chr_table_free(chr_tab);
  seed_table_free(seed_tab);
  
  return 0;
}

