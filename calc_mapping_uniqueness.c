
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>

#include "util.h"
#include "memutil.h"
#include "nuc.h"
#include "chr.h"
#include "chr_table.h"
#include "seed_table.h"
#include "seq.h"



int has_ambi_code(unsigned char *nucs, int len) {
  int i;
  for(i = 0; i < len; i++) {
    if(nucs[i] == NUC_N || nuc_is_ambi(nucs[i])) {
      return TRUE;
    }
  }
  return FALSE;
}


void map_reads(SeedTable *seed_tab, unsigned char *genome_nucs, 
	       long genome_len, long read_len) {
  long i, j, n_match;
  unsigned char *nucs;
  char *seed_str, *read_str;

  nucs = my_new(unsigned char, read_len);

  seed_str = my_new(char, seed_tab->seed_len);
  read_str = my_new(char, read_len);
  
  for(i = 0; i < genome_len - read_len + 1; i++) {
    if(has_ambi_code(&genome_nucs[i], read_len)) {
      /* TODO: fix this so that ambiguity codes can be used */
      continue;
    }

    /* copy nucleotide sequence */
    nucs = memcpy(nucs, &genome_nucs[i], read_len);

    read_str = nuc_ids_to_str(read_str, genome_nucs, read_len);
    fprintf(stderr, "read: %s\n", read_str);

    /* look for seed with lowest number of matches */
    for(j = 0; j < read_len - seed_tab->seed_len + 1; j++) {
      n_match = seed_table_n_match(seed_tab, &nucs[j]);
      nuc_ids_to_str(seed_str, &genome_nucs[j], seed_tab->seed_len);
      fprintf(stderr, "  seed %s has %ld matches\n", seed_str, n_match);
    }

    /* also try reverse strand of reads */
    nuc_ids_revcomp(nucs, read_len);
    for(j = 0; j < read_len - seed_tab->seed_len + 1; j++) {
      n_match = seed_table_n_match(seed_tab, &nucs[j]);
      nuc_ids_to_str(seed_str, &nucs[j], seed_tab->seed_len);
      fprintf(stderr, "  seed %s has %ld matches\n", seed_str, n_match);   
    }

    break;
  }

  my_free(read_str);
  my_free(seed_str);
  my_free(nucs);
}



/**
 * Reads all chromosomes into a single long concatenated array
 */
unsigned char *read_seqs(ChrTable *chr_tab, char **fasta_files, 
			 int n_fasta_files) {  
  int chr_idx, i, n_seq;
  long offset;
  unsigned char *nucs;
  Seq *seq;
  gzFile gzf;
  
  fprintf(stderr, "allocating %ld bytes of memory for genome sequence\n",
	  chr_tab->total_chr_len);

  /* allocate enough memory to hold all chromosomes, set to N */
  nucs = my_new(unsigned char, chr_tab->total_chr_len);
  memset(nucs, NUC_N, chr_tab->total_chr_len);

  seq = seq_new();
  n_seq = 0;
  for(i = 0; i < n_fasta_files; i++) {
    /* read sequence from fasta files, copy into nucleotide array */
    gzf = util_must_gzopen(fasta_files[i], "rb");
    while(seq_read_fasta_record(seq, gzf)) {
      fprintf(stderr, "%s\n", seq->name);
      chr_idx = chr_table_lookup(chr_tab, seq);

      offset = chr_tab->offset[chr_idx];
      memcpy(&nucs[offset], seq->sym, seq->len);
      n_seq += 1;
    }
    gzclose(gzf);
  }

  seq_free(seq);

  if(n_seq < chr_tab->n_chr) {
    my_warn("%s:%d: only read %ld sequences, but %ld are specified in "
	    "chromosome table", __FILE__, __LINE__, n_seq);
  }
  
  return nucs;
}




int main(int argc, char **argv) {
  char **fasta_files, *seed_index_file, *chrom_info_file, *output_dir;
  int n_fasta_files;
  SeedTable *seed_tab;
  ChrTable *chr_tab;
  unsigned char *genome_nucs;

  if(argc < 4) {
    fprintf(stderr, "usage: %s <seed_index_file> <chromInfo.txt> "
	    " <output_dir> <chr1.fa.gz> [<chr2.fa.gz [...]]\n",
	    argv[0]);
    exit(2);
  }
  
  seed_index_file = argv[1];
  chrom_info_file = argv[2];
  output_dir = argv[3];
  fasta_files = &argv[4];
  n_fasta_files = argc - 4;
  
  chr_tab = chr_table_read(chrom_info_file);

  fprintf(stderr, "reading seed index\n");
  seed_tab = seed_table_read(seed_index_file);

  fprintf(stderr, "reading genome sequence\n");
  genome_nucs = read_seqs(chr_tab, fasta_files, n_fasta_files);

  fprintf(stderr, "mapping reads\n");
  map_reads(seed_tab, genome_nucs, chr_tab->total_chr_len, 36);
  
  my_free(genome_nucs);
  seed_table_free(seed_tab);
  chr_table_free(chr_tab);

  return 0;
}
