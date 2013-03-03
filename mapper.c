
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>

#include "ambi.h"
#include "util.h"
#include "memutil.h"
#include "nuc.h"
#include "chr.h"
#include "chr_table.h"
#include "seed_table.h"
#include "seq.h"
#include "mapper.h"


/**
 * Given a read, identifies the seed with the fewest number of
 * matches, genome-wide. Returns the number of matches and sets the
 * index into the read and ptr to match offsets
 */
void get_best_seed(SeedTable *seed_tab, unsigned char *read_nucs,
		   unsigned int read_len, MapSeed *seed) {

  unsigned int i, n_match, lowest_n_match, lowest_idx;

  if(seed_tab->seed_len > read_len) {
    my_err("%s:%d: seed len (%d) > than read len (%d)", __FILE__, __LINE__,
	   seed_tab->seed_len, read_len);
  }

  /* start with seed at beginning of read */
  lowest_n_match = seed_table_get_matches(seed_tab, read_nucs, NULL);
  lowest_idx = 0;

  /* look for seed with lowest number of matches */
  for(i = 1; i < read_len - seed_tab->seed_len + 1; i++) {
    n_match = seed_table_get_matches(seed_tab, &read_nucs[i], NULL);

    if(n_match < lowest_n_match) {
      lowest_n_match = n_match;
      lowest_idx = i;
    }
  }

  /* retrieve genomic match locations for seed with lowest matches */
  seed->read_idx = lowest_idx;
  seed->n_match = seed_table_get_matches(seed_tab, &read_nucs[lowest_idx], 
					 &seed->matches);
}


void write_match(FILE *f, unsigned char *genome_nucs, MapRead *read, 
		 MapSeed *seed, unsigned int match_idx, char strand) {

  char genome_str[1024];
  char seed_str[1024];
  char fwd_read_str[1024];
  char rev_read_str[1024];
  long seed_offset,  read_genome_start, i;

  seed_offset = (long)seed->matches[match_idx];

  fprintf(stderr, "match_idx: %u, matches[%u]: %u, "
  	  "seed_offset: %ld, strand=%d\n",
  	  match_idx, match_idx, seed->matches[match_idx],
  	  seed_offset, strand);

  if(strand == STRAND_FWD) {
    nuc_ids_to_str(seed_str, &read->fwd_nucs[seed->read_idx], seed->len);
  } else {
    nuc_ids_to_str(seed_str, &read->rev_nucs[seed->read_idx], seed->len);
  }
  nuc_ids_to_str(fwd_read_str, read->fwd_nucs, read->len);
  nuc_ids_to_str(rev_read_str, read->rev_nucs, read->len);
  
  read_genome_start = seed_offset - seed->read_idx;
  nuc_ids_to_str(genome_str, &genome_nucs[read_genome_start], read->len);

  fprintf(f, "seed_offset: %ld, read_genome_start: %ld\n",
	  seed_offset, read_genome_start);
  fprintf(f, "fwd_read_str: %s\n", fwd_read_str);
  fprintf(f, "seed_str    : ");
  for(i = 0; i < seed->read_idx; i++) {
    fprintf(f, " ");
  }
  fprintf(f, "%s\n", seed_str);
  fprintf(f, "genome_str  : %s\n", genome_str);
  /* fprintf(f, "rev_read_str: %s\n", rev_read_str); */

}



void align_read(unsigned char *genome_nucs, long genome_len, 
		MapRead *read, char strand) {  
  long i, j, genome_seed_start, genome_seed_end;
  long genome_read_start, genome_read_end;
  unsigned char *read_nucs;
  MapSeed *seed;
  int perfect_match;


  if(strand == STRAND_FWD) {
    read_nucs = read->fwd_nucs;
    seed = &read->fwd_seed;
  }
  else if(strand == STRAND_REV) {
    read_nucs = read->rev_nucs;
    seed = &read->rev_seed;
  } else {
    read_nucs = NULL;
    seed = NULL;
    my_err("%s:%d: unknown strand", __FILE__, __LINE__);
  }

  /* try to align starting at each seed */
  /* currently we allow 0 mismatches */
  for(i = 0; i < seed->n_match; i++) {
    perfect_match = TRUE;

    /* get genomic offsets for start / end of read */
    genome_read_start = (long)seed->matches[i] - (long)seed->read_idx;
    genome_read_end   = genome_read_start + (long)read->len - 1;

    /* genomic offsets of start / end of seed match */
    genome_seed_start = (long)seed->matches[i];
    genome_seed_end   = (long)seed->matches[i] + (long)seed->len - 1;

    if((genome_read_start < 0) || (genome_read_end >= genome_len)) {
      /* read would overhang end of genome */
      continue;
    }
    
    /* check that left end of read (before seed) matches genome 
     * allowing for ambiguity codes
     */
    j = 0;
    while(perfect_match && (j < seed->read_idx)) {
      if(!ambi_nucs_match(genome_nucs[genome_read_start + j], read_nucs[j])) {
	perfect_match = FALSE;
      }
      j++;
    }

    /* check that right end of read (after seed) matches genome 
     * allowing for ambiguity codes
     */
    j = seed->read_idx + seed->len;
    while(perfect_match && (j < read->len)) {
      if(!ambi_nucs_match(genome_nucs[genome_read_start + j], read_nucs[j])) {
	perfect_match = FALSE;
      } 
      j++;
    }

    /* if(perfect_match) { */
    /*   fprintf(stderr, "MATCH!\n"); */
    /* } else { */
    /*   fprintf(stderr, "MISMATCH\n"); */
    /* } */
    /* write_match(stderr, genome_nucs, read, seed, i, strand); */


    if(perfect_match) {
      /* read matches! */
      if(read->map_code == MAP_CODE_NONE) {
	/* first match we have observed */
	read->map_code = MAP_CODE_UNIQUE;
	read->map_offset = genome_read_start;
	read->map_strand = strand;
      }
      else if(read->map_code == MAP_CODE_UNIQUE) {
	/* second match we have observed */
	read->map_code = MAP_CODE_MULTI;
	/* fprintf(stderr, "  read matches multiple times\n"); */
	break;
      }
      else {
	my_err("%s:%d: expected mapping code to be NONE or UNIQUE",
	       __FILE__, __LINE__);
      }
    }
  }
}



void mapper_map_one_read(SeedTable *seed_tab, unsigned char *genome_nucs,
			 long genome_len, MapRead *read) {
    
  /* find best seed (one with fewest matches) on fwd and rev strands */
  read->fwd_seed.len = seed_tab->seed_len;
  read->rev_seed.len = seed_tab->seed_len;
  read->map_code = MAP_CODE_NONE;

  /* try to align fwd strand of read to genome */
  get_best_seed(seed_tab, read->fwd_nucs, read->len,
		&read->fwd_seed);

  /* fprintf(stderr, "using fwd seed from position %u with %u matches\n", */
  /*  	  read->fwd_seed.read_idx, read->fwd_seed.n_match); */

  align_read(genome_nucs, genome_len, read, STRAND_FWD);

  if(read->map_code == MAP_CODE_NONE || read->map_code == MAP_CODE_UNIQUE) {
    /* try to align rev strand of read to genome */
    get_best_seed(seed_tab, read->rev_nucs, read->len,
		  &read->rev_seed);

    /* fprintf(stderr, "using rev seed from position %u with %u matches\n", */
    /* 	    read->rev_seed.read_idx, read->rev_seed.n_match); */

    align_read(genome_nucs, genome_len, read, STRAND_REV);
  }
}





/**
 * Reads all chromosomes into a single long concatenated array
 */
unsigned char *mapper_read_seqs(ChrTable *chr_tab, char **fasta_files, 
				int n_fasta_files) {  
  int chr_idx, i, n_seq;
  long offset;
  unsigned char *nucs;
  Seq *seq;
  gzFile gzf;
  char *mem_str;
  
  mem_str = util_long_to_comma_str(chr_tab->total_chr_len);
  fprintf(stderr, "allocating %s bytes of memory for genome sequence\n",
	  mem_str);
  my_free(mem_str);

  /* allocate enough memory to hold all chromosomes, set to N */
  nucs = my_new(unsigned char, chr_tab->total_chr_len);
  memset(nucs, NUC_N, chr_tab->total_chr_len);

  fprintf(stderr, "reading genome sequence\n");

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
    my_warn("%s:%d: only read %u sequences, but %u are specified in "
	    "chromosome table", __FILE__, __LINE__, n_seq, chr_tab->n_chr);
  }
  
  return nucs;
}


