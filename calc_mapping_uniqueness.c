
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

#define MAP_CODE_NONE 0   /* does not map */
#define MAP_CODE_UNIQUE 1 /* maps uniquely */
#define MAP_CODE_MULTI 2  /* maps multiple times */

typedef struct {
  unsigned int len;
  unsigned int read_idx;
  unsigned int n_match;
  unsigned int *matches;
} Seed;

typedef struct {
  /* fwd nucleotides of read */
  unsigned char *fwd_nucs;

  /* rev nucleotides of read */
  unsigned char *rev_nucs;

  /* length of read */
  unsigned int len;

  Seed fwd_seed;
  Seed rev_seed;

  /* offset and strand giving genomic location that read maps to */
  unsigned int map_offset;
  char map_strand;
  
  /* does this read map uniquely, multiple times, or not at all? */
  int map_code;
} SeqRead;



/**
 * returns true if provided array of nucleotide contains an ambiguity code
 */
int has_ambi_code(unsigned char *nucs, int len) {
  int i;
  for(i = 0; i < len; i++) {
    if(nucs[i] == NUC_N || nuc_is_ambi(nucs[i])) {
      return TRUE;
    }
  }
  return FALSE;
}



/**
 * Given a read, identifies the seed with the fewest number of
 * matches, genome-wide. Returns the number of matches and sets the
 * index into the read and ptr to match offsets
 */
void get_best_seed(SeedTable *seed_tab, unsigned char *read_nucs,
		   unsigned int read_len, Seed *seed) {

  unsigned int i, n_match, lowest_n_match, lowest_idx;
  unsigned int *lowest_matches, *matches;

  if(seed_tab->seed_len > read_len) {
    my_err("%s:%d: seed len (%d) > than read len (%d)", __FILE__, __LINE__,
	   seed_tab->seed_len, read_len);
  }

  /* start with seed at beginning of read */
  lowest_n_match = seed_table_get_matches(seed_tab, read_nucs, &matches);
  lowest_idx = 0;
  lowest_matches = matches;

  /* look for seed with lowest number of matches */
  for(i = 1; i < read_len - seed_tab->seed_len + 1; i++) {
    n_match = seed_table_get_matches(seed_tab, &read_nucs[i], &matches);

    if(n_match < lowest_n_match) {
      lowest_n_match = n_match;
      lowest_idx = i;
      lowest_matches = matches;
    }
  }

  seed->read_idx = lowest_idx;
  seed->matches = lowest_matches;
  seed->n_match = lowest_n_match;
}


void write_match(FILE *f, unsigned char *genome_nucs, SeqRead *read, 
		 Seed *seed, unsigned int match_idx, char strand) {

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
		SeqRead *read, char strand) {  
  long i, j, genome_seed_start, genome_seed_end;
  long genome_read_start, genome_read_end;
  unsigned char *read_nucs;
  Seed *seed;
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
      /* read overhangs end of genome */
      fprintf(stderr, "skipping alignment that would overhang edge "
	      "of genome sequence\n");
      continue;
    }
    
    /* check that left end of read (before seed) matches genome */
    j = 0;
    while(perfect_match && (j < seed->read_idx)) {
      if(genome_nucs[genome_read_start + j] != read_nucs[j]) {
	/* TODO: fix this to allow for ambiguity codes */
	/* does not match */
	perfect_match = FALSE;
      }
      j++;
    }

    /* check that right end of read (after seed) matches genome */
    j = seed->read_idx + seed->len;
    while(perfect_match && (j < read->len)) {
      if(genome_nucs[genome_read_start + j] != read_nucs[j]) {
	/* TODO: fix to allow for ambiguity codes */
	/* does not match */
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



void map_one_read(SeedTable *seed_tab, unsigned char *genome_nucs,
		  long genome_len, SeqRead *read) {
    
  /* find best seed (one with fewest matches) on fwd and rev strands */
  read->fwd_seed.len = seed_tab->seed_len;
  read->rev_seed.len = seed_tab->seed_len;

  get_best_seed(seed_tab, read->fwd_nucs, read->len,
		&read->fwd_seed);

  get_best_seed(seed_tab, read->rev_nucs, read->len,
		&read->rev_seed);

  /* fprintf(stderr, "using fwd seed from position %u with %u matches\n", */
  /*  	  read->fwd_seed.read_idx, read->fwd_seed.n_match); */

  /* fprintf(stderr, "using rev seed from position %u with %u matches\n", */
  /*  	  read->rev_seed.read_idx, read->rev_seed.n_match); */

  read->map_code = MAP_CODE_NONE;

  /* try to align fwd strand of read to genome */
  align_read(genome_nucs, genome_len, read, STRAND_FWD);

  if(read->map_code == MAP_CODE_NONE) {
    my_err("read should map once\n");
  }

  if(read->map_code == MAP_CODE_NONE || read->map_code == MAP_CODE_UNIQUE) {
    /* try to align rev strand of read to genome */
    align_read(genome_nucs, genome_len, read, STRAND_REV);
  }
}


/**
 * outputs information about a mapped read to stdout
 */
void write_read(ChrTable *chr_tab, SeqRead *read) {
  char read_str[1024];
  SeqCoord c;

  if((read->len + 1) > sizeof(read_str)) {
    my_err("%s:%d: read length too long for buffer", __FILE__, __LINE__);
  }

  if((read->map_code == MAP_CODE_MULTI) || 
     (read->map_code == MAP_CODE_UNIQUE)) {
    /* read mapped somewhere */
    if(read->map_strand == STRAND_FWD) {
      nuc_ids_to_str(read_str, read->fwd_nucs, read->len);
    } else if(read->map_strand == STRAND_REV) {
      nuc_ids_to_str(read_str, read->rev_nucs, read->len);
    } else {
      my_err("%s:%d: unknown strand", __FILE__, __LINE__);
    }
    
    /* convert offset to sequence coordinate */
    chr_table_offset_to_coord(chr_tab, read->map_offset, &c);

    fprintf(stdout, "%s %s %ld %c %d\n", 
	    read_str, c.chr->name, c.start, strand_to_char(read->map_strand),
	    read->map_code);
  } else {
    /* read did not map anywhere, use fwd strand of read when reporting */
    nuc_ids_to_str(read_str, read->fwd_nucs, read->len);

    fprintf(stdout, "%s NA NA . %d\n", read_str, read->map_code);
  }
}


void map_reads(ChrTable *chr_tab, SeedTable *seed_tab,
	       unsigned char *genome_nucs, 
	       long genome_len, long read_len) {
  long i;
  SeqRead read;
  
  read.fwd_nucs = my_new(unsigned char, read_len);
  read.rev_nucs = my_new(unsigned char, read_len);
  read.len = read_len;
  
  for(i = 0; i < genome_len - read_len + 1; i++) {
    if(has_ambi_code(&genome_nucs[i], read.len)) {
      /* TODO: fix this so that ambiguity codes can be used */
      continue;
    }

    /* copy nucleotide sequence and get revcomp of read */
    memcpy(read.fwd_nucs, &genome_nucs[i], read.len);
    memcpy(read.rev_nucs, &genome_nucs[i], read.len);
    nuc_ids_revcomp(read.rev_nucs, read.len);

    map_one_read(seed_tab, genome_nucs, genome_len, &read);

    write_read(chr_tab, &read);
  }

  my_free(read.fwd_nucs);
  my_free(read.rev_nucs);
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
  
  fprintf(stderr, "allocating %u bytes of memory for genome sequence\n",
	  chr_tab->total_chr_len);

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
  chr_table_write(stdout, chr_tab);

  fprintf(stderr, "reading seed index\n");
  seed_tab = seed_table_read(seed_index_file);

  fprintf(stderr, "reading genome sequence\n");
  genome_nucs = read_seqs(chr_tab, fasta_files, n_fasta_files);

  fprintf(stderr, "mapping reads\n");
  map_reads(chr_tab, seed_tab, genome_nucs, chr_tab->total_chr_len, 36);
  
  my_free(genome_nucs);
  seed_table_free(seed_tab);
  chr_table_free(chr_tab);

  return 0;
}
