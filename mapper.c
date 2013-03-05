
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
  lowest_n_match = seed_table_n_match(seed_tab, read_nucs);
  lowest_idx = 0;

  /* look for seed with lowest number of matches */
  for(i = 1; i < read_len - seed_tab->seed_len + 1; i++) {
    n_match = seed_table_n_match(seed_tab, &read_nucs[i]);

    if(n_match < lowest_n_match) {
      lowest_n_match = n_match;
      lowest_idx = i;
    }
  }

  /* retrieve genomic match locations for seed with lowest matches */
  seed->read_idx = lowest_idx;
  seed_table_lookup(seed_tab, &read_nucs[lowest_idx], &seed->match);
}


void write_match(FILE *f, unsigned char *genome_nucs, SeedTable *seed_tab,
		 MapRead *read, MapSeed *seed, unsigned int seed_offset, 
		 char strand) {

  char genome_str[1024];
  char seed_str[1024];
  char fwd_read_str[1024];
  char rev_read_str[1024];
  long read_genome_start;
  int i;

  fprintf(stderr, "seed_offset: %u, strand=%d\n", seed_offset, strand);

  if(strand == STRAND_FWD) {
    nuc_ids_to_str(seed_str, &read->fwd_nucs[seed->read_idx], seed->len);
  } else {
    nuc_ids_to_str(seed_str, &read->rev_nucs[seed->read_idx], seed->len);
  }
  nuc_ids_to_str(fwd_read_str, read->fwd_nucs, read->len);
  nuc_ids_to_str(rev_read_str, read->rev_nucs, read->len);
  
  read_genome_start = (long)seed_offset - seed->read_idx;
  nuc_ids_to_str(genome_str, &genome_nucs[read_genome_start], read->len);

  fprintf(f, "seed_offset: %u, read_genome_start: %ld\n",
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



void align_read(Mapper *mapper, MapRead *read, MapSeed *seed, char strand) {
  long k, match_offset, genome_seed_start, genome_seed_end;
  long genome_read_start, genome_read_end;
  unsigned char *read_nucs;
  unsigned int kmer_id, i, j;
  int n_mismatch;

  if(strand == STRAND_FWD) {
    read_nucs = read->fwd_nucs;
  }
  else if(strand == STRAND_REV) {
    read_nucs = read->rev_nucs;
  } else {
    read_nucs = NULL;
    my_err("%s:%d: unknown strand", __FILE__, __LINE__);
  }

  read->n_mismatch = 0;

  /* fprintf(stderr, "seed match has %u kmers\n", seed->match.n_kmer); */
  
  /* loop over every kmer for seed (can be >1 because of ambiguity codes) */
  for(i = 0; i < seed->match.n_kmer; i++) {

    /* loop over every genomic match location for this kmer */
    kmer_id = seed->match.kmer_ids[i];
    for(j = 0; j < mapper->seed_tab->n_match[kmer_id]; j++) {
      match_offset = (long)mapper->seed_tab->match[kmer_id][j];
      n_mismatch = 0;

      /* align from where kmer matched, allow up to max_mismatch mismatches */

      /* get genomic offsets for start / end of read */
      genome_read_start = match_offset - seed->read_idx;
      genome_read_end   = genome_read_start + read->len - 1;

      /* genomic offsets of start / end of kmer match */
      genome_seed_start = match_offset;
      genome_seed_end   = match_offset + (long)seed->len - 1;

      if((genome_read_start < 0) || (genome_read_end >= mapper->genome_len)) {
	/* read would overhang end of genome */
	continue;
      }
      
      /* check left end of read (before seed) matches genome, allowing
       * ambiguity codes
       */
      k = 0;
      while((n_mismatch <= mapper->max_mismatch) && (k < seed->read_idx)) {
	if(!ambi_nucs_match(mapper->genome_nucs[genome_read_start + k], 
			    read_nucs[k])) {
	  n_mismatch += 1;
	}
	k++;
      }

      /* check right end of read (after seed) matches genome, allowing
       * ambiguity codes
       */
      k = seed->read_idx + seed->len;
      while((n_mismatch <= mapper->max_mismatch) && (k < read->len)) {
	if(!ambi_nucs_match(mapper->genome_nucs[genome_read_start + k], 
			    read_nucs[k])) {
	  n_mismatch += 1;
	} 
	k++;
      }

      if(n_mismatch > mapper->max_mismatch) {
	/* too many mismatches */
	continue;
      }

      /* read matches! */

      /* write_match(stderr, genome_nucs, read, seed, match_offset, strand); */
      if(read->map_code == MAP_CODE_NONE) {
	/* first match we have observed */
	read->map_code = MAP_CODE_UNIQUE;
	read->map_offset = genome_read_start;
	read->map_strand = strand;
	read->n_mismatch = n_mismatch;
      }
      else if(read->map_code == MAP_CODE_UNIQUE) {
	/* second match we have observed, can quit now since we
	 * know this is not a uniquely-mapping read
	 */
	read->map_code = MAP_CODE_MULTI;
	return;
      }
      else {
	my_err("%s:%d: expected mapping code to be NONE or UNIQUE",
	       __FILE__, __LINE__);
      }
    }
  }
}



void mapper_map_perfect(Mapper *mapper, MapRead *read) {
  MapSeed seed;
  
  seed.len = mapper->seed_tab->seed_len;

  /* find best seed (one with fewest matches) on fwd strands */
  read->map_code = MAP_CODE_NONE;

  /* try to align fwd strand of read to genome */
  get_best_seed(mapper->seed_tab, read->fwd_nucs, read->len, &seed);

  /* fprintf(stderr, "using fwd seed from position %u with %u matches\n", */
  /*   	  seed.read_idx, seed.match.n_match); */

  align_read(mapper, read, &seed, STRAND_FWD);

  if(read->map_code == MAP_CODE_NONE || read->map_code == MAP_CODE_UNIQUE) {

    /* now try to align read to rev strand of read to genome */
    get_best_seed(mapper->seed_tab, read->rev_nucs, read->len, &seed);

    /* fprintf(stderr, "using rev seed from position %u with %u matches\n", */
    /*  	    seed.read_idx, seed.match.n_match); */

    align_read(mapper, read, &seed, STRAND_REV);
  }
}




void mapper_map_one_mismatch(Mapper *mapper, MapRead *read) {
  /** TODO **/

}



void mapper_map_one_read(Mapper *mapper, MapRead *read) {
  read->has_n = nuc_ids_have_n(read->fwd_nucs, read->len);
  if(read->has_n) {
    /* do not map reads that contain Ns */
    read->map_code = MAP_CODE_NONE;
    return;
  }

  if(mapper->max_mismatch == 0) {
    mapper_map_perfect(mapper, read);
  }
  else if(mapper->max_mismatch == 1) {
    mapper_map_one_mismatch(mapper, read);
  }
  else {
    my_err("%s:%d: mapping allowing %d mismatches not currently implemented",
	   __FILE__, __LINE__, mapper->max_mismatch);
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




/**
 * Allocates memory for and creates a new Mapper data structure
 */
Mapper *mapper_init(SeedTable *seed_tab, ChrTable *chr_tab, 
		    char **fasta_files, int n_fasta_files,
		    int max_mismatch) {
  Mapper *mapper;

  mapper = my_new(Mapper, 1);
  mapper->seed_tab = seed_tab;
  mapper->chr_tab = chr_tab;

  /* read entire genome sequence */
  fprintf(stderr, "reading genome sequence\n");
  mapper->genome_nucs = mapper_read_seqs(chr_tab, fasta_files, n_fasta_files);
  mapper->genome_len = chr_tab->total_chr_len;

  mapper->max_mismatch = max_mismatch;

  if(mapper->max_mismatch > 1) {
    my_err("%s:%d: only 0 or 1 mismatches currently implemented",
	   __FILE__, __LINE__);
  }

  return mapper;
}


/**
 * frees memory allocated for read mapper, including genome sequences
 * and seeds (but not including seed 
 */
void mapper_free(Mapper *mapper) {
  my_free(mapper->genome_nucs);
  my_free(mapper);
}
