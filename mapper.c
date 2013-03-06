
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



void align_read(Mapper *mapper, MapRead *read, MapSeed *seed, 
		char strand, int max_mismatch) {
  long k, match_offset, genome_seed_start, genome_seed_end;
  long genome_read_start, genome_read_end;
  unsigned char *read_nucs;
  unsigned int kmer_id, i, j;
  int n_mismatch;

  read_nucs = (strand == STRAND_FWD) ? read->fwd_nucs : read->rev_nucs;
  
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
      while((n_mismatch <= max_mismatch) && (k < seed->read_idx)) {
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
      while((n_mismatch <= max_mismatch) && (k < read->len)) {
	if(!ambi_nucs_match(mapper->genome_nucs[genome_read_start + k], 
			    read_nucs[k])) {
	  n_mismatch += 1;
	} 
	k++;
      }

      if(n_mismatch > max_mismatch) {
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



void mapper_map_perfect(Mapper *mapper, MapRead *read, char strand) {
  MapSeed seed;
  unsigned char *nucs;

  nucs = (strand == STRAND_FWD) ? read->fwd_nucs : read->rev_nucs;  
  seed.len = mapper->seed_tab->seed_len;

  /* find best seed (one with fewest genomic matches) */
  get_best_seed(mapper->seed_tab, nucs, read->len, &seed);
  align_read(mapper, read, &seed, strand, 0);
}




/**
 * This function is called once we have already failed to map
 * a read without mismatches. It tries to map the read a second
 * time, but allowing one mismatch.
 */
void mapper_map_one_mismatch(Mapper *mapper, MapRead *read, char strand) {
  MapSeed seed;
  unsigned int seed_end, orig_offset;
  int orig_map_code, orig_n_mismatch;
  unsigned char *nucs;

  nucs = (strand == STRAND_FWD) ? read->fwd_nucs : read->rev_nucs;
  seed.len = mapper->seed_tab->seed_len;
  
  if(read->len < seed.len*2) {
    my_err("%s:%d: read is too short (%d bp) (or seeds are too long "
	   "(%u bp) to map two independent seeds\n", __FILE__, __LINE__,
	   read->len, seed.len);
  }

  /* get best seed from first half of read */
  get_best_seed(mapper->seed_tab, nucs, read->len / 2, &seed);

  /* align using the first seed */
  align_read(mapper, read, &seed, strand, 1);

  if(read->map_code == MAP_CODE_MULTI) {
    /* read maps multiple times after using first seed, give up */
    return;
  }
  
  /* Now get best seed that does not overlap first and try to map again.
   * This is to allow for fact that there may be locations where read can 
   * map but there is a single mismatch located in the original seed
   */
  /* record mapping information from first attempt */
  orig_map_code = read->map_code;
  orig_offset = read->map_offset;
  orig_n_mismatch = read->n_mismatch;

  read->map_code = MAP_CODE_NONE;

  /* get new seed from second half of read */
  seed_end = seed.read_idx + seed.len;
  get_best_seed(mapper->seed_tab, &nucs[seed_end], 
		read->len - seed_end, &seed);

  /* update read_idx since we searched for seeds starting in middle of read */
  seed.read_idx += seed_end;

  /* try to map again */
  align_read(mapper, read, &seed, strand, 1);
  
  switch(read->map_code) {
  case(MAP_CODE_MULTI):
    /* read now maps multiple times, give up */
    break;
  case(MAP_CODE_NONE):
    if(orig_map_code == MAP_CODE_UNIQUE) {
      /* could not map with this seed, but original seed mapped */
      read->map_offset = orig_offset;
      read->map_code = MAP_CODE_UNIQUE;
      read->n_mismatch = orig_n_mismatch;
    }
    break;
  case(MAP_CODE_UNIQUE):
    if((orig_map_code == MAP_CODE_UNIQUE) && 
       (read->map_offset != orig_offset)) {
      /* read mapped to two different locations with different seeds */
      read->map_code = MAP_CODE_MULTI;
    }
    break;
  default:
    my_err("%s:%d: unknown map code");
  }
}




void mapper_map_one_read(Mapper *mapper, MapRead *read) {
  /* do not map reads that contain Ns */
  read->has_n = nuc_ids_have_n(read->fwd_nucs, read->len);
  if(read->has_n) {
    read->map_code = MAP_CODE_NONE;
    return;
  }

  /* try to align fwd strand of read to fwd strand of genome */
  read->map_code = MAP_CODE_NONE;
  mapper_map_perfect(mapper, read, STRAND_FWD);

  if(read->map_code == MAP_CODE_MULTI) {
    /* read mapped multiple times, give up */
    return;
  }

  /* now try to align reverse complement of read */
  mapper_map_perfect(mapper, read, STRAND_REV);
  
  if((read->map_code == MAP_CODE_NONE) && mapper->allow_mismatch) {
    /* 
     * read did not map with 0 mismatches, try again but allowing 1 mismatch
     */
    mapper_map_one_mismatch(mapper, read, STRAND_FWD);
    if(read->map_code == MAP_CODE_MULTI) {
      return;
    }

    mapper_map_one_mismatch(mapper, read, STRAND_REV);

    /* sanity check, should be single mismatch if we mapped this time */
    if(read->map_code == MAP_CODE_UNIQUE && read->n_mismatch != 1) {
      my_err("%s:%d: expected read to have 1 mismatch\n", __FILE__, __LINE__);
    }
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
		    int allow_mismatch) {
  Mapper *mapper;

  mapper = my_new(Mapper, 1);
  mapper->seed_tab = seed_tab;
  mapper->chr_tab = chr_tab;

  /* read entire genome sequence */
  fprintf(stderr, "reading genome sequence\n");
  mapper->genome_nucs = mapper_read_seqs(chr_tab, fasta_files, n_fasta_files);
  mapper->genome_len = chr_tab->total_chr_len;

  mapper->allow_mismatch = allow_mismatch;

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
