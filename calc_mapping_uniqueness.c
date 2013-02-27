
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>

#include "mapper.h"
#include "ambi.h"
#include "util.h"
#include "memutil.h"
#include "nuc.h"
#include "chr.h"
#include "chr_table.h"
#include "seed_table.h"
#include "seq.h"



int has_n(unsigned char *nucs, int len) {
  int i;
  for(i = 0; i < len; i++) {
    if(nucs[i] == NUC_N) {
      return TRUE;
    }
  }
  return FALSE;
}


gzFile get_out_file(const char *output_dir, Chromosome *chr) {
  gzFile out_file;
  char *out_path;

  if(util_str_ends_with(output_dir, "/")) {
    out_path = util_str_concat(output_dir, chr->name, ".wig.gz", NULL);
  } else {
    out_path = util_str_concat(output_dir, "/", chr->name, ".wig.gz", NULL);
  }


  if(util_file_exists(out_path)) {
    my_err("%s:%d: output file %s already exists", __FILE__, __LINE__,
	   out_path);
  }

  out_file = util_must_gzopen(out_path, "wb");

  fprintf(stderr, "writing to file %s\n", out_path);

  my_free(out_path);

  return out_file;
}



void map_reads(const char *output_dir, 
	       ChrTable *chr_tab, SeedTable *seed_tab,
	       unsigned char *genome_nucs, 
	       long genome_len, long read_len) {
  long read_start, read_end;
  MapRead read;
  int chr_idx;
  unsigned int chr_end_offset;
  Chromosome *chr;
  gzFile out_file;
  
  /* allocate memory to hold fwd and rev strand reads */
  read.fwd_nucs = my_new(unsigned char, read_len);
  read.rev_nucs = my_new(unsigned char, read_len);
  read.len = read_len;

  if(chr_tab->n_chr < 1) {
    return;
  }
  chr_idx = 0;
  chr = &chr_tab->chr_array[chr_idx];
  chr_end_offset = chr_tab->offset[chr_idx] + chr->len - 1;

  /* get first output file, write header */
  fprintf(stderr, "%s\n", chr->name);
  out_file = get_out_file(output_dir, chr);
  gzprintf(out_file, "fixedStep chrom=%s start=1 step=1\n",
	   chr->name);
  
  for(read_start = 0; read_start < genome_len; read_start++) {
    read_end = read_start + read_len - 1;
    
    if(read_start > chr_end_offset) {
      /* need to go to next chromosome */
      chr_idx += 1;
      if(chr_idx > chr_tab->n_chr) {
	my_err("%s:%d read start past end of last chromosome",
	       __FILE__, __LINE__);
      }
      chr = &chr_tab->chr_array[chr_idx];
      chr_end_offset = chr_tab->offset[chr_idx] + chr->len - 1;

      fprintf(stderr, "%s\n", chr->name);

      /* open new output file, write header */
      if(out_file) {
	gzclose(out_file);
      }
      out_file = get_out_file(output_dir, chr);
      gzprintf(out_file, "fixedStep chrom=%s start=1 step=1\n",
	       chr->name);
    }
    else if(read_end > chr_end_offset) {
      /* end of read past end of chromosome, report as not mapped */
      gzprintf(out_file, "%d\n", MAP_CODE_NONE);
    }
    else if (has_n(&genome_nucs[read_start], read.len)) {
      /* read contains N, report as not mapped  */
      gzprintf(out_file, "%d\n", MAP_CODE_NONE);
    }
    else {
      /* map the read to the genome */

      /* copy nucleotide sequence and get revcomp of read */
      memcpy(read.fwd_nucs, &genome_nucs[read_start], read.len);
      memcpy(read.rev_nucs, &genome_nucs[read_start], read.len);
      nuc_ids_revcomp(read.rev_nucs, read.len);

      mapper_map_one_read(seed_tab, genome_nucs, genome_len, &read);

      if((read.map_code == MAP_CODE_UNIQUE) && 
	 (read.map_offset != read_start)) {
	/* sanity check: if read maps uniquely, it must map back to 
	 * this location!
	 */
	my_err("%s:%d: read maps uniquely, but to wrong genomic offset "
	       "(%u != %u)", __FILE__, __LINE__, read_start, 
	       read.map_offset);
      }
      if(read.map_code == MAP_CODE_NONE) {
	/* samity check: read must map at least once! */
	my_err("%s:%d: read from genomic offset %u does not map back "
	       "to genome", __FILE__, __LINE__, read_start);	
      }

      gzprintf(out_file, "%d\n", read.map_code);
    }
  }

  if(out_file) {
    gzclose(out_file);
  }
  my_free(read.fwd_nucs);
  my_free(read.rev_nucs);
}



int main(int argc, char **argv) {
  char **fasta_files, *seed_index_file, *chrom_info_file, *output_dir;
  int n_fasta_files, read_len;
  SeedTable *seed_tab;
  ChrTable *chr_tab;
  unsigned char *genome_nucs;

  if(argc < 5) {
    fprintf(stderr, "usage: %s <seed_index_file> <chromInfo.txt> <read_len>"
	    " <output_dir> <chr1.fa.gz> [<chr2.fa.gz [...]]\n",
	    argv[0]);
    exit(2);
  }
  
  seed_index_file = argv[1];
  chrom_info_file = argv[2];
  read_len = util_parse_long(argv[3]);
  output_dir = argv[4];
  fasta_files = &argv[5];
  n_fasta_files = argc - 5;
  
  chr_tab = chr_table_read(chrom_info_file);
  
  fprintf(stderr, "reading seed index\n");
  seed_tab = seed_table_read(seed_index_file);

  fprintf(stderr, "reading genome sequence\n");
  genome_nucs = mapper_read_seqs(chr_tab, 
				 fasta_files, n_fasta_files);

  fprintf(stderr, "mapping reads\n");
  map_reads(output_dir, chr_tab, seed_tab, genome_nucs, 
	    chr_tab->total_chr_len, read_len);
  
  my_free(genome_nucs);
  seed_table_free(seed_tab);
  chr_table_free(chr_tab);

  return 0;
}
