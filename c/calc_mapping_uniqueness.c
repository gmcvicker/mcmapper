
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



gzFile get_out_file(const char *output_dir, Chromosome *chr) {
  gzFile out_file;
  char *out_path, *dir;

  if(util_str_ends_with(output_dir, "/")) {
    dir = util_str_dup(output_dir);
  } else {
    dir = util_str_concat(output_dir, "/", NULL);
  }

  out_path = util_str_concat(dir, chr->name, ".wig.gz", NULL);

  if(util_file_exists(out_path)) {
    my_err("%s:%d: output file %s already exists", __FILE__, __LINE__,
	   out_path);
  }

  out_file = util_must_gzopen(out_path, "wb");

  fprintf(stderr, "writing to file %s\n", out_path);

  /* write header */
  gzprintf(out_file, "fixedStep chrom=%s start=1 step=1\n", chr->name);

  my_free(dir);
  my_free(out_path);

  return out_file;
}


void map_read(gzFile out_file, Mapper *mapper, MapRead *read, 
	      long read_start) {

  mapper_map_one_read(mapper, read);
  
  if((read->map_code == MAP_CODE_UNIQUE) && 
     (read->map_offset != read_start)) {
    /* sanity check: if read maps uniquely, it must map back to 
     * this location
     */
    my_err("%s:%d: read maps uniquely, but to wrong genomic offset "
	   "(%u != %u)", __FILE__, __LINE__, read_start, 
	   read->map_offset);
  }
  if(read->map_code == MAP_CODE_NONE) {
    if(!read->has_n) {
      /* samity check: read must map at least once */
      my_err("%s:%d: read from genomic offset %u does not map back "
	     "to genome", __FILE__, __LINE__, read_start);	
    }
  }
  
  /* write result to wiggle file */
  gzprintf(out_file, "%d\n", read->map_code);
}


void map_reads(const char *output_dir, Mapper *mapper, int read_len) {
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

  if(mapper->chr_tab->n_chr < 1) {
    return;
  }
  chr_idx = 0;
  chr = &mapper->chr_tab->chr_array[chr_idx];
  chr_end_offset = mapper->chr_tab->offset[chr_idx] + chr->len - 1;

  /* get first output file, write header */
  fprintf(stderr, "%s\n", chr->name);
  out_file = get_out_file(output_dir, chr);
  
  for(read_start = 0; read_start < mapper->genome_len; read_start++) {
    read_end = read_start + read_len - 1;
    
    if(read_start > chr_end_offset) {
      /* need to go to next chromosome */
      chr_idx += 1;
      if(chr_idx > mapper->chr_tab->n_chr) {
	my_err("%s:%d read start past end of last chromosome",
	       __FILE__, __LINE__);
      }
      chr = &mapper->chr_tab->chr_array[chr_idx];
      chr_end_offset = mapper->chr_tab->offset[chr_idx] + chr->len - 1;

      fprintf(stderr, "%s\n", chr->name);

      /* close old, open new output file */
      gzclose(out_file);
      out_file = get_out_file(output_dir, chr);
    }

    if(read_end > chr_end_offset) {
      /* end of read is past end of chromosome, report as not mapped */
      gzprintf(out_file, "%d\n", MAP_CODE_NONE);
    } else {
      /* map the read to the genome */

      /* copy nucleotide sequence and get revcomp of read */
      memcpy(read.fwd_nucs, &mapper->genome_nucs[read_start], read.len);
      memcpy(read.rev_nucs, &mapper->genome_nucs[read_start], read.len);
      nuc_ids_revcomp(read.rev_nucs, read.len);

      map_read(out_file, mapper, &read, read_start);
    }
  }

  gzclose(out_file);
  my_free(read.fwd_nucs);
  my_free(read.rev_nucs);
}



int main(int argc, char **argv) {
  char **fasta_files, *seed_index_file, *chrom_info_file, *output_dir;
  int n_fasta_files, read_len, n_mismatch;
  SeedTable *seed_tab;
  ChrTable *chr_tab;
  Mapper *mapper;

  if(argc < 5) {
    fprintf(stderr, "usage: %s <seed_index_file> <chromInfo.txt> <read_len>"
	    " <n_mismatch> <output_dir> <chr1.fa.gz> [<chr2.fa.gz [...]]\n",
	    argv[0]);
    exit(2);
  }
  
  seed_index_file = argv[1];
  chrom_info_file = argv[2];
  read_len = util_parse_long(argv[3]);
  n_mismatch = util_parse_long(argv[4]);
  output_dir = argv[5];
  fasta_files = &argv[6];
  n_fasta_files = argc - 6;
  
  chr_tab = chr_table_read(chrom_info_file);
  
  fprintf(stderr, "reading seed index\n");
  seed_tab = seed_table_read(seed_index_file);

  mapper = mapper_init(seed_tab, chr_tab, fasta_files, n_fasta_files, 
		       read_len, n_mismatch);

  fprintf(stderr, "mapping reads\n");
  map_reads(output_dir, mapper, read_len);
  
  seed_table_free(seed_tab);
  chr_table_free(chr_tab);
  mapper_free(mapper);

  return 0;
}
