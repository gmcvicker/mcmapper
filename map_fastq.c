
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nuc.h"
#include "err.h"
#include "util.h"
#include "memutil.h"
#include "mapper.h"
#include "chr_table.h"
#include "fastq.h"
#include "seed_table.h"


/** 
 * opens one output file for each chromosome, returns aray of output files
 */
gzFile *open_out_files(const char *output_dir, ChrTable *chr_tab) {
  int i;
  char *filename;
  gzFile *out_files;

  out_files = my_new(gzFile, chr_tab->n_chr);
  
  for(i = 0; i < chr_tab->n_chr; i++) {
    if(util_str_ends_with(output_dir, "/")) {
      filename = util_str_concat(output_dir, chr_tab->chr_array[i].name,
				 ".mapped.txt.gz", NULL);
    } else {
      filename = util_str_concat(output_dir, "/", chr_tab->chr_array[i].name,
				 ".mapped.txt.gz", NULL);
    }
    
    if(util_file_exists(filename)) {
      my_err("%s:%d: output file already exists: %s",
	     __FILE__, __LINE__, filename);
    }
    
    out_files[i] = util_must_gzopen(filename, "wb");
  }

  return out_files;
}




/**
 * outputs information about a mapped read to stdout
 */
void write_read(gzFile *out_files, ChrTable *chr_tab, MapRead *read) {
  char read_str[FASTQ_MAX_LINE];
  SeqCoord c;
  gzFile out_f;
  int chr_idx;

  if((read->len + 1) > sizeof(read_str)) {
    my_err("%s:%d: read length too long for buffer", __FILE__, __LINE__);
  }

  if(read->map_strand == STRAND_FWD) {
    nuc_ids_to_str(read_str, read->fwd_nucs, read->len);
  } else if(read->map_strand == STRAND_REV) {
    nuc_ids_to_str(read_str, read->rev_nucs, read->len);
  } else {
    my_err("%s:%d: unknown strand", __FILE__, __LINE__);
  }
    
  /* convert offset to sequence coordinate */
  chr_idx = chr_table_offset_to_coord(chr_tab, read->map_offset, &c);
  

  fprintf(stderr, "writing to out file with index %d\n", chr_idx);

  /* each chromosome has it's own output file */
  out_f = out_files[chr_idx];

  /* write read sequence, and genomic coordinate */
  gzprintf(out_f, "%s %s %ld %c %d\n", 
	   read_str, c.chr->name, c.start, 
	   strand_to_char(read->map_strand), read->map_code);

  fprintf(stderr, "%s %s %ld %c %d\n", 
	  read_str, c.chr->name, c.start, 
	  strand_to_char(read->map_strand), read->map_code);



}


void read_from_fastq_record(MapRead *map_read, FastqRead *fastq_read) {    
  /* copy sequence from fastq record into read for mapping */
  map_read->len = fastq_read->read_len;
  nuc_str_to_ids(map_read->fwd_nucs, fastq_read->line2, map_read->len);
  
  /* create reverse complement of read */
  memcpy(map_read->rev_nucs, map_read->fwd_nucs, map_read->len);
  nuc_ids_revcomp(map_read->rev_nucs, map_read->len);
}




void map_reads(gzFile *output_files, gzFile reads_f,
	       ChrTable *chr_tab, SeedTable *seed_tab, 
	       unsigned char *genome_nucs) {
  FastqRead fastq_read;
  MapRead map_read;
  long warn_count, n_fastq_rec, n_fastq_err;
  long n_map_uniq, n_map_multi, n_map_none, line_num;

  map_read.fwd_nucs = my_new(unsigned char, FASTQ_MAX_LINE);
  map_read.rev_nucs = my_new(unsigned char, FASTQ_MAX_LINE);

  n_map_uniq = n_map_multi = n_map_none = 0;
  warn_count = n_fastq_err = n_fastq_rec = 0;

  line_num = 1;

  /* loop over all records in FASTQ file */
  while(TRUE) {
    long r = 0;

    /* read fastq record from file */
    r = fastq_parse_read(&fastq_read, reads_f);
    
    if(r == FASTQ_END) {
      /* we have reached the end of the file */
      break;
    }

    if(r == FASTQ_ERR) {
      /* this fastq record contains an error */

      if(warn_count < FASTQ_MAX_WARN) {
	warn_count += 1;
	my_warn("%s:%d: skipping invalid fastq record starting on line %ld:\n", 
		__FILE__, __LINE__, line_num);
	fprintf(stderr, "  %s\n  %s\n  %s\n  %s\n", fastq_read.line1, 
		fastq_read.line2, fastq_read.line3, fastq_read.line4);
      }
      n_fastq_err += 1;
    }
    else if(r == FASTQ_OK) {
      n_fastq_rec += 1;
      line_num += 4;

      read_from_fastq_record(&map_read, &fastq_read);

      if((n_fastq_rec % 1000000) == 0) {
	fprintf(stderr, ".");
      }

      /* try to map this read to genome */
      mapper_map_one_read(seed_tab, genome_nucs, chr_tab->total_chr_len, 
			  &map_read);

      if(map_read.map_code == MAP_CODE_NONE) {
	/* read does not map to genome */
	n_map_none += 1;
      }
      else if(map_read.map_code == MAP_CODE_MULTI) {
	/* read maps to multiple genomic locations */
	n_map_multi += 1;
      }
      else if(map_read.map_code == MAP_CODE_UNIQUE) {
	/* read maps to single genomic location */
	n_map_uniq += 1;

	/* output mapped read to file */
	write_read(output_files, chr_tab, &map_read);
      }
      else {
	my_err("%s:%d: unknown mapping code", __FILE__, __LINE__);
      }
    } else {
      my_err("%s:%d: unknown fastq status", __FILE__, __LINE__);
    }
  }

  fprintf(stderr, "fastq errors: %ld\n", n_fastq_err);
  fprintf(stderr, "fastq records (without errors): %ld\n", n_fastq_rec);
  fprintf(stderr, "unmapped reads: %ld\n", n_map_none);
  fprintf(stderr, "uniquely mapping reads: %ld\n", n_map_uniq);
  fprintf(stderr, "multiply mapping reads: %ld\n", n_map_multi);

  my_free(map_read.fwd_nucs);
  my_free(map_read.rev_nucs);
}



int main(int argc, char **argv) {
  char **fasta_files, *reads_file, *seed_index_file;
  char *chrom_info_file, *output_dir;
  int n_fasta_files, i;
  SeedTable *seed_tab;
  ChrTable *chr_tab;
  unsigned char *genome_nucs;
  gzFile *out_files, reads_f;

  if(argc < 5) {
    fprintf(stderr, "usage: %s <seed_index_file> <chromInfo.txt> "
	    "<input_reads.fq.gz> <output_dir> <ref_chr1.fa.gz> [<ref_chr2.fa.gz [...]]\n",
	    argv[0]);
    exit(2);
  }
  
  seed_index_file = argv[1];
  chrom_info_file = argv[2];
  reads_file = argv[3];
  output_dir = argv[4];
  fasta_files = &argv[5];
  n_fasta_files = argc - 5;
  
  chr_tab = chr_table_read(chrom_info_file);

  reads_f = util_must_gzopen(reads_file, "rb");
  out_files = open_out_files(output_dir, chr_tab);

  fprintf(stderr, "reading seed index\n");
  seed_tab = seed_table_read(seed_index_file);

  fprintf(stderr, "reading genome sequence\n");
  genome_nucs = mapper_read_seqs(chr_tab, fasta_files, n_fasta_files);

  fprintf(stderr, "mapping reads\n");
  map_reads(out_files, reads_f, chr_tab, seed_tab, genome_nucs);

  for(i = 0; i < chr_tab->n_chr; i++) {
    gzclose(out_files[i]);
  }

  gzclose(reads_f);
  my_free(genome_nucs);
  seed_table_free(seed_tab);
  chr_table_free(chr_tab);

  return 0;
}
