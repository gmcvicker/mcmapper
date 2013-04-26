
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
gzFile *open_out_files(const char *output_dir, ChrTable *chr_tab,
		       gzFile *unmapped_out_file, gzFile *multi_out_file) {
  int i;
  char *filename, *dir;
  gzFile *out_files;

  out_files = my_new(gzFile, chr_tab->n_chr);

  /* add trailing '/' if not present */
  if(util_str_ends_with(output_dir, "/")) {
    dir = util_str_dup(output_dir);
  } else {
    dir = util_str_concat(output_dir, "/", NULL);
  }

  /* open file for each chromosome for reads that map uniquely */
  for(i = 0; i < chr_tab->n_chr; i++) {
    filename = util_str_concat(dir, chr_tab->chr_array[i].name,
			       ".mapped.txt.gz", NULL);
    
    out_files[i] = util_check_gzopen(filename);
    my_free(filename);
  }


  /* open files for multiply-mapped and unmapped reads */
  filename = util_str_concat(dir, "unmapped.txt.gz", NULL);
  *unmapped_out_file = util_check_gzopen(filename);
  my_free(filename);

  filename = util_str_concat(dir, "multi_mapped.txt.gz", NULL);
  *multi_out_file = util_check_gzopen(filename);
  my_free(filename);

  my_free(dir);

  return out_files;
}


/**
 * outputs information about unmapped read 
 */
void write_unmapped_read(gzFile out_file, MapRead *read) {
  char read_str[FASTQ_MAX_LINE];

  if((read->len + 1) > sizeof(read_str)) {
    my_err("%s:%d: read length too long for buffer", __FILE__, __LINE__);
  }

  nuc_ids_to_str(read_str, read->fwd_nucs, read->len);

  /* write read sequence, and genomic coordinate */
  gzprintf(out_file, "%s . . . %d 0\n", read_str, read->map_code);
}


/**
 * outputs information about a mapped read to file
 */
void write_read(gzFile *out_files, ChrTable *chr_tab, MapRead *read,
		int separate_files) {
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

  if(separate_files) {
    /* each chromosome has its own output file */
    out_f = out_files[chr_idx];
  } else {
    /* only one output file */
    out_f = out_files[0];
  }

  /* write read sequence, and genomic coordinate */
  gzprintf(out_f, "%s %s %ld %c %d %d\n", 
	   read_str, c.chr->name, c.start, 
	   strand_to_char(read->map_strand), read->map_code,
	   read->n_mismatch);
}


void read_from_fastq_record(MapRead *map_read, FastqRead *fastq_read) {    
  /* copy sequence from fastq record into read for mapping */
  map_read->len = fastq_read->read_len;
  nuc_str_to_ids(map_read->fwd_nucs, fastq_read->line2, map_read->len);
  
  /* create reverse complement of read */
  memcpy(map_read->rev_nucs, map_read->fwd_nucs, map_read->len);
  nuc_ids_revcomp(map_read->rev_nucs, map_read->len);
}




void map_reads(gzFile *output_files, gzFile multi_out_file, 
	       gzFile unmapped_out_file, gzFile reads_f, Mapper *mapper) {
  FastqRead fastq_read;
  MapRead map_read;
  long warn_count, n_fastq_rec, n_fastq_err;
  long n_map_uniq, n_map_multi, n_map_none;

  map_read.fwd_nucs = my_new(unsigned char, FASTQ_MAX_LINE);
  map_read.rev_nucs = my_new(unsigned char, FASTQ_MAX_LINE);

  n_map_uniq = n_map_multi = n_map_none = 0;
  warn_count = n_fastq_err = n_fastq_rec = 0;

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
	my_warn("%s:%d: skipping invalid fastq record:\n", 
		__FILE__, __LINE__);
	fprintf(stderr, "  %s\n  %s\n  %s\n  %s\n", fastq_read.line1, 
		fastq_read.line2, fastq_read.line3, fastq_read.line4);
      }
      n_fastq_err += 1;
    }
    else if(fastq_read.read_len != mapper->seed_finder_fwd->read_len) {
      /* check that read length is correct */
	warn_count += 1;
	my_warn("%s:%d: specified read length is %u, but got %d, "
		"skipping read\n",  __FILE__, __LINE__, 
		mapper->seed_finder_fwd->read_len, fastq_read.read_len);
	n_fastq_err += 1;
    }
    else if(r == FASTQ_OK) {
      n_fastq_rec += 1;

      read_from_fastq_record(&map_read, &fastq_read);

      if((n_fastq_rec % 1000000) == 0) {
	fprintf(stderr, ".");
      }

      /* try to map this read to genome */
      mapper_map_one_read(mapper, &map_read);

      if(map_read.map_code == MAP_CODE_NONE) {
	/* read does not map to genome */
	n_map_none += 1;
	write_unmapped_read(unmapped_out_file, &map_read);
      }
      else if(map_read.map_code == MAP_CODE_MULTI) {
	/* read maps to multiple genomic locations */
	n_map_multi += 1;
	write_read(&multi_out_file, mapper->chr_tab, &map_read, FALSE);
      }
      else if(map_read.map_code == MAP_CODE_UNIQUE) {
	/* read maps to single genomic location, output mapped read to file */
	n_map_uniq += 1;
	write_read(output_files, mapper->chr_tab, &map_read, TRUE);
      }
      else {
	my_err("%s:%d: unknown mapping code", __FILE__, __LINE__);
      }
    } else {
      my_err("%s:%d: unknown fastq status", __FILE__, __LINE__);
    }
  }

  fprintf(stderr, "\ndone\n");
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
  int n_fasta_files, n_mismatch, read_len, i;
  SeedTable *seed_tab;
  ChrTable *chr_tab;
  Mapper *mapper;
  gzFile *out_files, reads_f, multi_out_file, unmapped_out_file;

  if(argc < 5) {
    fprintf(stderr, "usage: %s <seed_index_file> <chromInfo.txt> "
	    "<read_len> <n_mismatch> "
	    "<input_reads.fq.gz> <output_dir> "
	    "<ref_chr1.fa.gz> [<ref_chr2.fa.gz [...]]\n",
	    argv[0]);
    exit(2);
  }
  
  seed_index_file = argv[1];
  chrom_info_file = argv[2];
  read_len = util_parse_long(argv[3]);
  n_mismatch = util_parse_long(argv[4]);
  reads_file = argv[5];
  output_dir = argv[6];
  fasta_files = &argv[7];
  n_fasta_files = argc - 7;
  
  chr_tab = chr_table_read(chrom_info_file);

  reads_f = util_must_gzopen(reads_file, "rb");
  out_files = open_out_files(output_dir, chr_tab, 
			     &unmapped_out_file, &multi_out_file);

  fprintf(stderr, "reading seed index\n");
  seed_tab = seed_table_read(seed_index_file);

  mapper = mapper_init(seed_tab, chr_tab, fasta_files, n_fasta_files, 
		       read_len, n_mismatch);

  fprintf(stderr, "mapping reads\n");
  map_reads(out_files, multi_out_file, unmapped_out_file, reads_f, mapper);

  for(i = 0; i < chr_tab->n_chr; i++) {
    gzclose(out_files[i]);
  }
  gzclose(multi_out_file);
  gzclose(unmapped_out_file);

  gzclose(reads_f);
  seed_table_free(seed_tab);
  chr_table_free(chr_tab);
  mapper_free(mapper);

  return 0;
}
