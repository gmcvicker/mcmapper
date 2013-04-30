
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "nuc.h"
#include "err.h"
#include "util.h"
#include "memutil.h"
#include "mapper.h"
#include "chr_table.h"
#include "fastq.h"
#include "seed_table.h"

#define READS_FORMAT_FASTQ 1
#define READS_FORMAT_QSEQ 2

#define OUTPUT_TYPE_SINGLE 1
#define OUTPUT_TYPE_MULTI 2

typedef struct {
  /* path to input file containing seed match table */
  char *seed_index_file;

  /* path to input file containing chromosome names/lengths */
  char *chrom_info_file;

  /* length of reads in input file */
  int read_len;

  /* number of mismatches allowed */
  int n_mismatch;

  /* path to input file containing reads */
  char *reads_file;

  /* flag specifying format of read file (FASTQ or QSEQ) */
  int reads_format;

  /* paths to fastq files containing reference genome sequence */
  char **fasta_files;
  int n_fasta_files;

  /* directory to write output files to */
  char *output_dir;

  /* write output to a single file (ordered same as input file)
   * or to separate unique, unmapped, multimapped files?
   */
  int output_type;

} MapOptions;



void write_map_options(FILE *f, MapOptions *opt) {
  int i;
  
  fprintf(f, "# seed_index_file: %s\n", opt->seed_index_file);
  fprintf(f, "# chrom_info_file: %s\n", opt->chrom_info_file);
  fprintf(f, "# read_len: %d\n", opt->read_len);
  fprintf(f, "# n_mismatch: %d\n", opt->n_mismatch);
  fprintf(f, "# reads_file: %s\n", opt->reads_file);
  if(opt->reads_format == READS_FORMAT_QSEQ) {
    fprintf(f, "# reads_format: qseq\n");
  } else if(opt->reads_format == READS_FORMAT_FASTQ) {
    fprintf(f, "# reads_format: fastq\n");
  } else {
    fprintf(f, "# reads_format: ????\n");
  }

  fprintf(f, "# output_dir: %s\n", opt->output_dir);
  
  if(opt->output_type == OUTPUT_TYPE_SINGLE) {
    fprintf(f, "# output_type: single\n");
  } else if(opt->output_type == OUTPUT_TYPE_MULTI) {
    fprintf(f, "# output_type: multi\n");
  } else {
    fprintf(f, "# output_type: ????\n");
  }
  
  fprintf(f, "# fastq_files:");
  for(i = 0; i < opt->n_fasta_files; i++) {
    fprintf(f, " %s", opt->fasta_files[i]);
  }
  fprintf(f, "\n");
}


void usage() {
  fprintf(stderr, "usage: map_fastq OPTIONS "
	  "<ref_fastq1.fa> [<ref_fasta2.fa> ...]\n");
  fprintf(stderr, "  optional flags:\n");
  fprintf(stderr, "    --fastq  read input file is fastq format (default)\n");
  fprintf(stderr, "    --qseq   read input file is qseq format\n");

  fprintf(stderr, "    --out_single  write output to single file and"
	  " keep reads in the same\n"
	          "                  order as input file (default)\n");
  fprintf(stderr, "    --out_multi   write mapped reads to "
	  "separate chromosome files and\n"
	          "                  write unmapped and multi-mapped "
	  "reads to their own files\n");

  fprintf(stderr, "  required arguments:\n");
  fprintf(stderr, "    --chrom_info <path>      "
	  "file with chromosome names and lengths\n");
  fprintf(stderr, "    --read_len <integer>     "
	  "length of reads in input file\n");
  fprintf(stderr, "    --max_mismatch <integer> "
	  "maximum number of read/genome mismatches\n");
  fprintf(stderr, "    --out_dir <path>         "
	  "directory to write output to\n");
  fprintf(stderr, "    --seed_index <path>      "
	  "path to seed index for reference genome\n");
  fprintf(stderr, "    --reads_file <path>      "
	  "path to fastq or qseq file with reads to map\n\n");
}




void check_map_options(MapOptions *opt) {
  if(opt->reads_file == NULL) {
    usage();
    fprintf(stderr, "must provide --reads_file argument\n");
    exit(2);
  }
  if(opt->output_dir == NULL) {
    usage();
    fprintf(stderr, "must provide --output_dir argument\n");
    exit(2);
  }
  if(opt->seed_index_file == NULL) {
    usage();
    fprintf(stderr, "must provide --seed_index argument\n");
    exit(2);    
  }
  if(opt->read_len < 1) {
    usage();
    fprintf(stderr, "must provide --read_len argument that is >= 1\n");
    exit(2);
  }
  if(opt->n_mismatch < 0) {
    usage();
    fprintf(stderr, "must provide --max_mismatch argument that is >= 0\n");
    exit(2);
  }
}


void parse_map_options(MapOptions *opt,  int argc, char **argv) {
  int c;

  static int output_type_flag = OUTPUT_TYPE_SINGLE;
  static int read_format_flag = READS_FORMAT_FASTQ;

  opt->chrom_info_file = NULL;
  opt->read_len = 0;
  opt->n_mismatch = -1;
  opt->reads_file = NULL;
  opt->fasta_files = NULL;
  opt->n_fasta_files = 0;
  opt->output_dir = NULL;

    
  while (1) {
    static struct option long_options[] = {
	/* flags for single or multiple output files */
	{"out_single",  no_argument, 
	 &output_type_flag, OUTPUT_TYPE_SINGLE},
	{"out_multi",   no_argument, 
	 &output_type_flag, OUTPUT_TYPE_MULTI},

	/* flags for read input file format */
	{"fastq", no_argument, &read_format_flag, READS_FORMAT_FASTQ},
	{"qseq",  no_argument, &read_format_flag, READS_FORMAT_QSEQ},

	{"chrom_info",    required_argument, 0, 'c'},
	{"read_len",      required_argument, 0, 'l'},
	{"max_mismatch",  required_argument, 0, 'm'},
	{"out_dir",       required_argument, 0, 'o'},
	{"seed_index",    required_argument, 0, 's'},
	{"reads_file",    required_argument, 0, 'r'},
	{0, 0, 0, 0}
      };

    /* getopt_long stores the option index here. */
    int option_index = 0;
     
    c = getopt_long(argc, argv, "c:l:m:o:s:r:", 
		    long_options, &option_index);
     
    /* Detect the end of the options. */
    if (c == -1)
      break;
     
    switch (c) {
      case 0:
	/* this is a flag, don't have to do anything */
	break;
      case 'c':
	opt->chrom_info_file = util_str_dup(optarg);
	break;
	
      case 'l':
	opt->read_len = util_parse_long(optarg);
	break;
     
      case 'm':
	opt->n_mismatch = util_parse_long(optarg);
	break;
	
      case 'o':
	opt->output_dir = util_str_dup(optarg);
	break;
	
      case 's':
	opt->seed_index_file = util_str_dup(optarg);
	break;

      case 'r':
	opt->reads_file = util_str_dup(optarg);
	break;
	
      case '?':
	/* getopt_long already printed an error message. */
	usage();
	exit(2);
	break;
     
      default:
	fprintf(stderr, "unknown command line option\n");
	usage();
	exit(2);
    }
  }

  opt->output_type = output_type_flag;
  opt->reads_format = read_format_flag;
  
  /* Print any remaining command line arguments (not options). */
  opt->n_fasta_files = argc - optind;
  
  if(opt->n_fasta_files < 1) {
    usage();
    fprintf(stderr, "must provide at least one reference fasta file\n");
    exit(2);
  }
  opt->fasta_files = &argv[optind];

  check_map_options(opt);
}



gzFile open_single_out_file(const char *output_dir, 
			    const char *input_filename) {
  char *prefix, *out_path;
  char *dot;
  gzFile out_f;
    
  /* want to strip off leading dirs, so get ptr to last '/' */
  prefix = strrchr(input_filename, '/');
  if(prefix == NULL) {
    /* there were no '/' chars */
    prefix = util_str_dup(prefix);
  } else {
    /* take off leading directories */
    prefix = util_str_dup(&prefix[1]);
  }
  
  /* look for first '.', strip off everything after that */
  dot = strchr(prefix, '.');
  if(dot != NULL) {
    /* terminate string at '.' */
    dot[0] = '\0';
  }

  if(strlen(prefix) == 0) {
    fprintf(stderr, "Could not parse prefix from input file. "
	    "Using 'mapped' instead\n");
    my_free(prefix);
    out_path = util_str_concat(output_dir, "mapped.txt.gz", NULL);
  } else {
    out_path = util_str_concat(output_dir, prefix, ".txt.gz", NULL);
    my_free(prefix);
  }
  fprintf(stderr, "writing to output file '%s'\n", out_path);
  out_f = util_check_gzopen(out_path);


  my_free(out_path);

  return out_f;
}



/** 
 * opens one output file for each chromosome, returns aray of output files
 */
gzFile *open_multi_out_files(const char *output_dir,
			     ChrTable *chr_tab, gzFile *unmapped_out_file, 
			     gzFile *multi_out_file) {
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
	       gzFile unmapped_out_file, gzFile reads_f, Mapper *mapper,
	       int reads_format, int output_type) {
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
    if(reads_format == READS_FORMAT_FASTQ) {
      r = fastq_parse_read(&fastq_read, reads_f);
    }
    else if(reads_format == READS_FORMAT_QSEQ) {
      r = fastq_parse_qseq_read(&fastq_read, reads_f);
    } 
    else {
      my_err("%s:%d: unknown read format", __FILE__, __LINE__);
    }
    
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
	if(output_type == OUTPUT_TYPE_SINGLE) {
	  write_unmapped_read(output_files[0], &map_read);
	} else {
	  write_unmapped_read(unmapped_out_file, &map_read);
	}
      }
      else if(map_read.map_code == MAP_CODE_MULTI) {
	/* read maps to multiple genomic locations */
	n_map_multi += 1;

	if(output_type == OUTPUT_TYPE_SINGLE) {
	  write_read(output_files, mapper->chr_tab, &map_read, FALSE);
	} else {
	  write_read(&multi_out_file, mapper->chr_tab, &map_read, FALSE);
	}
      }
      else if(map_read.map_code == MAP_CODE_UNIQUE) {
	/* read maps to single genomic location */
	n_map_uniq += 1;
	
	if(output_type == OUTPUT_TYPE_SINGLE) {
	  write_read(output_files, mapper->chr_tab, &map_read, FALSE);
	} else {
	  write_read(output_files, mapper->chr_tab, &map_read, TRUE);
	}
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
  MapOptions opt;
  int i;
  SeedTable *seed_tab;
  ChrTable *chr_tab;
  Mapper *mapper;
  gzFile *out_files, reads_f, multi_out_file, unmapped_out_file, out_file;

  parse_map_options(&opt, argc, argv);
  write_map_options(stderr, &opt);
  
  chr_tab = chr_table_read(opt.chrom_info_file);

  reads_f = util_must_gzopen(opt.reads_file, "rb");

  if(opt.output_type == OUTPUT_TYPE_MULTI) {
    out_files = open_multi_out_files(opt.output_dir, chr_tab, 
				     &unmapped_out_file, &multi_out_file);
  } 
  else if(opt.output_type == OUTPUT_TYPE_SINGLE) {
    out_file = open_single_out_file(opt.output_dir, opt.reads_file);
    out_files = &out_file;
    unmapped_out_file = NULL;
    multi_out_file = NULL;
  }
  else {
    out_files = NULL;
    multi_out_file = unmapped_out_file = NULL;
    my_err("%s:%d: unknown output type\n", __FILE__, __LINE__);
  }

  fprintf(stderr, "reading seed index\n");
  seed_tab = seed_table_read(opt.seed_index_file);

  mapper = mapper_init(seed_tab, chr_tab, opt.fasta_files, 
		       opt.n_fasta_files, 
		       opt.read_len, opt.n_mismatch);

  fprintf(stderr, "mapping reads\n");
  map_reads(out_files, multi_out_file, unmapped_out_file, reads_f, mapper,
	    opt.reads_format, opt.output_type);

  if(opt.output_type == OUTPUT_TYPE_MULTI) {
    for(i = 0; i < chr_tab->n_chr; i++) {
      gzclose(out_files[i]);
    }
    gzclose(multi_out_file);
    gzclose(unmapped_out_file);
  } else {
    gzclose(out_files[0]);
  }

  gzclose(reads_f);
  seed_table_free(seed_tab);
  chr_table_free(chr_tab);
  mapper_free(mapper);

  return 0;
}
