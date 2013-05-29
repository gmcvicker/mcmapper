#include <stdio.h>
#include <string.h>

#include "nuc.h"
#include "memutil.h"
#include "util.h"
#include "seq.h"

#include "ambi.h"

#define MAX_SNP_LINE 1000000


void markup_snps(Seq *seq, gzFile snp_gzf) {
  char *line;
  char allele1[101], allele2[101];
  long line_num, pos;
  unsigned char nuc1, nuc2, ambi_nuc;

  line = my_new(char, MAX_SNP_LINE);

  /* read header */
  gzgets(snp_gzf, line, MAX_SNP_LINE);

  line_num = 1;
  if(!util_str_ends_with(line, "\n")) {
    my_err("%s:%d: header line exceeeds max line length\n",
	   __FILE__, __LINE__);
  }

  while(gzgets(snp_gzf, line, MAX_SNP_LINE)) {
    line_num += 1;

    if((line_num % 100000) == 0) {
      fprintf(stderr, ".");
    }

    if(!util_str_ends_with(line, "\n")) {
      my_err("%s:%d: line %ld exceeeds max line length (%ld)\n",
	     __FILE__, __LINE__, line_num, MAX_SNP_LINE);
    }
    
    sscanf(line, "%*s %ld %100s %100s", &pos, allele1, allele2);
    
    if(pos < 1 || pos > seq->len) {
      my_warn("skipping SNP at position %ld, which is outside of "
	      "chromosome range 1-%ld", pos, seq->len);
      continue;
    }

    if((strlen(allele1) != 1) || (strlen(allele2) != 1)) {
      /* this is an indel--ignore */
      continue;
    }
    if((allele1[0] == '-') || (allele2[0] == '-')) {
      /* this is also an indel--ignore */
      continue;
    }

    /* lookup ambiguity code */
    nuc1 = nuc_char_to_id(allele1[0]);
    nuc2 = nuc_char_to_id(allele2[0]);
    ambi_nuc = ambi_from_nucs(nuc1, nuc2);
    
    if((nuc1 != seq->sym[pos-1]) && (nuc2 != seq->sym[pos-1])) {
      my_warn("neither allele (%c/%c) of SNP at position %ld matches "
	      "reference sequence (%c)", allele1[0], allele2[0], pos,
	      nuc_id_to_char(seq->sym[pos-1]));
    }

    /* update sequence */
    seq->sym[pos-1] = ambi_nuc;
  }  
  fprintf(stderr, "\n");

  my_free(line);
}


int main(int argc, char **argv) {
  char *in_fasta_file, *out_fasta_file, *snp_file;
  gzFile fasta_gzf, snp_gzf;
  Seq *seq;
 

  if(argc != 4) {
    fprintf(stderr, "usage: %s <input.fa.gz> <output.fa.gz> "
	    "<snp_file.legend.gz>\n", argv[0]);
    exit(2);
  }

  in_fasta_file = argv[1];
  out_fasta_file = argv[2];
  snp_file = argv[3];

  if(util_file_exists(out_fasta_file)) {
    my_err("%s:%d: output file '%s' already exists--remove "
	   "file and run again", __FILE__, __LINE__, out_fasta_file);
  }

  /* read input sequence */
  fprintf(stderr, "reading sequence from file %s\n", in_fasta_file);
  fasta_gzf = util_must_gzopen(in_fasta_file, "rb");
  seq = seq_new();
  seq_read_fasta_record(seq, fasta_gzf);
  gzclose(fasta_gzf);

  /* open output file */
  fasta_gzf = util_must_gzopen(out_fasta_file, "wb");

  /* open SNP input file */
  fprintf(stderr, "reading snps from file %s\n", snp_file);
  snp_gzf = util_must_gzopen(snp_file, "rb");

  /* markup SNPs in sequence with ambiguity codes */
  markup_snps(seq, snp_gzf);
  
  /* write new fasta file */
  fprintf(stderr, "writing new sequence to file %s\n", out_fasta_file);
  seq_write_fasta_record(seq, fasta_gzf);

  seq_free(seq);
  gzclose(snp_gzf);
  gzclose(fasta_gzf);

  return 0;
}
