
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

#include "util.h"
#include "seq.h"


int main(int argc, char **argv) {
  char **fasta_files;
  int n_fasta_files, i;
  Seq *seq;
  gzFile gzf;

  if(argc < 2) {
    fprintf(stderr, "usage: %s chr1.fa.gz [chr2.fa.gz [...]] "
	    "> chromInfo.txt\n", argv[0]);
    exit(2);
  }

  n_fasta_files = argc-1;
  fasta_files = &argv[1];
  seq = seq_new();
  
  for(i = 0; i < n_fasta_files; i++) {
    gzf = util_must_gzopen(fasta_files[i], "rb");
    while(seq_read_fasta_record(seq, gzf)) {
      fprintf(stdout, "%s %ld\n", seq->name, seq->len);
    }
    gzclose(gzf);
  }
  
  seq_free(seq);

  return 0;
}
