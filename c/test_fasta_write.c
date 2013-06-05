#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>

#include "util.h"
#include "memutil.h"
#include "seq.h"
#include "nuc.h"

#define MAX_AMBI 64

int main(int argc, char **argv) {
  Seq *seq;
  char *seq_str;
  gzFile f;

  if(argc != 2) {
    fprintf(stderr, "usage: %s <nuc_str>\n", argv[0]);
    exit(2);
  }

  seq = seq_new();
  seq_read_str(seq, argv[1]);

  fprintf(stderr, "seq->len=%ld\n", seq->len);
  seq_str = seq_get_seqstr(seq);
  fprintf(stderr, "%s\n", seq_str);
  my_free(seq_str);

  fprintf(stderr, "reading from fasta\n");
  f = gzopen("test.fa.gz", "wb");
  seq_write_fasta_record(seq, f);
  gzclose(f);

  fprintf(stderr, "writing to fasta\n");
  f = gzopen("test.fa.gz", "rb");
  seq_read_fasta_record(seq, f);
  gzclose(f);

  fprintf(stderr, "seq->len=%ld\n", seq->len);
  seq_str = seq_get_seqstr(seq);
  fprintf(stderr, "%s\n", seq_str);
  my_free(seq_str);


  seq_free(seq);
  
  return 0;
}
