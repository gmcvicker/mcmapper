#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "kmer.h"
#include "nuc.h"

int main(int argc, char **argv) {
  long i;
  int kmer_len;
  unsigned char nucs[14];

  if(argc != 2) {
    fprintf(stderr, "usage: %s <kmer>\n", argv[0]);
    exit(2);
  }
  char *nuc_str = argv[1];

  for(i = 1; i <= 15; i++) {
    fprintf(stderr, "kmer_len: %ld, n_kmer: %u\n",
	    i, kmer_num_kmers(i));
  }

  for(i = 1; i <= 100000000; i++) {
    kmer_num_kmers(14);
  }

  kmer_len = strlen(nuc_str);
  fprintf(stderr, "%s (%d)\n", nuc_str, kmer_len);
  nuc_str_to_ids(nucs, nuc_str, kmer_len);

  fprintf(stderr, "kmer_id for %s is %u\n",
	  nuc_str, kmer_nucs_to_id(nucs, kmer_len));


  for(i = 1; i <= 100000000; i++) {
    kmer_nucs_to_id(nucs, kmer_len);
  }

  return 0;
}
