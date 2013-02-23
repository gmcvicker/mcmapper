#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "memutil.h"
#include "kmer.h"
#include "nuc.h"
#include "ambi.h"

#define MAX_AMBI 64

int main(int argc, char **argv) {
  int i, len, n_unambig;
  unsigned char nucs[1023];
  unsigned char **unambig_nucs;
  char unambig_nuc_str[1024];

  if(argc != 2) {
    fprintf(stderr, "usage: %s <nuc_str>\n", argv[0]);
    exit(2);
  }
  char *nuc_str = argv[1];

  len = strlen(nuc_str);
  if(len > sizeof(nucs)) {
    len = sizeof(nucs);
  }
  fprintf(stderr, "len=%d\n", len);
  unambig_nucs = my_new(unsigned char *, MAX_AMBI);
  for(i = 0; i < MAX_AMBI; i++) {
    unambig_nucs[i] = my_new(unsigned char, len);
  }

  nuc_str_to_ids(nucs, nuc_str, len);
  fprintf(stderr, "original nucs:\n");
  nuc_ids_to_str(unambig_nuc_str, nucs, len);
  fprintf(stderr, "  %s\n", unambig_nuc_str);  

  n_unambig = ambi_resolve(nucs, len, unambig_nucs, MAX_AMBI);

  if(n_unambig == 0) {
    fprintf(stderr, "too many possible unambiguous seqs\n");
  } else {
    fprintf(stderr, "there are %d unambiguous seqs\n", n_unambig);

    for(i = 0; i < n_unambig; i++) {
      nuc_ids_to_str(unambig_nuc_str, unambig_nucs[i], len);
      fprintf(stderr, "  %s\n", unambig_nuc_str);
    }
  }

  for(i = 0; i < MAX_AMBI; i++) {
    my_free(unambig_nucs[i]);
  }
  my_free(unambig_nucs);
  
  return 0;
}
