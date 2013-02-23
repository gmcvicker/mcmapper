
#include <stdlib.h>
#include <string.h>

#include "ambi.h"
#include "err.h"
#include "nuc.h"
#include "memutil.h"

int ambi_resolve_recurs(unsigned char *nucs, int len,
			unsigned char **unambig_nucs, int n,
			int max, int nuc_idx);


void ambi_get_nucs(unsigned char ambi_code, unsigned char *nuc1, 
		   unsigned char *nuc2) {
  switch(ambi_code) {
  case(NUC_M):
    *nuc1 = NUC_A;
    *nuc2 = NUC_C;
    break;
  case(NUC_K):
    *nuc1 = NUC_G;
    *nuc2 = NUC_T;
    break;
  case(NUC_R):
    *nuc1 = NUC_A;
    *nuc2 = NUC_G;
    break;
  case(NUC_Y):
    *nuc1 = NUC_C;
    *nuc2 = NUC_T;
    break;
  case(NUC_W):
    *nuc1 = NUC_A;
    *nuc2 = NUC_T;
    break;
  case(NUC_S):
    *nuc1 = NUC_G;
    *nuc2 = NUC_C;
    break;
  default:
    *nuc1 = NUC_N;
    *nuc2 = NUC_N;
    my_err("%s:%d: unknown ambiguity code %d", __FILE__, __LINE__, ambi_code);
  }
}
 


int ambi_resolve_recurs(unsigned char *nucs, int len,
			unsigned char **unambig_nucs, int n,
			int max, int start) {
  unsigned char nuc1, nuc2;
  int i, old_n;

  for(i = start; i < len; i++) {
    if(nuc_is_ambi(nucs[i])) {
      /* this is an ambiguous nucleotide */
      ambi_get_nucs(nucs[i], &nuc1, &nuc2);

      /* set position in nucleotide array to one of two unambiguous nucs */
      unambig_nucs[n-1][i] = nuc1;
      old_n = n;
      n = ambi_resolve_recurs(nucs, len, unambig_nucs, n, max, i+1);
    
      /* copy the nucleotide array, set this position 
       * to the other unambiguous nuc, and make recursive 
       * call to this function
       */
      if(n == 0 || n == max) {
	/* we are out of arrays */
	return 0;
      }

      if(i > 0) {
	memcpy(unambig_nucs[n], unambig_nucs[old_n-1], i);
      }

      unambig_nucs[n][i] = nuc2;
      return ambi_resolve_recurs(nucs, len, unambig_nucs, n+1, max, i+1);
    } else {
      unambig_nucs[n-1][i] = nucs[i];
    }
  }

  return n;
}


/**
 * "resolves" ambiguity codes in the provided nucleotide array
 * populating the pre-allocated buffer provided with all possible
 * combinations of non-ambiguous nucleotides, (up to max).
 * Returns number of non-ambiguous arrays that were populated
 * or 0 on failure.
 */
int ambi_resolve(unsigned char *nucs, int len, 
		 unsigned char **unambig_nucs, int max) {
  return ambi_resolve_recurs(nucs, len, unambig_nucs, 1, max, 0);
}
