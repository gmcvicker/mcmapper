
#include <stdlib.h>
#include <string.h>

#include "ambi.h"
#include "err.h"
#include "nuc.h"
#include "memutil.h"
#include "util.h"



/**
 * return TRUE if provided nucleotides contain an
 * ambiguous nucleotide
 */
int ambi_has_ambi(unsigned char *nucs, int len) {
  int i;
  for(i = 0; i < len; i++) {
    if(nuc_is_ambi(nucs[i])) {
      return TRUE;
    }
  }
  return FALSE;
}



/**
 * Given an ambiguous nucleotide, provides two non-ambiguous nucleotides
 * that it corresponds to
 */
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
    
      if(n == 0 || n == max) {
	/* we are out of arrays */
	return 0;
      }

      /* copy the nucleotide array, set this position 
       * to the other unambiguous nuc, and make recursive 
       * call to this function
       */
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




/**
 * returns TRUE if provided nucleotides match, taking into account
 * that either one can be an ambiguity code.
 */
int ambi_nucs_match(unsigned char nuc1, unsigned char nuc2) {
  unsigned char n1, n2;

  if(nuc_is_ambi(nuc1)) {
    if(nuc_is_ambi(nuc2)) {
      /* Both are ambiguous. Since we only allow ambi codes representing
       * two nucleotides, ambiguity codes must match 
       */
      return nuc1 == nuc2;
    }

    /* nuc1 is ambiguous, check if either possible base matches nuc2 */
    ambi_get_nucs(nuc1, &n1, &n2);

    return ((n1 == nuc2) || (n2 == nuc2));
  } else if(nuc_is_ambi(nuc2)) {
    ambi_get_nucs(nuc2, &n1, &n2);

    /* nuc2 is ambiguous, check if either possible base matches nuc1 */
    return ((n1 == nuc1) || (n2 == nuc1));
  }

  /* neither nucleotide is ambiguous */
  return nuc1 == nuc2;  
}
