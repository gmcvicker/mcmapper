#include <stdio.h>
#include <stdlib.h>

#include "util.h"
#include "err.h"
#include "nuc.h"
#include "kmer.h"


static unsigned int KMER_POW_ARRAY[KMER_MAX_LEN];


/**
 * Returns the number of kmers there are of length kmer_len
 * by using an array of pre-computed NUM_NUCS^x values.
 */
unsigned int kmer_num_kmers(int kmer_len) {
  static int compute_kmers = TRUE;
  int i;

  if(kmer_len < 0 || kmer_len > KMER_MAX_LEN) {
    my_err("%s:%d: kmer length must be between 0 and %d", 
	   __FILE__, __LINE__, KMER_MAX_LEN-1);
  }

  if(compute_kmers) {
    fprintf(stderr, "computing kmers\n");
    KMER_POW_ARRAY[0] = 1;
    
    for(i = 1; i < KMER_MAX_LEN; i++) {
      KMER_POW_ARRAY[i] = KMER_POW_ARRAY[i-1] * NUM_REAL_NUCS;
    }
    
    compute_kmers = FALSE;
  }
  
  return KMER_POW_ARRAY[kmer_len];
}



/**
 * Given an array of NUC ids representing a sequence of k-nucleotides
 * returns an identifier that can be used to index kmer counts in an
 * array.
 */
unsigned int kmer_nucs_to_id(unsigned char *nucs, int kmer_len) {
  unsigned int id, i;

  id = 0;
  for(i = 0; i < kmer_len; i++) {

    if(!nuc_is_real(nucs[i])) {
      my_err("%s:%d: can only create kmer ids for unambiguous nucleotides",
	     __FILE__, __LINE__);
    }

    id += kmer_num_kmers(i) * nucs[i];
  }

  return id;  
}



/**
 * Sets the values in an already allocated buffer to the the NUC ids
 * that make up the kmer with the provided id. The provided buffer
 * should be sufficiently large to hold kmer_len ids.
 */
unsigned char *kmer_id_to_nucs(unsigned int id, unsigned char *nuc_buf, 
			       int kmer_len) {
  int i, j;
  unsigned int digit_val;
  
  /* start with the high-order digits and move down */
  for(i = kmer_len-1; i >= 0; i--) {
    /* figure out nucleotide value of this digit by taking one which
     * results in largest value that is smaller than or equal to 
     * the entire remaining number
     */
    for(j = NUM_NUCS-1; j >= 0; j--) {
      digit_val = kmer_num_kmers(i) * j;
      if(digit_val <= id) {
	nuc_buf[i] = j;
	id -= digit_val;
	break;
      }
    }
  }

  return nuc_buf;  
}



/**
 * Sets the character string in the already allocated buffer to the
 * the symbols that make up the kmer with the provided id. The provided
 * buffer should be sufficiently large to hold kmer_len+1 characters.
 */
char *kmer_id_to_str(unsigned int id, char *kmer_buf, int kmer_len) {
  int i, j;
  unsigned int digit_val;
  
  kmer_buf[kmer_len] = '\0';

  /* start with the high-order digits and move down */
  for(i = kmer_len-1; i >= 0; i--) {

    /* figure out nucleotide value of this digit by taking one which
     * results in largest value that is smaller than or equal to 
     * the entire remaining number
     */
    for(j = NUM_NUCS-1; j >= 0; j--) {
      digit_val = kmer_num_kmers(i) * j;
      if(digit_val <= id) {
	kmer_buf[i] = nuc_id_to_char(j);
	id -= digit_val;
	break;
      }
    }
  }
  return kmer_buf;
}



