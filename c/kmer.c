#include <stdio.h>
#include <stdlib.h>

#include "util.h"
#include "err.h"
#include "nuc.h"
#include "kmer.h"


/**
 * Returns number of possible kmers that are of provided length
 */
unsigned int kmer_num_kmers(int kmer_len) {
  if(kmer_len < 1 || kmer_len > KMER_MAX_LEN) {
    my_err("%s:%d: kmer length (%d) must be between 1 and %d", 
	   __FILE__, __LINE__, kmer_len, KMER_MAX_LEN);
  }

  /* return 4^kmer_len */
  return 1 << (kmer_len << 1);
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
    if(nucs[i] > 3) {
      my_err("%s:%d: can only create kmer ids for unambiguous nucleotides",
	     __FILE__, __LINE__);
    }

    id += nucs[i] << (i<<1);
  }

  return id;  
}



/**
 * Sets the values in an already allocated buffer to the the NUC ids
 * that make up the kmer with the provided id. The provided buffer
 * should be sufficiently large to hold kmer_len ids.
 *
 * TODO: this could probably be made much faster using bitwise operations
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



