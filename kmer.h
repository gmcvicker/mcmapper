#ifndef __KMER_H__
#define __KMER_H__


#define KMER_MAX_LEN 15

unsigned int kmer_num_kmers(int kmer_len);
unsigned int kmer_nucs_to_id(unsigned char *nucs, int kmer_len);
char *kmer_id_to_str(unsigned int id, char *str_buf, int kmer_len);
unsigned char *kmer_id_to_nucs(unsigned int id, unsigned char *nucs_buf, 
			       int kmer_len);

#endif
