
#ifndef __AMBI_H__
#define __AMBI_H__


int ambi_has_ambi(unsigned char *nucs, int len);


void ambi_get_nucs(unsigned char ambi_code, unsigned char *nuc1, 
		   unsigned char *nuc2);
 

int ambi_resolve(unsigned char *nucs, int len, 
		 unsigned char **unambig_nucs, int max);


int ambi_nucs_match(unsigned char nuc1, unsigned char nuc2);


#endif
