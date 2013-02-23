
#ifndef __AMBI_H__
#define __AMBI_H__


void ambi_get_nucs(unsigned char ambi_code, unsigned char *nuc1, 
		   unsigned char *nuc2);
 

int ambi_resolve(unsigned char *nucs, int len, 
		 unsigned char **unambig_nucs, int max);


#endif
