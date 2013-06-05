
#
# edit the following line to point to the location of the libgenome
# library. The code for this library can be downloaded from
# github: https://github.com/gmcvicker/genome
# 
LIBGENOME=../../genome/c/lib

CFLAGS=-Wall -O3 -I$(LIBGENOME) -g
LFLAGS=-lz -lgenome -L$(LIBGENOME)
CC=gcc

# GLIB_FLAGS=`pkg-config --cflags --libs glib-2.0`

objects=chr_table.o kmer.o seed_table.o ambi.o mapper.o fastq.o seed_finder.o

default: $(objects) build_seed_index calc_mapping_uniqueness map_fastq mark_snps make_chrom_info

$(objects): %.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@


mark_snps.o: mark_snps.c
	$(CC) -c $(CFLAGS) -o mark_snps.o mark_snps.c

mark_snps: $(objects) mark_snps.o
	$(CC) $(LFLAGS) -o mark_snps mark_snps.o $(objects)


make_chrom_info.o: make_chrom_info.c
	$(CC) -c $(CFLAGS) -o make_chrom_info.o make_chrom_info.c

make_chrom_info: $(objects) make_chrom_info.o
	$(CC) $(LFLAGS) -o make_chrom_info make_chrom_info.o $(objects)


build_seed_index.o: build_seed_index.c
	$(CC) -c $(CFLAGS) -o build_seed_index.o build_seed_index.c

build_seed_index: $(objects) build_seed_index.o
	$(CC) $(LFLAGS) -o build_seed_index build_seed_index.o $(objects)


calc_mapping_uniqueness.o: calc_mapping_uniqueness.c
	$(CC) -c $(CFLAGS) -o calc_mapping_uniqueness.o calc_mapping_uniqueness.c

calc_mapping_uniqueness: $(objects) calc_mapping_uniqueness.o
	$(CC) $(LFLAGS) -o calc_mapping_uniqueness calc_mapping_uniqueness.o $(objects)


map_fastq.o: map_fastq.c
	$(CC) -c $(CFLAGS) -o map_fastq.o map_fastq.c

map_fastq: $(objects) map_fastq.o
	$(CC) $(LFLAGS) -o map_fastq map_fastq.o $(objects)


test_kmer.o: test_kmer.c
	$(CC) -c $(CFLAGS) -o test_kmer.o test_kmer.c

test_kmer: $(objects) test_kmer.o
	$(CC) $(LFLAGS) -o test_kmer test_kmer.o $(objects)


test_ambi.o: test_ambi.c
	$(CC) -c $(CFLAGS) -o test_ambi.o test_ambi.c

test_ambi: $(objects) test_ambi.o
	$(CC) $(LFLAGS) -o test_ambi test_ambi.o $(objects)

test_fasta_write.o: test_fasta_write.c
	$(CC) -c $(CFLAGS) -o test_fasta_write.o test_fasta_write.c

test_fasta_write: $(objects) test_fasta_write.o
	$(CC) $(LFLAGS) -o test_fasta_write test_fasta_write.o $(objects)


clean:
	rm -f $(objects) build_seed_index test_kmer test_ambi calc_mapping_uniqueness map_fastq mark_snps
