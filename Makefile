
CFLAGS=-Wall -O3 -I../genome/c/lib -g
LFLAGS=-lz -lgenome -L../genome/c/lib
CC=gcc

# not sure if GLIB will be needed yet or not...
GLIB_FLAGS=`pkg-config --cflags --libs glib-2.0`

objects=chr_table.o kmer.o seed_table.o

default: $(objects) build_seed_index calc_mapping_uniqueness test_kmer

$(objects): %.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@


build_seed_index.o: build_seed_index.c
	$(CC) -c $(CFLAGS) -o build_seed_index.o build_seed_index.c

build_seed_index: $(objects) build_seed_index.o
	$(CC) $(LFLAGS) -o build_seed_index build_seed_index.o $(objects)


calc_mapping_uniqueness.o: calc_mapping_uniqueness.c
	$(CC) -c $(CFLAGS) -o calc_mapping_uniqueness.o calc_mapping_uniqueness.c

calc_mapping_uniqueness: $(objects) calc_mapping_uniqueness.o
	$(CC) $(LFLAGS) -o calc_mapping_uniqueness calc_mapping_uniqueness.o $(objects)


test_kmer.o: test_kmer.c
	$(CC) -c $(CFLAGS) -o test_kmer.o test_kmer.c

test_kmer: $(objects) test_kmer.o
	$(CC) $(LFLAGS) -o test_kmer test_kmer.o $(objects)




clean:
	rm -f $(objects) build_seed_index
