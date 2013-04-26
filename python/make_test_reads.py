
import genome.seq

def write_fastq(read_seq, name):
    print "@%s" % name
    print read_seq
    print "+"
    print "".join(["B"] * len(read_seq))


# 28bp:    
# read_seq = "ATGTGTCGGGGGTACAGGGCCAAGCGGA"

# 50 bp:
read_seq = "ATGTGTCGGGGGTACAGGGCCAAGCGGAGCAAGATCCCTTGTGTCCTGTG"
read_array = list(read_seq)

write_fastq(read_seq, "no_mismatch")

# write reads with a single mismatch at each position
for i in range(len(read_seq)):
    # add mismatch at this position
    read_array[i] = genome.seq.revcomp(read_array[i])

    new_seq = "".join(read_array)

    write_fastq(new_seq, "mismatch_at_%d" % (i+1))
    
    # set back to old base
    read_array[i] = read_seq[i]


    
# repeat, but with reverse complement read
read_seq = genome.seq.revcomp(read_seq)
read_array = list(read_seq)

write_fastq(read_seq, "RC_no_mismatch")

for i in range(len(read_seq)):
    # add mismatch at this position
    read_array[i] = genome.seq.revcomp(read_array[i])

    new_seq = "".join(read_array)

    write_fastq(new_seq, "RC_mismatch_at_%d" % (i+1))
    
    # set back to old base
    read_array[i] = read_seq[i]


# repeat but introduce second mismatch
# read_seq = "ATtTGTCGGGGGTACAGGGCCAAGCGGA"
read_seq = "ATtTGTCGGGGGTACAGGGCCAAGCGGAGCAAGATCCCTTGTGTCCTGTG"
read_array = list(read_seq)

for i in range(len(read_seq)):
    # add mismatch at this position
    read_array[i] = genome.seq.revcomp(read_array[i])

    new_seq = "".join(read_array)

    write_fastq(new_seq, "mismatch_at_3_and_%d" % (i+1))
    
    # set back to old base
    read_array[i] = read_seq[i]




# repeat but with ambiguity code and one mismatch
# read_seq = "ATSTGTCGGGGGTACAGGGCCAAGCGGA"
read_seq = "ATSTGTCGGGGGTACAGGGCCAAGCGGAGCAAGATCCCTTGTGTCCTGTG"
read_array = list(read_seq)

for i in range(len(read_seq)):
    # add mismatch at this position
    read_array[i] = genome.seq.revcomp(read_array[i])

    new_seq = "".join(read_array)

    write_fastq(new_seq, "ambi_and_mismatch%d" % (i+1))
    
    # set back to old base
    read_array[i] = read_seq[i]
