
import sys
import os
import gzip

import genome.seq
import genome.fasta
import genome.db


SNP_UNDEF = -1
SNP_TRACK_NAME = "impute2/snps"
SNP_INDEX_TRACK_NAME = "impute2/snp_index"
SNP_REF_MATCH_TRACK_NAME = "impute2/snp_match_ref"
DEL_TRACK_NAME = "impute2/deletions"



def is_indel(snp):
    if (len(snp['allele1']) != 1) or (len(snp['allele2']) != 1) or \
        (snp['allele1'] == '-') or (snp['allele2'] == '-'):
        return True


def get_ambi_code(allele1, allele2):
    if allele1 == "A":
        if allele2 == "C":
            return "M"
        elif allele2 == "G":
            return "R"
        elif allele2 == "T":
            return "W"
        else:
            sys.stderr.write("ambiguity code for alleles %s/%s unknown"
                             (allele1, allele2))
            return "N"
        
    elif allele1 == "C":
        if allele2 == "A":
            return "M"
        elif allele2 == "G":
            return "S"
        elif allele2 == "T":
            return "Y"
        else:            
            sys.stderr.write("ambiguity code for alleles %s/%s unknown\n" %
                             (allele1, allele2))
            return "N"
L
    elif allele1 == "G":
        if allele2 == "A":
            return "R"
        elif allele2 == "C":
            return "S"
        elif allele2 == "T":
            return "K"
        else:
            sys.stderr.write("ambiguity code for alleles %s/%s unknown\n" %
                             (allele1, allele2))
            return "N"
        
    elif allele1 == "T":
        if allele2 == "A":
            return "W"
        elif allele2 == "C":
            return "Y"
        elif allele2 == "G":
            return "K"
        else:
            sys.stderr.write("ambiguity code for alleles %s/%s unknown\n" %
                             (allele1, allele2))
            return "N"
    else:
        sys.stderr.write("ambiguity code for alleles %s/%s unknown\n" %
                         (allele1, allele2))
        return "N"
        
    
def main():
    gdb = genome.db.GenomeDB()

    if len(sys.argv) != 2:
        sys.stderr.write("usage: %s <output_dir>\n" % sys.argv[0])
        exit(2)

    out_dir = sys.argv[1]
    

    seq_track = gdb.open_track("seq")
    snp_track = gdb.open_track(SNP_TRACK_NAME)

    for chrom in gdb.get_all_chromosomes():
        filename = "%s/%s.fa.gz" % (out_dir, chrom.name)
        
        if os.path.exists(filename):
            sys.stderr.write("WARNING: skipping chromosome %s " 
                             "output file already exists\n" % chrom.name)
            continue
        
        out_f = gzip.open(filename, "wb")
        
        sys.stderr.write("%s\n" % chrom.name)

        node_name = "/%s" % chrom.name
        if node_name in snp_track.h5f:
            snp_tab = snp_track.h5f.getNode(node_name)
        else:
            snp_tab = []
            sys.stderr.write("WARNING: no SNPs for chromosome %s\n" % chrom.name)
            
        seq_array = seq_track.get_nparray(chrom)

        for snp in snp_tab:
            if is_indel(snp):
                continue

            allele1 = snp['allele1']
            allele2 = snp['allele2']
            
            idx = snp['pos'] - 1
            ref_nuc = chr(seq_array[idx])
            
            if allele1 != ref_nuc:
                # allele1 does not match reference, just set to N
                ambi_code = "N"
            else:
                ambi_code = get_ambi_code(allele1, allele2)

            # update nucleotide, by setting to ambi code
            seq_array[idx] = ord(ambi_code)

        sys.stderr.write("making sequence string\n")
        chr_seq = genome.seq.from_nparray(seq_array)

        sys.stderr.write("writing sequence to fasta file\n")
        
        genome.fasta.write_fasta(out_f, chrom.name, chr_seq)
        out_f.close()
    
main()            


