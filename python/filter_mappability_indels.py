
import sys
import os.path
import gzip
import numpy as np

import genome.db


SNP_UNDEF = -1
SNP_TRACK_NAME = "impute2/snps"
SNP_INDEX_TRACK_NAME = "impute2/snp_index"
SNP_REF_MATCH_TRACK_NAME = "impute2/snp_match_ref"
DEL_TRACK_NAME = "impute2/deletions"


def guess_chromosome(filename, chrom_dict):
    words = filename.split(".")

    for word in words:
        if word in chrom_dict:
            return chrom_dict[word]

    raise ValueError("Could not guess chromosome name from filename '%s'" %
                     filename)




def is_indel(snp):
    if (len(snp['allele1']) != 1) or (len(snp['allele2'])) != 1:
        return True




    
def filter_reads(out_f, read_str, read_len, chrom, snp_index_array, 
                 snp_tab, snp_ref_match, del_array, map_array,
                 filter_counts):

    words = read_str.rstrip().split()

    seq_str = words[0]
    if len(seq_str) != read_len:
        raise ValueError("expected read len to be %d but got %d" %
                         (read_len, len(seq_str)))
    chrom_name = words[1]
    if chrom_name != chrom.name:
        raise ValueError("expected reads to be mapped to chr %s, "
                         "but got %s" % (chrom_name, chrom.name))

    
    start = int(words[2])
    end = start + read_len - 1

    # TODO: check for indels overlapping read

    if map_array[start-1] != 1:
        # reads from at least one allele do not map uniquely
        filter_counts['MAP'] += 1
        return

    if np.any(del_array[start-1:end]):
        # read overlaps a base that is deleted in non-reference
        filter_counts['DELETION'] += 1
        return

    idx = snp_index_array[start-1:end]

    # get offsets within read of any overlapping SNPs
    read_offsets = np.where(idx != SNP_UNDEF)[0]
    n_overlapping_snps = read_offsets.size

    for offset in read_offsets:
        match_ref = snp_ref_match[idx[offset]]
        if not match_ref:
            # read overlaps a SNP that does not match reference
            filter_counts['SNP_MISMATCH_REF'] += 1
            return

        snp = snp_tab[idx[offset]]
        
        if is_indel(snp):
            # read overlaps an indel (with inserted bases in non-reference)
            filter_counts['INSERTION'] += 1
            return


    filter_counts['KEPT'] += 1

    out_f.write(read_str)

    

        
            
            


        
        
        
        
    



def get_outfile(out_dir, read_path):
    if out_dir.endswith("/"):
        pass
    else:
        out_dir = out_dir + "/"

    read_filename = read_path.split("/")[-1]
        
    out_path = out_dir + read_filename

    if out_path.endswith("gz"):
        pass
    else:
        out_path = out_path + ".gz"

    if os.path.exists(out_path):
        raise IOError("output file already exists: %s\n" % out_path)

    return gzip.open(out_path, "wb")




def main():
    if len(sys.argv) < 4:
        sys.stderr.write("usage: %s <read_len> <out_dir> "
                         "<read_file1> [<read_file2> ...]\n" 
                         % sys.argv[0])
        exit(2)

    read_len = int(sys.argv[1])
    out_dir = sys.argv[2]
    read_filenames = sys.argv[3:]

    gdb = genome.db.GenomeDB()
    chrom_dict = gdb.get_chromosome_dict()

    trackname = "mappability/%dbp/mappability_0_mismatch" % read_len
    map_track = gdb.open_track(trackname)

    snp_track = gdb.open_track(SNP_TRACK_NAME)
    snp_index_track = gdb.open_track(SNP_INDEX_TRACK_NAME)
    ref_match_track = gdb.open_track(SNP_REF_MATCH_TRACK_NAME)
    del_track = gdb.open_track(DEL_TRACK_NAME)

    filter_counts = {'MAP' : 0,
                     'SNP_MISMATCH_REF' : 0,
                     'DELETION' : 0,
                     'INSERTION' : 0,
                     'KEPT' : 0}
    
    for read_filename in read_filenames:
        chrom = guess_chromosome(read_filename, chrom_dict)

        sys.stderr.write("%s\n" % chrom.name)
        
        f = gzip.open(read_filename)

        out_f = get_outfile(out_dir, read_filename)
        
        sys.stderr.write("fetching SNPs\n")
        snp_tab = snp_track.h5f.getNode("/%s" % chrom.name)
        snp_index_array = snp_index_track.get_nparray(chrom)
        snp_ref_match = ref_match_track.h5f.getNode("/%s" % chrom.name)
        
        del_array = del_track.get_nparray(chrom)

        map_array = map_track.get_nparray(chrom)

        sys.stderr.write("filtering reads\n")
        for line in f:
            filter_reads(out_f, line, read_len, chrom, 
                         snp_index_array, snp_tab, snp_ref_match,
                         del_array, map_array, filter_counts)

        f.close()
        out_f.close()

    sys.stderr.write("filtered reads: \n")
    sys.stderr.write("  non-unique mapping: %d\n" % filter_counts['MAP'])
    sys.stderr.write("  overlap SNP that does not match ref: %d\n" % 
                     filter_counts["SNP_MISMATCH_REF"])
    sys.stderr.write("  overlap deletion: %d\n" %
                     filter_counts["DELETION"])
    sys.stderr.write("  overlap insertion: %d\n" % 
                     filter_counts["INSERTION"])
    sys.stderr.write("kept reads: %d\n" % filter_counts["KEPT"])
    

    snp_track.close()
    snp_index_track.close()
    ref_match_track.close()
    del_track.close()
    map_track.close()    
                    

            
main()
