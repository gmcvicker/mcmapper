#
#
#

INPUT_DIR=/data/share/10_IND_v2/mapper/input_files
READ_LEN=36
N_MISMATCH=1

OUT_DIR=/data/share/10_IND_v2/mapper/mappability/${READ_LEN}bp/${N_MISMATCH}_mismatch 

mkdir -p OUT_DIR

./calc_mapping_uniqueness $INPUT_DIR/seed_idx.13bp.hg18.gz $INPUT_DIR/chromInfo.hg18.txt $READ_LEN $N_MISMATCH $OUT_DIR /data/share/10_IND_v2/mapper/input_files/hg18_with_snps/chr*.fa.gz
