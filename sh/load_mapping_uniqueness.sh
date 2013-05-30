


SCRIPT=$HOME/proj/genome/python/script/db/create_track.py

DATA_DIR=/data/share/10_IND_v2/mapper/mappability


for READ_LEN in 28 36 
do
    for N_MISMATCH in 0
    do
	echo "read_len:$READ_LEN n_mismatch:$N_MISMATCH" >&2
	TRACK_NAME=mappability/${READ_LEN}bp/mappability_${N_MISMATCH}_mismatch
        python $SCRIPT --dtype uint8 --format wiggle $TRACK_NAME $DATA_DIR/${READ_LEN}bp/${N_MISMATCH}_mismatch/chr*.wig.gz
    done
done
