#!/bin/sh
#
#
# preserve shell environment vars:
#$ -V
#
# run a job-array consisting of following range of task ids
#$ -t 1-62
#
# write output and err files to following dirs
#$ -o /mnt/lustre/home/gmcvicker/sge/out
#$ -e /mnt/lustre/home/gmcvicker/sge/err
#
# jobs need 20GB of mem:
#$ -l h_vmem=20g


DATA_DIR=/data/share/10_IND_v2

INPUT_DIR=$DATA_DIR/mapper/input_files

TASK_FILE=$DATA_DIR/fastq/file_list.txt

DATA_TYPE=`sed -n ${SGE_TASK_ID}p $TASK_FILE | awk '{print $1}'`
FILE_NAME=`sed -n ${SGE_TASK_ID}p $TASK_FILE | awk '{print $2}'`

FASTQ_PATH=$DATA_DIR/fastq/$DATA_TYPE/$FILE_NAME

LANE_ID=`echo $FILE_NAME | sed s/.fastq.gz//`
OUTPUT_DIR=$DATA_DIR/mapper/mapped_reads/$DATA_TYPE/$LANE_ID
mkdir -p $OUTPUT_DIR

MAPPER=$HOME/proj/mapper/map_fastq

SEED_FILE=$INPUT_DIR/seed_idx.13bp.hg18.gz
CHROM_FILE=$INPUT_DIR/chromInfo.hg18.txt
FASTA_DIR=$INPUT_DIR/hg18_with_snps


# Align reads, output suffix array (SA) coordinates
echo "mapping reads" 1>&2
echo "$MAPPER $SEED_FILE $CHROM_FILE $FASTQ_PATH $OUTPUT_DIR $FASTA_DIR/chr*.fa.gz" 1>&2
$MAPPER $SEED_FILE $CHROM_FILE $FASTQ_PATH $OUTPUT_DIR $FASTA_DIR/chr*.fa.gz
if [ "$?" -ne "0" ]; then
    echo "read mapping failed" 1>&2
    exit 1
fi

echo "done" 1>&2

