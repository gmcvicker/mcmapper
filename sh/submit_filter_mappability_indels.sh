#!/bin/sh
#
#
# preserve shell environment vars:
#$ -V
#
# run a job-array consisting of following range of task ids
#$ -t 1-50
#
# write output and err files to following dirs
#$ -o /mnt/lustre/home/gmcvicker/sge/out
#$ -e /mnt/lustre/home/gmcvicker/sge/err
#
# jobs need 8GB of mem:
#$ -l h_vmem=8g

DATA_DIR=/data/share/10_IND_v2

TASK_FILE=$DATA_DIR/10_IND_tasklist.txt

DATA_TYPE=`sed -n ${SGE_TASK_ID}p $TASK_FILE | awk '{print $2}'`
INDIVIDUAL=`sed -n ${SGE_TASK_ID}p $TASK_FILE | awk '{print $3}'`
READ_LEN=`sed -n ${SGE_TASK_ID}p $TASK_FILE | awk '{print $4}'`


INPUT_DIR=$DATA_DIR/filtered/rmdup/$DATA_TYPE/$INDIVIDUAL

OUTPUT_DIR=$DATA_DIR/filtered/mappability/$DATA_TYPE/$INDIVIDUAL
mkdir -p $OUTPUT_DIR

SCRIPT=$HOME/proj/mapper/python/filter_mappability_indels.py


# Align reads, output suffix array (SA) coordinates
echo "filtering reads" 1>&2
echo "python $SCRIPT $READ_LEN $OUTPUT_DIR $INPUT_DIR/$INDIVIDUAL*.txt.gz"  1>&2
python $SCRIPT $READ_LEN $OUTPUT_DIR $INPUT_DIR/$INDIVIDUAL*.txt.gz
if [ "$?" -ne "0" ]; then
    echo "read filtering failed" 1>&2
    exit 1
fi

echo "done" 1>&2

