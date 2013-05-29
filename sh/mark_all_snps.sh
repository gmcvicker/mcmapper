#!/bin/sh

# directory where reference genome  *.fa.gz files reside:
FASTA_INPUT_DIR=$HOME/data/seq/hg19_nb

# directory with IMPUTE legend files, giving SNP positions and alleles
SNP_INPUT_DIR=/data/share/10_IND_v2/IMPUTE/ALL_1000G_phase1integrated_v3_impute_macGT1

# prefix for IMPUTE legend files (what comes before chromosome name)
SNP_PREFIX=ALL_1000G_phase1integrated_v3_
# postfix for IMPUTE legend files (what comes a
SNP_POSTFIX=_impute_macGT1.legend.gz

# output directory to write marked-up fasta files to:
FASTA_OUTPUT_DIR=/data/share/10_IND_v2/mapper/input_files/hg19_with_snps

# path to binary that marks up fasta files
MARK_SNPS=../c/mark_snps

mkdir -p $FASTA_OUTPUT_DIR

# loop over all .fa.gz files in input directory
for FASTA_IN in `ls $FASTA_INPUT_DIR/chr*.fa.gz`
do
    # extract chromosome name from FASTA file
    CHR_NAME=`echo $FASTA_IN | perl -ne '/.*\/(.*).fa.gz$/; print "$1\n"'`
    echo $CHR_NAME

    SNP_FILE=$SNP_INPUT_DIR/${SNP_PREFIX}${CHR_NAME}${SNP_POSTFIX}

    if [ -f $SNP_FILE ]; then
	FASTA_OUT=$FASTA_OUTPUT_DIR/${CHR_NAME}.fa.gz

	# run the program to make a new fasta file with SNPs
	$MARK_SNPS $FASTA_IN $FASTA_OUT $SNP_FILE

	if [ $? -gt 0 ]; then
	    # something went wrong
	    echo "mark_snps failed for $CHR_NAME: giving up" >&2
	    exit
	fi

    else
	# the random chromosomes etc. do not have SNP files
	echo "no SNP file found for $CHR_NAME: skipping" >&2
    fi
done

