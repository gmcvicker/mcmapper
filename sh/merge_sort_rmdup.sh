head=/mnt/lustre/data/share/10_IND_v2

for type in PolII H3K4me3 H3K4me1 H3K27ac H3K27me3
  do
  for ind in 18505 18507 18508 18516 18522 19141 19193 19204 19238 19239
    do
    for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
      do
      echo "zcat ${head}/mapper/mapped_reads/${type}/${type}_H_${ind}*/chr${chr}.mapped.txt.gz | sort -k 3n > ${head}/filtered/${type}/${ind}.chr${chr}.temp; \
python ../../rmdup.py ${head}/filtered/${type}/${ind}.chr${chr}.temp ${head}/filtered/${type}/${ind}.chr${chr}.txt; \
gzip -f ${head}/filtered/${type}/${ind}.chr${chr}.txt; \
rm ${head}/filtered/${type}/${ind}.chr${chr}.temp " | qsub -l h_vmem=4g -cwd -V
    done
  done
done
