module load samtools

SAMPLELIST="tr_fl_${SLURM_ARRAY_TASK_ID}.txt"

for SAMPLE in `cat $SAMPLELIST`; do

## Convert to bam file for storage (including all the mapped reads)
  # -F 4 means retain only mapped reads
        samtools view -bS -F 4 ../02-alignment/${SAMPLE}_aligned.bam  > ./mapped/${SAMPLE}_mapped.bam

## Filter the mapped reads (to only retain reads with high mapping quality)
  # Filter bam files to remove poorly mapped reads (non-unique mappings and mappings with a quality score < 20)
        samtools view -h -q 20 ./mapped/${SAMPLE}_mapped.bam | samtools view -buS - | samtools sort -o ./sorted/${SAMPLE}_sorted_minq20.bam
        samtools index ./sorted/${SAMPLE}_sorted_minq20.bam

  # output depth files
        samtools depth -aa ./sorted/${SAMPLE}_sorted_minq20.bam | cut -f 3 | gzip > ./depth_reps/${SAMPLE}_sorted_minq20.bam.depth.gz

done
