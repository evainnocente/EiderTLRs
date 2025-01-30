module load bwa
module load samblaster
module load samtools

# $SLURM_ARRAY_TASK_ID is the id of the current script copy

SAMPLELIST="tr_fl_${SLURM_ARRAY_TASK_ID}.txt"

for SAMPLE in `cat $SAMPLELIST`; do

echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

bwa mem -t 8 ../refgenome/ref_genome.fa ../01-adapter_trimming/${SAMPLE}_R1.fastq ../01-adapter_trimming/${SAMPLE}_R2.fastq | samblaster --removeDups | samtools view -h -b -@8 -o ./${SAMPLE}_aligned.bam

done
