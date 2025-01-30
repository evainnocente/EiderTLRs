module load fastp

mkdir fastp_reps

# one way to pass sample list = first passed parameter on command line
# SAMPLELIST=$1

# or, just set the file name directly
SAMPLELIST="file_split${SLURM_ARRAY_TASK_ID}.txt"

# loop through sample list
while IFS=$'\t' read -r -a array; do
fastp -r -w 16 -i ../00-raw/${array[0]}_R1.fastq.gz -I ../00-raw/${array[0]}_R2.fastq.gz \
-o ./${array[0]}_trimmed_R1.fastq.gz -O ./${array[0]}_trimmed_R2.fastq.gz \
--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-h ./fastp_reps/${array[0]}.html
done < $SAMPLELIST
