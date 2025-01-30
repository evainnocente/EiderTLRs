module load samtools

touch depth_${SLURM_ARRAY_TASK_ID}.txt

# Get number of sites, provided they are equal across all files
# NSITES=$(awk 'NR==1 {print $NF}' sites.txt) #1195581529 in this case

SAMPLELIST="tr_fl_${SLURM_ARRAY_TASK_ID}.txt"
for SAMPLE in `cat $SAMPLELIST`; do
echo -n "$SAMPLE:" >> depth_${SLURM_ARRAY_TASK_ID}.txt
samtools depth -a ./sorted/${SAMPLE}_sorted_minq20.bam | awk '{sum+=$3} END {print "Average = ",sum/1195581529}' >> depth_${SLURM_ARRAY_TASK_ID}.txt
done
