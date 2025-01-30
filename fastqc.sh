module load fastqc

unalias ls
ls ./*.fastq.gz -1 | xargs -n 8 fastqc -o ./trimmed_reports
