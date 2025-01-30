module load samtools
touch sites.txt

for file in $(ls -1 ./sorted/*.bam); do
echo -n "$file: " >> sites.txt
samtools view -H ${file} | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}' >> sites.txt
done

# Check if the number of sites are equal across files
awk 'NR>1{a[$NF]++} END{for(b in a) print b}' sites.txt
