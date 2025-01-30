# remove any unusual ls behaviours

unalias ls

echo "Number of files: " > check.txt
ls -la | wc -l >> check.txt

# loop through all checksums
for hash in $(ls *.md5 -1); do
  md5sum -c $hash >> check.txt
done
