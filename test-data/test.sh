md5=$(which md5sum 2>/dev/null || which md5)
seq=$(which seq 2>/dev/null || which gseq)
tests=20
ref=DENV_1.fa
for i in $($seq 1 $tests); do
    python ./create_contigs.py $ref > rand.fa;
   ../gapped_contig_joiner.py -c rand.fa -r $ref -o gapped_contig_joiner_out | sreformat fasta - > out.fa;
   for f in $ref out.fa; do
      grep -v '^>' $f  | $md5;
   done  | sort -u | wc -l | grep -cq 1  || exit 1; # should only output one md5
   echo "$i/$tests passed"
done

