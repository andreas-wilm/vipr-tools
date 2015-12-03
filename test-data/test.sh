#!/bin/bash

# http://www.redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'


md5=$(which md5sum 2>/dev/null || which md5)
seq=$(which seq 2>/dev/null || which gseq)
NUM_TESTS=10
JOINER=../simple_contig_joiner.py
log=$(mktemp -t test.XXXXXX)

echo "Will exit 1 at any stage if tests fail" | tee -a $log
echo "See $log for more info" | tee -a $log
for ref in DENV_1.fa HIV_1.fa; do
    echo "$NUM_TESTS times generating random contigs from $ref and stitching them back together." | tee -a $log
    for i in $($seq 1 $NUM_TESTS); do
        python ./create_contigs.py $ref > rand.fa;
        $JOINER -c rand.fa -r $ref -o joiner_out -v -v > out.fa 2>>$log;
        for f in $ref out.fa; do
            grep -v '^>' $f | tr -d '\n'  | $md5;
        done  | sort -u | wc -l | grep -cq 1  || exit 1; # should only output one md5
        echo "$ref: $i/$NUM_TESTS passed" | tee -a $log
    done
done
echo "All tests passed"
#rm $log

