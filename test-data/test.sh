#!/bin/bash

# http://www.redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail
IFS=$'\n\t'

NUM_TESTS=10
JOINER=../simple_contig_joiner.py

md5=$(which md5sum 2>/dev/null || which md5)
seq=$(which seq 2>/dev/null || which gseq)
log=$(mktemp -t test.XXXXXX)

SUCCESS_MSG="Done: all tests passed"
echo "Test only successful if final message says: $SUCCESS_MSG" | tee -a $log
echo "Might exit 1 at any stage on failure" | tee -a $log
echo "See $log for more info" | tee -a $log
for ref in DENV_1.fa HIV_1.fa; do
    randfa="rand.fa"
    outfa="out.fa"
    echo "$NUM_TESTS times generating random contigs from $ref and stitching them back together." | tee -a $log
    for i in $($seq 1 $NUM_TESTS); do
        python ./create_contigs.py $ref > $randfa;
        rm $outfa 2>/dev/null
        cmd="$JOINER -c $randfa -r $ref -o $outfa -v -v"
        echo "DEBUG Executing $cmd" >> $log
        eval $cmd 2>>$log;
        # should only output one md5
        num_sqs=$(for f in $ref $outfa; do grep -v '^>' $f | tr -d '\n'  | $md5; done | sort -u | wc -l)
        if [ $(echo $num_sqs) -ne 1 ]; then
            echo "ERROR: seq in $outfa and $ref differ" 1>&2
            echo "ERROR: command was $cmd" 1>&2
            exit 1
        fi
        echo "$ref: $i/$NUM_TESTS passed" | tee -a $log
    done
done
echo $SUCCESS_MSG | tee -a $log

