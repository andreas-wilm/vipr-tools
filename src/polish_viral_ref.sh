#!/bin/bash
#
# Polish an initial viral (small!) reference by iteratively mapping
# variant calling and integration. For high coverage samples
# downsampling is a good idea
#
# Requires bwa, samtools, lofreq, bcfools and tabix

set -eou pipefail

threads=4


usage() {
    echo "Usage: $(basename $0) -1 R1.fastq.gz -2 R2.fastq.gz -r ref.fa -o out.fa [-t threads=$threads]" 1>&2;
    exit 1;
}


fq1=""
fq2=""
startfa=""
outfa=""
while getopts ":1:2:o:r:t:" o; do
    case "${o}" in
        1)
            fq1=${OPTARG}
            ;;
        2)
            fq2=${OPTARG}
            ;;
        r)
            startfa=${OPTARG}
            ;;
        o)
            outfa=${OPTARG}
            ;;
        t)
            threads=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${fq1}" ] || [ -z "${fq2}" ] || [ -z "${startfa}" ] || [ -z "{outfa}" ]; then
    usage
fi

if [ ! -e "${fq1}" ] || [ ! -e "${fq2}" ] || [ ! -e "${startfa}" ]; then
    echo "FATAL: At least one of the input files does not exist" 1>&2
    exit 1
fi

if [ -e "${outfa}" ]; then
    echo "FATAL: refusing to overwite existing output file ${outfa}" 1>&2
    exit 1
fi

# check for tools
which bwa >/dev/null || exit 1
which samtools >/dev/null || exit 1
which lofreq >/dev/null || exit 1
which bcftools >/dev/null || exit 1
which tabix >/dev/null || exit 1

maxiter=${maxiter:-10}
tmpdir=$(mktemp -d)
echo "Workdir is $tmpdir"

# first pass consensus from given msa
iter=0
consfa=$(mktemp -p $tmpdir cons.XXXXXX.${iter}.fa)
cp $startfa $consfa

# now map, call variants, insert into ref if >50% and repeat if changes were made
while [ 1 ]; do
    let iter=iter+1
    echo "Iteration ${iter} $(date +%Y%m%d-%H%M): indexing"
    samtools faidx $consfa
    bwa index $consfa >& ${consfa}.log

    echo "Iteration ${iter} $(date +%Y%m%d-%H%M): aligning reads"
    bam=$tmpdir/cons.${iter}.bam
    bwa mem -t $threads $consfa $fq1 $fq2 2>${bam}.log | \
        samtools sort -o $bam -T ${bam}.tmp -

    echo "Iteration ${iter} $(date +%Y%m%d-%H%M): processing BAM"
    # FIXME: try split and join for speedup
    proc_bam=$tmpdir/cons_proc.${iter}.bam
    lofreq viterbi -f $consfa $bam | \
        lofreq alnqual -u - $consfa | \
        lofreq indelqual --dindel -f $consfa - | \
        samtools sort -o $proc_bam -T ${proc_bam}.tmp -
    samtools index $proc_bam
    
    echo "Iteration ${iter} $(date +%Y%m%d-%H%M): calling variants"
    vcf=${proc_bam%.bam}.lf.vcf.gz
    lofreq call-parallel --pp-threads ${threads} -f $consfa \
           --call-indels -o $vcf $proc_bam >& ${vcf}.log

    echo "Iteration ${iter} $(date +%Y%m%d-%H%M): filtering variants"
    proc_vcf=${vcf%.vcf.gz}.flt.vcf.gz
    bcftools view -i 'AF>0.5' $vcf -O z -o $proc_vcf

    # stop if filtered vcf is empty
    nlines=$(zcat $proc_vcf | sed  '/^#/d' | wc -l)
    echo "Found $nlines variants >50% AF"
    if [ $nlines -eq 0 ]; then
        echo "Breaking at iteration $iter"
        break
    fi
    
    echo "Iteration ${iter} $(date +%Y%m%d-%H%M): incorporating filtered variants into consensus"
    tabix $proc_vcf
    newconsfa=$(mktemp -p $tmpdir cons.XXXXXX.${iter}.fa)
    rm $newconsfa    
    bcftools consensus -f $consfa $proc_vcf > $newconsfa
    consfa=$newconsfa
    
    if [ $iter == "$maxiter" ]; then
        echo "Maximum number of iterations reached"
        break
    fi
done

cp $consfa $outfa
rm -rf $tmpdir
echo "Created $outfa after $iter iterations"
