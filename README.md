# Tools for ViPR

For use with ViPR (version>=4) only

## polish-viral-ref.sh

Polish an initial viral reference by iteratively mapping
variant calling and integration. For high coverage samples
downsampling is a good idea.

Requirements:
- bwa
- samtools
- lofreq
- bcfools
- tabix

Limitations:
- Small genomes only!

## simple_contig_joiner.py

Stitch contigs together by means of Mummer's nucmer and fill gaps
with given reference sequence. 

Requirements:
- Python 2.7 or 3
- [Mummer3](http://mummer.sourceforge.net/)
- Testing requires Biopython

Limitations:
- Currently for small genomes only only (reads all sequences into memory)
- No support for multiple chromosomes


## vipr_af_vs_cov_html.py and vipr_af_vs_cov_pdf.py

Plot variant AFs and coverage vs genome in same plot.

Requirements:
- Python 3
- vipr_af_vs_cov_html.py: plotly
- vipr_af_vs_cov_pdf.py: matplotlib-2



## vipr_gaps_to_n.py

Make a polished assembly submission ready by masking positions without coverage to Ns.

Requirements:
- Python 3

