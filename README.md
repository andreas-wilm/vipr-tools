simple-contig-joiner
====================

Stitch contigs together by means of Mummer's nucmer and fill gaps
with given reference sequence. 

Requirements:
- Python 2.7 or 3.4 (might work on other versions as well)
- Mummer3: http://mummer.sourceforge.net/
- Testing requires Biopython

Limitations (see TODO.md):
- Currently for small genomes only only (reads all sequences into memory)
- No support for multiple chromosomes
