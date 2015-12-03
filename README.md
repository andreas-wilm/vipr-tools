simple-contig-joiner
====================

Stitch contigs together by means of Mummer's nucmer and fill gap
with given reference sequence. 

Requirements:
- Python 2.7
- Mummer3: http://mummer.sourceforge.net/
- Testing requires Biopython

Limitations (see TODO.md):
- Currently for small genomes only only (reads all sequences into memory)
- No support for multiple chromosomes
