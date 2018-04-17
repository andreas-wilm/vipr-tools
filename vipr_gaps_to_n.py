#!/usr/bin/env python3
"""Plotting of ViPR3 coverage and AF values
"""

#--- standard library imports
#
import os
import gzip
from collections import OrderedDict
from itertools import groupby
import argparse
import textwrap

#--- third-party imports
#
#/

#--- project specific imports
#
# /


__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2018 Genome Institute of Singapore"
__license__ = "The MIT License (MIT)"


# oroginal version in vipr_af_vs_cov.py.  should really use module,
# but that makes software setup for vipr messy
def parse_genomecov(genomecov_gzfile, offset=1):
    """Parse gzipped bedtools genomecov output

    WARNING: offset in file is depends on genomecov option (e.g. -d =
    1 based; -dz = 0 based).  Assuming -d = 1 as default so offsetting
    to 0 here.
    """
    genomecov = OrderedDict()
    with gzip.open(genomecov_gzfile) as fh:
        for line in fh:
            sq, pos, cov = line.decode().rstrip().split("\t")
            pos = int(pos)-offset
            cov = int(float(cov))# float for support of scientific notation
            if sq not in genomecov:
                genomecov[sq] = OrderedDict()
            genomecov[sq][pos] = cov
    return genomecov


def fasta_iter(fh):
    """Given a fasta file. Yield tuples of header, sequence.

    Modified form of Brent Pedersen's version
    https://www.biostars.org/p/710/
    """

    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.__next__()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        yield (header, seq)


def gaps_to_n(fasta_fn, gencov_fn, mincov=1):
    """replaces bases with below coverage threshold with Ns
    """
    genomecov = parse_genomecov(gencov_fn)
    assert len(genomecov) == 1

    with open(fasta_fn) as fh:
        seqs = list(fasta_iter(fh))
    assert len(seqs) == 1, ("Expecting one sequence in %s but got %s" % (
        fasta_fn, len(seqs)))
    seqid, seqstr = seqs[0]

    gencov_seqid = list(genomecov.keys())[0]

    # seqid in gencov and fasta not necessarily 100% the same
    assert gencov_seqid.split(' ')[0] == seqid.split(' ')[0]

    assert min(genomecov[gencov_seqid].values()) >= 0
    assert len(seqstr) < max(genomecov[gencov_seqid].values())

    outseql = list(seqstr)
    for pos in range(len(seqstr)):
        # depending on how coverage was computed we might not have a value in the dict
        cov = genomecov[gencov_seqid].get(pos, 0)
        if cov < mincov:
            outseql[pos] = 'N'
    outseqstr = ''.join(outseql)

    print(">%s\n%s\n" % (seqid, textwrap.fill(outseqstr, width=80)))


def main():
    """main function"""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--fasta", required=True,
                        help="Input fasta file")
    parser.add_argument("-c", "--gencov", required=True,
                        help="bedtools genomecov file (gzipped)")
    args = parser.parse_args()

    assert os.path.exists(args.fasta)
    assert os.path.exists(args.gencov)

    gaps_to_n(args.fasta, args.gencov)


if __name__ == "__main__":
    main()
