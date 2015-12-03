#!/usr/bin/env python
"""Stitch contigs together by means of Mummer's nucmer and fill gap
with given reference sequence

"""

#--- standard library imports
#
import sys
import os
import logging
import argparse
import subprocess
import tempfile
import shutil
from itertools import groupby
import string
from collections import namedtuple

#--- third-party imports
#
# /

#--- project specific imports
#
# /


__author__ = "Andreas Wilm"
__version__ = "0.2"
__email__ = "wilma@gis.a-star.edu.sg"
__license__ = "WTFPL http://www.wtfpl.net/"


# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.INFO,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


TilingContig = namedtuple('TilingContig', [
    'ref_name', 'ref_start', 'ref_end', 'gap2next',
    'len', 'aln_cov', 'perc_ident', 'ori', 'name'])


def rev_comp(dna):
    """compute reverse complement for dna

    No support for RNA.

    >>> rev_comp("AaGgCcTtNn")
    'nNaAgGcCtT'
    """

    # maketrans doc: "Don't use strings derived from lowercase and
    # uppercase as arguments; in some locales, these don't have the
    # same length. For case conversions, always use str.lower() and
    # str.upper()."
    old_chars = "ACGTN"
    old_chars += str.lower(old_chars)
    replace_chars = "TGCAN"
    replace_chars += str.lower(replace_chars)
    trans = string.maketrans(old_chars, replace_chars)
    return dna.translate(trans)[::-1]


def fasta_iter(fasta_name):
    """
    Given a fasta file. yield tuples of header, sequence

    Author: Brent Pedersen
    https://www.biostars.org/p/710/
    """
    fa_fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fa_fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq
    fa_fh.close()


def cmdline_parser():
    """
    creates argparse instance
    """

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-v', '--verbose',
                        action='count', default=0)
    parser.add_argument('-q', '--quiet',
                        action='count', default=0)
    parser.add_argument("-c", "--c",
                        required=True,
                        dest="fcontigs",
                        help="Input file containing contigs to join (fasta)")
    parser.add_argument("-n", "--dont-fill-with-ref",
                        dest="dont_fill_with_ref",
                        action="store_true",
                        help="Don't fill gaps with reference (keep Ns)")
    parser.add_argument("-r", "--ref",
                        required=True,
                        dest="fref",
                        help="Reference sequence file (fasta)")
    parser.add_argument("-o", "--output",
                        required=True,
                        dest="fout",
                        help="output file (fasta)")
    parser.add_argument("--keep-tmp-files",
                        dest="keep_temp_files",
                        action="store_true",
                        help="Don't delete temp (nucmer etc) files")
    parser.add_argument("--tmp-dir",
                        dest="tmp_dir", # type="string|int|float"
                        help="directory to save temp files in")

    return parser


def run_nucmer(fref, fcontigs, out_prefix, nucmer="nucmer"):
    """Run's nucmer on given reference and query (contigs).
    Returns path to delta file"""

    fdelta = out_prefix + ".delta"
    cmd = [nucmer, fref, fcontigs, '-p', out_prefix]
    LOG.debug("Calling {}".format(cmd))
    try:
        o = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except (subprocess.CalledProcessError, OSError):
        LOG.fatal("The following command failed: {}\n".format(" ".join(cmd)))
        raise
    assert os.path.exists(fdelta), (
        "Expected file '{}' missing. Command was: '{}'. Output was '{}'".format(
            fdelta, ' '.join(cmd), o))
    return fdelta


def run_showtiling(fdelta):
    """Run mummer's show-tiling which creates a pseudo molecule from
    the contigs, filled with Ns. Return path to pseudo molecule
    and tiling file
    """

    fpseudo = fdelta + ".pseudo.fa"
    ftiling = fdelta + ".tiling.txt"
    cmd = ['show-tiling', '-p', fpseudo, fdelta]
    LOG.debug("Calling {}".format(cmd))
    try:
        o = subprocess.check_output(cmd, stderr=subprocess.STDOUT)
    except (subprocess.CalledProcessError, OSError):
        LOG.fatal("The following command failed: {}\n".format(" ".join(cmd)))
        raise
    assert os.path.exists(fpseudo), (
        "Expected file '{}' missing. Command was: '{}'. Output was '{}'".format(
            fpseudo, ' '.join(cmd), o))
    with open(ftiling, 'w') as tile_fh:
        tile_fh.write(o)
    LOG.info("Pseudo molecule written to '{}'".format(fpseudo))
    LOG.info("Tiling info written to '{}'".format(ftiling))
    return fpseudo, ftiling


def parse_tiling(tiling_file):
    """Parse Nucmer's tiling file and yields TilingContig's

    function returns refstart and refstart zero-based half open (i.e.
    in slice notation)

    From Mummers doc: Standard output has an 8 column list per
    mapped contig, separated by the FastA headers of each
    reference sequence. These columns are as follows: [1] start in
    the reference [2] end in the reference [3] gap between this
    contig and the next [4] length of this contig [5] alignment
    coverage of this contig [6] average percent identity of this
    contig [7] contig orientation [8] contig ID.

    """

    with open(tiling_file) as tile_fh:
        for line in tile_fh:
            if line.startswith(">"):
                ref_name = line[1:].split()[0].strip()
                continue

            (ref_start, ref_end, gap2next, contig_len, aln_cov,
             perc_ident, ori, name) = line.strip().split("\t")
            assert ori in "+-"
            (ref_start, ref_end, gap2next, contig_len) = (
                int(x) for x in (ref_start, ref_end, gap2next, contig_len))
            (aln_cov, perc_ident) = (float(x) for x in (aln_cov, perc_ident))
            ref_start -= 1# refstart: 0-based, refend: exclusive
            yield TilingContig._make([
                ref_name, ref_start, ref_end, gap2next, contig_len,
                aln_cov, perc_ident, ori, name])


def merge_contigs_and_ref(contig_seqs, ref_seq, tiling_file, out_fh=sys.stdout):
    """Merged contigs and reference based tiling data
    """

    out_fh.write(">joined\n")
    last_refend = 0# exclusive
    contig = None
    for contig in parse_tiling(tiling_file):
        assert contig.ref_name in ref_seq, (
            "Tiling reference name '{}' not found in refseqs".format(
                contig.ref_name))

        # if there was a gap before this contig, fill with reference
        if contig.ref_start > last_refend:
            LOG.debug("ref {}+1:{}".format(last_refend, contig.ref_start))
            sq = ref_seq[contig.ref_name][last_refend:contig.ref_start]
            out_fh.write(sq)

        # if there's overlap with the next contig we clip the current
        # one (assumes all contigs are equally good)
        printto = None
        if contig.gap2next < 0:
            printto = contig.gap2next
        if contig.ori == '+':
            LOG.debug("con+ {}+1:{}".format(0, printto))
            sq = contig_seqs[contig.name][:printto]
        elif contig.ori == '-':
            LOG.debug("con- {}+1:{}".format(0, printto))
            sq = rev_comp(contig_seqs[contig.name])[:printto]
        else:
            raise ValueError(contig.ori)
        out_fh.write(sq)
        last_refend = contig.ref_end

    if contig is not None:
        if last_refend < len(ref_seq[contig.ref_name]):
            LOG.debug("ref {}+1:".format(last_refend))
            sq = ref_seq[contig.ref_name][last_refend:]
            out_fh.write(sq)


def main():
    """
    The main function
    """

    parser = cmdline_parser()
    args = parser.parse_args()

    logging_level = logging.WARN + 10*args.quiet - 10*args.verbose
    LOG.setLevel(logging_level)

    for fname in [args.fref, args.fcontigs]:
        if not os.path.exists(fname):
            LOG.fatal("file '{}' does not exist.".format(fname))
            sys.exit(1)
    for fname in [args.fout]:
        if os.path.exists(fname):
            LOG.fatal("Refusing to overwrite existing file {}'.".format(fname))
            sys.exit(1)


    tmp_files = []

    # run mummer's nucmer
    #
    out_prefix = tempfile.NamedTemporaryFile(
        delete=False, dir=args.tmp_dir).name
    fdelta = run_nucmer(args.fref, args.fcontigs, out_prefix)
    tmp_files.append(fdelta)
    LOG.info("Delta written to '{}'".format(fdelta))

    fpseudo, ftiling = run_showtiling(fdelta)
    tmp_files.extend([fpseudo, ftiling])

    if args.dont_fill_with_ref:
        LOG.info("Not replacing gaps with ref. Copying to '{}'".format(args.fout))
        shutil.copyfile(fpseudo, args.fout)
        if not args.keep_temp_files:
            for f in tmp_files:
                os.unlink(f)
        else:
            LOG.info("Not deleting temp files")
        sys.exit(0)


    # load reference and contigs
    #
    ref_seq = dict((x[0].split()[0], x[1])
                   for x in fasta_iter(args.fref))
    assert len(ref_seq) == 1, ("Only one reference sequence supported for N filling")
    contigs = dict((x[0].split()[0], x[1])
                   for x in fasta_iter(args.fcontigs))

    merge_contigs_and_ref(contigs, ref_seq, ftiling, out_fh=sys.stdout)




if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
