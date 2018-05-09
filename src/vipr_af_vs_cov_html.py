#!/usr/bin/env python3
"""Plotting of ViPR3 coverage and AF values
"""

#--- standard library imports
#
from argparse import ArgumentParser
from collections import OrderedDict
import gzip
from math import ceil, log10
from os.path import exists

#--- third-party imports
#
from plotly.graph_objs import Figure, Layout, Scatter
from plotly.offline import plot

#--- project specific imports
#
# /


__author__ = "Andreas WILM, LIEW Jun Xian"
__email__ = "wilma@gis.a-star.edu.sg, liewjx@gis.a-star.edu.sg"
__copyright__ = "2018 Genome Institute of Singapore"
__license__ = "The MIT License (MIT)"


def parse_genomecov(genomecov_gzfile):
    """Parse gzipped bedtools genomecov output"""
    genomecov = OrderedDict()
    with gzip.open(genomecov_gzfile) as fh:
        for line in fh:
            sq, pos, cov = line.decode().rstrip().split("\t")
            pos = int(pos) - 1
            # float for support of scientific notation
            cov = int(float(cov))
            if sq not in genomecov:
                genomecov[sq] = OrderedDict()
            genomecov[sq][pos] = cov
    return genomecov


def af_from_vcf(vcf_gz, snps_only=False):
    """Parse AF from INFO fields of LoFreq vcf's"""
    afs = OrderedDict()
    with gzip.open(vcf_gz) as fh:
        for line in fh:
            if line.decode().startswith("#"):
                continue
            sq, pos, _, ref, alt, _, _, info = line.decode().rstrip().split("\t")[:8]
            if snps_only:
                if len(ref) > 1 or len(alt) > 1:
                    continue
            pos = int(pos) - 1
            af = [float(x[4:]) for x in info.split(";") if x.startswith("AF=")][0]
            if sq not in afs:
                afs[sq] = OrderedDict()
            # keep max af if multi-allelic
            af2 = afs[sq].get(pos)
            if af2 and af2 > af:
                af = af2
            afs[sq][pos] = af
    return afs


def interactive_plot(cov_gzfile, vcf_gzfile, plot_file):
    """read cov and vcf and creates plot"""
    genomecov = parse_genomecov(cov_gzfile)
    afs = af_from_vcf(vcf_gzfile, snps_only=True)

    # only one genome allowed
    assert len(genomecov) == 1
    # only one sq is supported for plotting but could be 0 if vcf is empty
    assert len(afs) < 2
    sq = list(genomecov.keys())[0]

    if len(afs):
        assert sq in afs

    x1 = list(genomecov[sq].keys())
    y1 = list(genomecov[sq].values())

    if len(afs):
        x2 = list(afs[sq].keys())
        y2 = list(afs[sq].values())
    else:
        x2 = []
        y2 = []

    trace1 = Scatter(
        line=dict(color="blue", width=5),
        # marker=dict(size = 10, color="blue"),
        mode="lines", name="Coverage",
        x=x1, y=y1
    )

    trace2 = Scatter(
        marker=dict(size=10, color="red"),
        mode="markers", name="AF",
        x=x2, y=y2,
        yaxis="y2"
    )

    data = [trace1, trace2]

    layout = Layout(
        font=dict(color="black", family="Arial, monospace", size=20),
        #        legend=dict(x = 1.05, y = 1.05),
        title=sq,
        xaxis=dict(
            autorange=True,
            dtick=2000,
            exponentformat="none",
            showexponent="none",
            showgrid=True,
            showline=True,
            title="Position",
            titlefont=dict(
                color="black"
            )
        ),
        yaxis=dict(
            dtick=1,
            exponentformat="power",
            range=[0, ceil(log10(max(y1)))],
            showexponent="all",
            showgrid=False,
            showline=False,
            tick0=1,
            title="Coverage",
            titlefont=dict(
                color="blue"
            ),
            type="log",
        ),
        yaxis2=dict(
            dtick=0.2,
            exponentformat="none",
            overlaying="y",
            range=[0, 1.0],
            showexponent="none",
            showgrid=False,
            showline=True,
            side="right",
            tick0=0,
            tickfont=dict(
                color="black"
            ),
            title="AF",
            titlefont=dict(
                color="red"
            )
        )
    )

#    plot(Figure(data = data, layout = layout), auto_open = False, filename = plot_file + ".html", show_link = False, image = "png", image_filename = plot_file, image_height = 1080, image_width = 1920)
    plot(Figure(data=data, layout=layout), auto_open=False,
         filename=plot_file, show_link=False)


def main():
    """main function"""
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--vcf", required=True,
                        help="LoFreq vcf files (gzipped)")
    parser.add_argument("-c", "--cov", required=True,
                        help="bedtools genomecov file (gzipped)")
    parser.add_argument("-p", "--plot", required=True,
                        help="html plot filename")
    args = parser.parse_args()

    assert not exists(args.plot)
    interactive_plot(args.cov, args.vcf, args.plot)


if __name__ == "__main__":
    main()
