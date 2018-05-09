# https://packaging.python.org/tutorials/distributing-packages/
import sys
from distutils.core import setup

python_requires='>=3',

PACKAGE_NAME="ViPR-Tools"
PACKAGE_TARNAME="vipr-tools"
PACKAGE_VERSION="0.1"
PACKAGE_BUGREPORT="wilma@gis.a-star.edu.sg"

setup(name = PACKAGE_NAME,
      version = PACKAGE_VERSION,
      description="Tools for the viral amplicon analytis pipeline ViPR",
      author="Andreas Wilm",
      author_email=PACKAGE_BUGREPORT,
      requires = ['plotly (>=2.5.1)'],
      scripts = [
          "src/polish_viral_ref.sh",
          "src/simple_contig_joiner.py",
          "src/vipr_af_vs_cov_html.py",
          "src/vipr_af_vs_cov_pdf.py",
          "src/vipr_gaps_to_n.py",
      ],
      # http://pypi.python.org/pypi?%3Aaction=list_classifiers
      classifiers=['Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'Natural Language :: English',
                   'Operating System :: Unix',
                   'Programming Language :: Python :: 3',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   ],
      keywords='bioinformatics'
      )
