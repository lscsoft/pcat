#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Setup the PCAT package
"""

from __future__ import print_function

import glob
import os.path

from setuptools import (setup, find_packages)

# set basic metadata
PACKAGENAME = 'pcat'
DISTNAME = 'pcat'
AUTHOR = 'Daniele Trifiro'
AUTHOR_EMAIL = 'brethil@phy.olemiss.edu'
LICENSE = 'GPLv3'

# -- dependencies -------------------------------------------------------------

install_requires = [
]
requires = [
    'numpy',
    'scipy'
    'matplotlib',
    'glue',
    'gwpy',
    'pylal',
    'sklearn',
]

# -- run setup ----------------------------------------------------------------

packagenames = find_packages()
print(packagenames)
scripts = glob.glob(os.path.join('bin', '*'))
print(scripts)

setup(name=DISTNAME,
      provides=[PACKAGENAME],
      description=None,
      long_description=None,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      url='https://github.com/lscsoft/pcat',
      packages=packagenames,
      scripts=scripts,
      install_requires=install_requires,
      requires=requires,
      dependency_links=[
          'http://software.ligo.org/lscsoft/source/glue-1.49.1.tar.gz'
              '#egg=glue-1.49.1',
      ],
      use_2to3=True,
      classifiers=[
          'Programming Language :: Python',
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Developers',
          'Natural Language :: English',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Scientific/Engineering :: Physics',
          'Operating System :: POSIX',
          'Operating System :: Unix',
          'Operating System :: MacOS',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
      ],
      )
