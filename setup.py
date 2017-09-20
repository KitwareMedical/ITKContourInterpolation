# -*- coding: utf-8 -*-
from __future__ import print_function
from os import sys

try:
    from skbuild import setup
except ImportError:
    print('scikit-build is required to build from source.', file=sys.stderr)
    print('Please run:', file=sys.stderr)
    print('', file=sys.stderr)
    print('  python -m pip install scikit-build')
    sys.exit(1)

setup(
    name='itk-morphologicalcontourinterpolation',
    version='1.0.0',
    author='Dženan Zukić',
    author_email='community@itk.org',
    packages=['itk'],
    package_dir={'itk': 'itk'},
    download_url=r'https://github.com/KitwareMedical/ITKMorphologicalContourInterpolation',
    description=r'Interslice interpolations of label images based on their morphology .',
    long_description='ITK is an open-source, cross-platform library that '
                     'provides developers with an extensive suite of software '
                     'tools for image analysis. This package implements morphological contour interpolation. '
                     'For more information, see  '
                     'Zukić D., Vicory J., McCormick M., Wisse L., Gerig G., Yushkevich P., Aylward S. '
                     '"ND Morphological Contour Interpolation" '
                     'http://insight-journal.org/browse/publication/977 '
                     'http://hdl.handle.net/10380/3563 ',
    classifiers=[
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python",
        "Programming Language :: C++",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Medical Science Apps.",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Software Development :: Libraries",
        "Operating System :: Android",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Operating System :: Unix",
        "Operating System :: MacOS"
        ],
    license='Apache',
    keywords='ITK InsightToolkit',
    url=r'https://itk.org/',
    install_requires=[
        r'itk'
    ]
    )
