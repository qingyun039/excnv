from os.path import dirname, join
from setuptools import setup

setup_args = {}

DIR = (dirname(__file__) or '.')
with open(join(DIR, 'cnvsq', '_version.py')) as handle:
    VERSION = handle.readline().split('=')[-1].strip().replace('"', '')

install_requires = [
    'CNVkit >= 0.9.9',
    'numpy >= 1.22.3',
    'pandas >= 1.4.2',
    'pybedtools >= 0.9.0',
    'pyBigWig >= 0.3.18',
    'pyfaidx >= 0.6.4',
    'pysam >= 0.19.0',
]
    

setup_args.update(
        name = 'cnvsq',
        version = VERSION,
        description = __doc__,
        author = 'sl chen',
        author_email = 'no email',
        url = 'no url',
        packages = ['cnvsq', 'cnvsq.classifycnv'],
        include_package_data=True,
        entry_points = {'console_scripts': ['dacnvseq = cnvsq.commands:main'], },
        install_requires = install_requires,
        classifiers = [
            "Development Status :: 4 - Beta",
            "Environment :: Console",
            "Intended Audience :: Science/Research",
            "Intended Audience :: Healthcare Industry",
            "Licence :: OSI Approved :: Apache Software Licence",
            "Operating System :: POSIX",
            "Operating System :: POSIX :: Linux",
            "Programming Language :: Python",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Scientific/Engineering :: Medical Science Apps"
        ]
)

setup(**setup_args)
