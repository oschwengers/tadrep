
from os import path
from setuptools import setup, find_packages
import tadrep


# Get the long description from the README file
setup_dir = path.abspath(path.dirname(__file__))
with open(path.join(setup_dir, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='tadrep',
    version=tadrep.__version__,
    description='TaDRep: Targeted Detection and Reconstruction of Plasmids',
    keywords=['bioinformatics', 'bacteria', 'plasmids'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    license='GPLv3',
    author='Oliver Schwengers',
    author_email='oliver.schwengers@computational.bio.uni-giessen.de',
    url='',
    packages=find_packages(include=['tadrep', 'tadrep.*', 'database', 'database.*']),
    python_requires='>=3.6',
    include_package_data=False,
    zip_safe=False,
    install_requires=[
        'biopython >= 1.72',
        'xopen >= 1.1.0'
    ],
    entry_points={
        'console_scripts': [
            'tadrep=tadrep.main:main',
            'tadrep_db=database.main:main'
        ]
    },
    classifiers=[
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Development Status :: 4 - Beta ',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Natural Language :: English'
    ],
    project_urls={
        'Bug Reports': '',
        'Source': ''
    },
)
