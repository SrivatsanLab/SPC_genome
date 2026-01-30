from setuptools import setup, find_packages

setup(
    name="mutation_profiling",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'matplotlib',
        'pyfaidx',
        'cyvcf2',
        'tqdm',
        'scipy',
        'biopython',
    ],
    author="Dustin Mullane",
    author_email="",
    description="A library for mutation spectrum analysis in bulk and single-cell sequencing data",
    long_description=open('mutation_profiling/README.md').read(),
    long_description_content_type="text/markdown",
    python_requires='>=3.7',
)
