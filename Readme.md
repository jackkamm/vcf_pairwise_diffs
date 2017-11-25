This is just a simple module to compute the number of pairwise
differences between all samples in a VCF file.
It computes the number of pairwise differences in C
using htslib, then uses Cython and pysam to pass the result
back to Python.
