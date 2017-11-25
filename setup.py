import os
import numpy
from setuptools import setup, Extension
import pysam

wrapper_mod = Extension("vcf_pairwise_diffs.wrapper",
                        sources=["vcf_pairwise_diffs/wrapper.pyx",
                                 "src/vcfpairdiffs.c"],
                        include_dirs=pysam.get_include() + [
                            "src", numpy.get_include()],
                        extra_link_args=pysam.get_libraries())

setup(name="vcf_pairwise_diffs",
	  packages=["vcf_pairwise_diffs"],
      ext_modules=[wrapper_mod],
      install_requires=['numpy', 'pysam'],
      setup_requires=['numpy', 'pysam', 'cython'])
