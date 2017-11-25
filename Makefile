.PHONY: build clean

build:
	python setup.py build_ext --inplace

clean:
	rm -rf build
	rm vcf_pairwise_diffs/*c
	rm vcf_pairwise_diffs/*so
