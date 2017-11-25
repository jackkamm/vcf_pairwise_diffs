import os
import ast
import pprint
import vcf_pairwise_diffs

# test.vcf was copied from htslib/test/test-vcf-api.out
pairwise_diffs = vcf_pairwise_diffs.pairwise_diffs("test.vcf")
for k in ('n_comparisons', 'n_differences'):
    pairwise_diffs[k] = pairwise_diffs[k].tolist()

pprint.pprint(pairwise_diffs)

with open("test.out") as f:
    assert ast.literal_eval(f.read()) == pairwise_diffs
