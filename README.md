# MutFinder
# Python program: Finds mutations in focal samples using all others as bait (impurity filters)
# Author: Jobran Chebib, March 2024 
# Build for: Python 2.7.5
# Requires: Blas        $ sudo yum install blas
# Requires: Cython      $ sudo yum install Cython
# Requires: cyvcf2      $ pip install cyvcf2
# Requires: NumPy       $ pip install numpy
# get dependencies
import argparse,sys
from cyvcf2 import VCF
import numpy as np
np.set_printoptions(threshold=np.inf)

# get argument values
parser = argparse.ArgumentParser(description="ID mutations present in a VCF of a single starting genoytpe of MA", usage="MutFinder_cyvcf2_mus_bait.py [options] input_file.vcf(.gz)")
parser.add_argument("vcf_file",help="The name of the vcf file")
parser.add_argument("--min_QUAL", dest="min_QUAL", default=90, type=int, help="The min QUAL for a site to be considered as a candidate mutational site")
parser.add_argument("--min_DP", dest="min_DP", default=10, type=int, help="The min coverage of depth for a site to be considered as a candidate mutational site")
parser.add_argument("--avg_DP", dest="avg_DP", default=30, type=int, help="The average coverage of depth, for which a site is excluded as a candidate mutational site if the site coverage is greater than 2 times the average depth")
parser.add_argument("--min_sample_coverage", dest="min_sample_coverage", default=1.0, type=float, help="The min proportion of samples that must meet the coverage of depth at a site")
parser.add_argument("--max_impurity", dest="max_impurity", default=42, type=int, help="The max number of ALT alleles in combined non-focal genotypes")
parser.add_argument("--allele_balance_ratio", dest="allele_balance_ratio", default=0.25, type=float, help="The min/max ratio of ALT to REF alleles in a heterozygote")
parser.add_argument("--chromosome", dest="chromosome", default=1, type=str, help="The chromosome of the desired sequences")
parser.add_argument("--start", dest="start", default=191000000, type=int, help="The start position in the range of desired sequences")
parser.add_argument("--finish", dest="finish", default=191002000, type=int, help="The end position in the range of desired sequences")

**EXAMPLE USAGE:**
> python MutFinder_allbait.py --chromosome NC_000067.7 --start 20000001 --finish 40000000 --max_impurity 49 --allele_balance_ratio 0.25 input.vcf.gz > output.txt
