##############################################################################################################
# Python program: Finds number of callable sites that pass QUAL and DP filters 
# Author: Jobran Chebib, December 2022 (companion of MutFinder_nobait.py August 2022)
# Build for: Python 2.7.5
# Requires: Blas        $ sudo yum install blas
# Requires: Cython      $ sudo yum install Cython
# Requires: cyvcf2      $ pip install cyvcf2
# Requires: NumPy       $ pip install numpy
##############################################################################################################
# vcf fields:   vcf.samples, vcf.seqnames
# vcf row fields: v.REF, v.ALT (REF='A', ALT=['C', 'T']), v.CHROM, v.start,
#               v.end, v.ID, v.FILTER, v.QUAL, v.genotypes,
#               v.gt_types (array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT), v.gt_depths,
#               v.gt_ref_depths, v.gt_alt_depths, v.gt_phases, v.gt_quals,
#               v.gt_bases, v.call_rate, v.INFO.get('MQ'), v.INFO.get('DP')
#               v.is_deletion, v.is_indel, v.is_snp, v.is_sv, v.is_transition
##############################################################################################################
# get dependencies
import argparse,sys
from cyvcf2 import VCF
import numpy as np
np.set_printoptions(threshold=np.inf)

# get argument values
parser = argparse.ArgumentParser(description="ID mutations present in a VCF of a single starting genoytpe of MA", usage="MutFinder_cyvcf2_mus_bait.py [options] input_file.vcf(.gz)")
#parser = argparse.ArgumentParser(description="ID mutations present in a VCF of a single starting genoytpe of MA", usage="MutFinder_cyvcf2_mus_bait.py [options] output_file.txt")
parser.add_argument("vcf_file",help="The name of the vcf file")
#parser.add_argument("output_file",help="The name of the output file")
parser.add_argument("--min_QUAL", dest="min_QUAL", default=90, type=int, help="The min QUAL for a site to be considered as a candidate mutational site")
parser.add_argument("--min_DP", dest="min_DP", default=10, type=int, help="The min coverage of depth for a site to be considered as a candidate mutational site")
parser.add_argument("--avg_DP", dest="avg_DP", default=30, type=int, help="The average coverage of depth, for which a site is excluded as a candidate mutational site if the site coverage is greater than 2 times the average depth")
parser.add_argument("--min_sample_coverage", dest="min_sample_coverage", default=1.0, type=float, help="The min proportion of samples that must meet the coverage of depth at a site")
parser.add_argument("--chromosome", dest="chromosome", default=2, type=str, help="The chromosome of the desired sequences")
parser.add_argument("--start", dest="start", default=10179965, type=int, help="The start position in the range of desired sequences")
parser.add_argument("--finish", dest="finish", default=10179966, type=int, help="The end position in the range of desired sequences")

# get data and argument values
args = parser.parse_args()
vcf_data = VCF(args.vcf_file)
#output_file = args.output_file
min_DP = args.min_DP
avg_DP = args.avg_DP
min_QUAL = args.min_QUAL
min_sample_coverage = args.min_sample_coverage
start = args.start
finish = args.finish
chromosome = args.chromosome
#vcf_data = VCF('GenotypeGVCFs_gen8_norep_chr2.vcf.gz')
#vcf_data = VCF('GenotypeGVCFs_gen8_2_10-11M.vcf.gz')
#vcf_data = VCF('GenotypeGVCFs_8refSorted_1_0-20M.vcf.gz')
#vcf_data = VCF('GenotypeGVCFs_gen8_2_50Kz.vcf.gz')

# initialize variables
f_nbr_passed_sites=0 # passed, but founders excluded
nbr_total_sites=0
sites=[]

#f = open(output_file,"w+")
##f = open("founder_vcf_data.txt","w+")
num_founders = 2   # e.g. number of founders in MA mouse lines
num_samples = len(vcf_data.samples)   # e.g. number of samples in MA mouse line (including founders)

# show sample names
#print(vcf_data.samples)
#print("Number of samples: ", num_samples)

##### iterate through founder vcf to filter sites #####
seq_range = str(chromosome) + ':' + str(start) + '-' + str(finish)
for v in vcf_data(seq_range):
#for v in vcf_data('4:95999990-95999999'):
#for v in vcf_data:
	nbr_total_sites += 1
	if (v.QUAL is not None and v.QUAL < min_QUAL): continue
	if ((np.sum(v.gt_depths >= min_DP)/float(len(v.gt_depths)) < min_sample_coverage)): continue # check certain proportion of founders have depth > min_DP
	if ((np.sum(v.gt_depths)/float(len(v.gt_depths)) > (2*avg_DP))): continue # check average depth < 2 times the average genome coverage
#	if(np.any(v.gt_depths < min_DP)): continue
#	if(np.any(v.gt_depths > (2*avg_DP))): continue
#	if(not(v.is_snp)): continue   # only variant SNPs
#	sites.append(v.start+1)
	print(str(v.CHROM) + '\t' + str(v.start+1))

#print("Total sites: ", nbr_total_sites, "QUAL/DEPTH pass sites: ", len(np.unique(sites)))
#print("QUAL >", min_QUAL, min_sample_coverage, "proportion have depth >", min_DP, "Average genome depth:", avg_DP)
