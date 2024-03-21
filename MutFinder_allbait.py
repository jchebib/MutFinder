##############################################################################################################
# Python program: Finds mutations in focal samples using all others as bait (impurity filters) 
# Author: Jobran Chebib, March 2024 (based on August 2022; March 2022; based on October 2019 version)
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
parser.add_argument("--max_impurity", dest="max_impurity", default=42, type=int, help="The max number of ALT alleles in combined non-focal genotypes")
parser.add_argument("--allele_balance_ratio", dest="allele_balance_ratio", default=0.25, type=float, help="The min/max ratio of ALT to REF alleles in a heterozygote")
parser.add_argument("--chromosome", dest="chromosome", default=1, type=str, help="The chromosome of the desired sequences")
parser.add_argument("--start", dest="start", default=191000000, type=int, help="The start position in the range of desired sequences")
parser.add_argument("--finish", dest="finish", default=191002000, type=int, help="The end position in the range of desired sequences")

# get data and argument values
args = parser.parse_args()
vcf_data = VCF(args.vcf_file)
#output_file = args.output_file
min_DP = args.min_DP
avg_DP = args.avg_DP
min_QUAL = args.min_QUAL
min_sample_coverage = args.min_sample_coverage
max_impurity = args.max_impurity
allele_balance_ratio = args.allele_balance_ratio
start = args.start
finish = args.finish
chromosome = args.chromosome
#vcf_data = VCF('GenotypeGVCFs_gen8_norep_chr2.vcf.gz')
#vcf_data = VCF('GenotypeGVCFs_gen8_2_10-11M.vcf.gz')
#vcf_data = VCF('GenotypeGVCFs_8refSorted_1_0-20M.vcf.gz')
#vcf_data = VCF('GenotypeGVCFs_gen8_2_50Kz.vcf.gz')

# initialize variables
nbr_passed_sites=0
f_nbr_passed_sites=0 # passed, but founders excluded
nbr_total_sites=0
passed_sites=[]
f_passed_sites=[] # passed, but founders excluded
sites=[]
nbr_phased_candidate_sites=0
nbr_homo_candidate_sites=0
nbr_hetero_candidate_sites=0
candidate_sites=[]
hetero_candidate_sites=[]
homo_candidate_sites=[]
phased_candidate_sites=[]
unbalanced_candidate_sites=[]
f_hetero_candidate_sites=[] # variant in descendents and at least one founder
f_homo_candidate_sites=[] # variant in descendents and at least one founder
f_phased_candidate_sites=[] # variant in descendents and at least one founder
f_unbalanced_candidate_sites=[] # variant in descendents and at least one founder

#f = open(output_file,"w+")
##f = open("founder_vcf_data.txt","w+")
num_founders = 10   # e.g. number of founders in MA mouse lines
num_samples = len(vcf_data.samples)   # e.g. number of samples in MA mouse line (including founders)

# show sample names
print(vcf_data.samples)
print("Number of samples: ", num_samples)

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
#       if(sum(v.gt_alt_depths) > max_impurity or sum(v.gt_alt_depths) < 0): continue
#	if(not(v.is_snp)): continue   # only variant SNPs
	sites.append(v.start+1)
#	f.write(str(v.CHROM) + ':' + str(v.start+1) + ':' + str(v.gt_types.tolist()) + ':' + str(v.gt_depths.tolist()) + ':' + str(v.gt_alt_depths.tolist()) + ':' + str(v.gt_phases.tolist()) + ':' + str(v.is_snp) + ':' + str(v.is_indel) + '\n')
	for p in range(num_samples):
		if(v.gt_types[p]!=2): # check focal sample not unknown
			if(v.gt_types[p]!=0 and np.all(np.delete(v.gt_types, p)==0)):  # check that only focal sample is not Homo REF
				if((sum(np.delete(v.gt_alt_depths, p)) <= max_impurity)): # check non-focal samples combined have fewer ALT reads than max_impurity
					passed_sites.append(v.start+1)
#        				print(str(v.CHROM) + ':' + str(v.start+1) + ':' + str(v.gt_types.tolist()) + ':' + str(v.gt_depths.tolist()) + ':' + str(v.gt_alt_depths.tolist()) + ':' + str(v.gt_phases.tolist()) + ':' + str(v.is_snp) + ':' + str(v.is_indel) + '\n')
					if(v.gt_types[p]==1): # if focal sample is HET
						#if(v.gt_phases[p]==1):
						if(np.any(v.gt_phases == 1)):
							phased_candidate_sites.append(v.start+1)
							print("Phased hetero candidate site for sample", vcf_data.samples[p], v.CHROM, v.start+1, "SNP: ", str(v.is_snp), v.gt_types, v.gt_alt_depths, v.gt_depths, v.gt_phases, v.QUAL, "INDEL: ", str(v.is_indel))
						elif(v.gt_alt_depths[p]/float(v.gt_depths[p]) <= allele_balance_ratio or v.gt_alt_depths[p]/float(v.gt_depths[p])>=(1-allele_balance_ratio)):
							unbalanced_candidate_sites.append(v.start+1)
							print("Unbalanced hetero candidate site for sample", vcf_data.samples[p], v.CHROM, v.start+1, "SNP: ", str(v.is_snp), v.gt_types, v.gt_alt_depths, v.gt_depths, v.gt_phases, v.QUAL, "INDEL: ", str(v.is_indel))
						#elif(v.gt_phases[p]==0):
						elif(np.all(v.gt_phases == 0)):
							hetero_candidate_sites.append(v.start+1)
							print("Hetero candidate site for sample", vcf_data.samples[p], v.CHROM, v.start+1, "SNP: ", str(v.is_snp), v.gt_types, v.gt_alt_depths, v.gt_depths, v.gt_phases, v.QUAL, "INDEL: ", str(v.is_indel))
					if(v.gt_types[p]==3): # if focal sample is HOMO ALT
						if(np.any(v.gt_phases == 1)):
							phased_candidate_sites.append(v.start+1)
							print("Phased homo ALT candidate site for sample", vcf_data.samples[p], v.CHROM, v.start+1, "SNP: ", str(v.is_snp), v.gt_types, v.gt_alt_depths, v.gt_depths, v.gt_phases, v.QUAL, "INDEL: ", str(v.is_indel))
						#elif(v.gt_phases[p]==0):
						elif(np.all(v.gt_phases == 0)):
							homo_candidate_sites.append(v.start+1)
							print("Homo ALT candidate site for sample", vcf_data.samples[p], v.CHROM, v.start+1, "SNP: ", str(v.is_snp), v.gt_types, v.gt_alt_depths, v.gt_depths, v.gt_phases, v.QUAL, "INDEL: ", str(v.is_indel))

####################### perform same filters with founders exempt to capture variants in both founders AND descendents 
			if(p>=num_founders):  # not founders
				if(v.gt_types[p and (0 or 1 or 2 or 3 or 4 or 5 or 6 or 7 or 8 or 9)]!=0 and np.all(np.delete(v.gt_types, [0,1,2,3,4,5,6,7,8,9,p])==0)):  # check only focal sample and either founder are not Homo REF
					if((sum(np.delete(v.gt_alt_depths, [0,1,2,3,4,5,6,7,8,9,p])) <= max_impurity)): # check non-focal samples combined have fewer ALT reads than max_impurity
						f_passed_sites.append(v.start+1)

#        			print(str(v.CHROM) + ':' + str(v.start+1) + ':' + str(v.gt_types.tolist()) + ':' + str(v.gt_depths.tolist()) + ':' + str(v.gt_alt_depths.tolist()) + ':' + str(v.gt_phases.tolist()) + ':' + str(v.is_snp) + ':' + str(v.is_indel) + '\n')
						if(v.gt_types[p]==1): # if focal sample is HET
							if(np.any(v.gt_phases == 1)):
								f_phased_candidate_sites.append(v.start+1)
								print("FOUNDER: Phased hetero candidate site for sample", vcf_data.samples[p], v.CHROM, v.start+1, "SNP: ", str(v.is_snp), v.gt_types, v.gt_alt_depths, v.gt_depths, v.gt_phases, v.QUAL, "INDEL: ", str(v.is_indel))
							elif(v.gt_alt_depths[p]/float(v.gt_depths[p]) <= allele_balance_ratio or v.gt_alt_depths[p]/float(v.gt_depths[p])>=(1-allele_balance_ratio)):
								f_unbalanced_candidate_sites.append(v.start+1)
								print("FOUNDER: Unbalanced hetero candidate site for sample", vcf_data.samples[p], v.CHROM, v.start+1, "SNP: ", str(v.is_snp), v.gt_types, v.gt_alt_depths, v.gt_depths, v.gt_phases, v.QUAL, "INDEL: ", str(v.is_indel))
							#elif(v.gt_phases[p]==0):
							elif(np.all(v.gt_phases == 0)):
								f_hetero_candidate_sites.append(v.start+1)
								print("FOUNDER: Hetero candidate site for sample", vcf_data.samples[p], v.CHROM, v.start+1, "SNP: ", str(v.is_snp), v.gt_types, v.gt_alt_depths, v.gt_depths, v.gt_phases, v.QUAL, "INDEL: ", str(v.is_indel))
						if(v.gt_types[p]==3): # if focal sample is HOMO ALT
							if(np.any(v.gt_phases == 1)):
								f_phased_candidate_sites.append(v.start+1)
								print("FOUNDER: Phased homo ALT candidate site for sample", vcf_data.samples[p], v.CHROM, v.start+1, "SNP: ", str(v.is_snp), v.gt_types, v.gt_alt_depths, v.gt_depths, v.gt_phases, v.QUAL, "INDEL: ", str(v.is_indel))
							elif(np.all(v.gt_phases == 0)):
								f_homo_candidate_sites.append(v.start+1)
								print("FOUNDER: Homo ALT candidate site for sample", vcf_data.samples[p], v.CHROM, v.start+1, "SNP: ", str(v.is_snp), v.gt_types, v.gt_alt_depths, v.gt_depths, v.gt_phases, v.QUAL, "INDEL: ", str(v.is_indel))

###### OUTPUT ##################################################				
unique_passed_sites = np.unique(passed_sites)
nbr_passed_sites = len(unique_passed_sites)
f_unique_passed_sites = np.unique(f_passed_sites)
f_nbr_passed_sites = len(f_unique_passed_sites)

#	f.write(str(v.CHROM) + ':' + str(v.start+1) + ':' + str(v.gt_types.tolist()) + ':' + str(v.gt_depths.tolist()) + ':' + str(v.gt_alt_depths.tolist()) + ':' + str(v.gt_phases.tolist()) + ':' + str(v.is_snp) + ':' + str(v.is_indel) + '\n')
#f.close()

print("List passed sites: ", unique_passed_sites)
print("List passed sites (founders excluded): ", f_unique_passed_sites)
print("Total sites: ", nbr_total_sites, "QUAL/DEPTH pass sites: ", len(np.unique(sites)), "Passed sites: ", nbr_passed_sites)
print("QUAL >", min_QUAL, min_sample_coverage, "proportion have depth >", min_DP, "Average genome depth:", avg_DP, "Max impurity:", max_impurity, ". Allele balance ratio: ", allele_balance_ratio)
print("Passed sites in descendents only (founders excluded from filters): ", f_nbr_passed_sites)


