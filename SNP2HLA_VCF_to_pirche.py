import sys
import pandas as pd
import gzip
from collections import defaultdict

vcf_filename = "./SNP2HLA_Imputation/Penn.PTI/chr6.dose.vcf.gz"
# vcf_df = pandas.read_csv(vcf_filename,delimiter='\t',compression="gzip")

sample_IDs = []
loci_ID = ['A', 'C', 'B', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1']
sample_ID_HLA_geno = defaultdict(lambda: defaultdict(dict)) # store the HLA allele-level genotype for all locus
with gzip.open(vcf_filename, 'rt') as file:
    line_index = 0
    for line in file:

        line = str(line)
        if line.startswith('##'):
            continue  # Skip header lines starting with '#'

        if line.startswith('#'):
            fields = line.strip().split('\t')
            print ("Number of columns: " + str(len(fields)))
            # print (line)
            for x in range(9, len(fields)):   # 896 for PMBB, 1163 for UPenn
                print (x)
                print (fields[x])
                sample_IDs.append(str(fields[x]))
            continue

        # print (sample_IDs)

        # print (line)
        # print ("Number of rows: " + str(len(line)))  # 40K characters per line
        fields = line.strip().split('\t')

        # CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	05-01308_05-01308
        # 6	29910248	HLA_A*01:01	A	T	.	PASS	AF=0.10860;MAF=0.10860;R2=0.98821;IMPUTED	GT:DS:HDS:GP	0|0:0.006:0.000,0.006:0.994,0.006,0.000	0|0:0:0,0:1,0,0	0|0:0:0,0:1,0,0

        # 1163 columns
        # print ("Number of columns: " + str(len(fields)))

        chrom = fields[0]  # always chromosome 6
        pos = int(fields[1])  # physical map for SNPs - might be start of the HLA-A gene?
        id = fields[2]  # allele name - parse this 

        # nonHLA : rs9277489, SNPS_DPB1_9913_33053731_intron4, SNPS_DPB1_9944_33053762_intron4_A

        # select rows for HLA specificities and reformat allele name
        hla_specificity = ""
        if (id.startswith("HLA")):

            (locus,hla_allele) = id.split("*")
            if (":" in hla_allele):
                (hla_prefix,hla_locus) = locus.split("_")
                hla_specificity = hla_locus + "*" + hla_allele
                
            else:  # skip first-field HLA variants
                continue

        else: # skip non-HLA variants
            continue

        # print (id)


        # if (hla_locus != "A"):
        #    continue
        # print (hla_locus)
        # print (hla_specificity)

        ref = fields[3]  # nucleotide for major allele e.g. A
        alt = fields[4]  # nucleotide for minor allele e.g. T
        qual = fields[5] # quality score - always just a "."

        # print (qual)

        filter = fields[6] # filter - always a pass for HLA

        # print (filter)

        info = fields[7] # AF=0.10860;MAF=0.10860;R2=0.98821;IMPUTED

        # skip info and derive the allele frequency later
        # every HLA variant is a minor allele (MAF = AF)

        format = fields[8]  # always "GT:DS:HDS:GP"

        # print (format)

        # genotype is indexes 9 through 1162 - sample IDs are in the header row
        genotype_dict = defaultdict(lambda: defaultdict(dict))
        genotype_GT_dict = defaultdict(lambda: defaultdict(dict))
        for i in range(9, len(fields)):
            sample_ID_index = i - 9
            sample_ID = sample_IDs[sample_ID_index]
            # print (sample_ID)
            genotype = fields[i]
            (geno_GT,geno_DS,geno_HDS,geno_GP) = genotype.split(',')
            # TBD - determine what the genotype string definition is - look at SNP2HLA documentation
            (allele_presence_geno,float_unknown_1,float_unknown2) = geno_GT.split(":")
            (allele1_presence,allele2_presence) = allele_presence_geno.split("|")
            for loci in loci_ID:
                if hla_locus != loci:
                    continue
                # print ("This is the loci:", loci)
                if allele_presence_geno != "0|0":  # only look at genotypes where the allele was present - either position can be positive
                    print("Sample ID: " + sample_ID)
                    print("GT_1: " + allele1_presence)
                    if allele1_presence == "1":
                        print(hla_locus)
                        if hla_locus == loci:
                            if loci not in sample_ID_HLA_geno[sample_ID]:  # put in allele if no previous allele
                                sample_ID_HLA_geno[sample_ID][loci] = hla_specificity
                            else:  # else make into a genotype
                                sample_ID_HLA_geno[sample_ID][loci] = str(sample_ID_HLA_geno[sample_ID][loci]) + "+" + hla_specificity
                            print("A locus type present 1")
                            # print(sample_ID_HLA_geno)
                    print("GT_2: " + allele2_presence)
                    if allele2_presence == "1":
                        if hla_locus == loci:
                            if loci not in sample_ID_HLA_geno[sample_ID]:  # put in allele if no previous allele
                                sample_ID_HLA_geno[sample_ID][loci] = hla_specificity
                            else:  # else make into a genotype
                                sample_ID_HLA_geno[sample_ID][loci] = str(sample_ID_HLA_geno[sample_ID][loci]) + "+" + hla_specificity
                            print("A locus type present 2")
                            # print(sample_ID_HLA_geno)
                    print(sample_ID_HLA_geno[sample_ID])
                    # print ("GT: " + geno_GT)    # 1|0:1.000:1.000    0|1:0.949:0.000
                    # print ("DS: " + geno_DS)
            genotype_dict[sample_ID][loci] = genotype
            genotype_GT_dict[sample_ID][loci] = geno_GT

            # print (genotype)

        # print (genotype_dict)

        # ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        # ##FORMAT=<ID=DS,Number=1,Type=Float,Description="Estimated Alternate Allele Dosage : [P(0/1)+2*P(1/1)]">
        # ##FORMAT=<ID=HDS,Number=2,Type=Float,Description="Estimated Haploid Alternate Allele Dosage ">
        # ##FORMAT=<ID=GP,Number=3,Type=Float,Description="Estimated Posterior Probabilities for Genotypes 0/0, 0/1 and 1/1 ">

        # GT stands for "Genotype" and represents the genotype of the sample at the variant locus. The genotype is typically represented as a combination of alleles, 
        #   where each allele is represented by a number or a character. For diploid organisms like humans, the genotype is typically represented by two alleles 
        #   separated by a delimiter (e.g., / or |). Examples of genotypes include 0/0 (homozygous reference), 0/1 (heterozygous), or 1/1 (homozygous alternate).
        # DS stands for "Genotype Dosage" and represents the estimated number of copies of the alternate allele in the sample's genotype. It provides a continuous value ranging from 0 to 2, 
        #   where 0 indicates no copies of the alternate allele, 1 indicates one copy, and 2 indicates two copies.
        # HDS stands for "Hard Depth Scaled" and represents the allele depth for each allele in the genotype. It provides the number of supporting reads for each allele in the genotype call. 
        #   This information is useful in assessing the confidence of the genotype call and can be used for filtering or quality control purposes.
        # GP stands for "Genotype Probability" and represents the genotype likelihoods in Phred-scaled format. It provides the probabilities of different possible genotypes given the observed data. 
        #   The values are typically expressed as negative log10 likelihoods, with higher values indicating higher likelihoods. The genotype probabilities are often used for genotype calling 
        #   and variant quality assessment.




        # print ("Line:" + line)
print (sample_ID_HLA_geno)  # should have a HLA genotype for each ID for locus A

# Make into a dataframe to make it easier to use
PIRCHE = pd.DataFrame.from_dict(sample_ID_HLA_geno, orient='index')

# the ones with n=1 have missing happairs
PIRCHE[["A_1", "A_2"]] = PIRCHE.A.str.split('+', expand=True)
PIRCHE[["C_1", "C_2"]] = PIRCHE.C.str.split('+', n=1, expand=True)
PIRCHE[["B_1", "B_2"]] = PIRCHE.B.str.split('+', n=1, expand=True)
PIRCHE[["DRB1_1", "DRB1_2"]] = PIRCHE.DRB1.str.split('+', n=1, expand=True)
PIRCHE[["DQA1_1", "DQA1_2"]] = PIRCHE.DQA1.str.split('+', n=1, expand=True)
PIRCHE[["DQB1_1", "DQB1_2"]] = PIRCHE.DQB1.str.split('+', expand=True)
PIRCHE[["DPA1_1", "DPA1_2"]] = PIRCHE.DPA1.str.split('+', n=1, expand=True)
PIRCHE[["DPB1_1", "DPB1_2"]] = PIRCHE.DPB1.str.split('+', expand=True)

PIRCHE = PIRCHE.drop(columns=['A', 'C', 'B', 'DRB1', 'DQA1', 'DQB1', 'DPA1', 'DPB1'])   # not needed anymore

# pd.set_option('display.max_columns', None)
# print(subjHLA.head())
# pirche_file = pd.read_csv(pirche_filename, sep=',', header=None)
# pairID_filename = './SNP2HLA_Imputation/Penn.PTI/Penn.PTI.DR.pairs.IDs.csv'
# pairID_file = pd.read_csv(pairID_filename, sep=',')

# CSV file
subjHLA_filename = 'Subj_HLA_' + "PennPTI" + ".csv"
PIRCHE.to_csv(subjHLA_filename, header=False)

# '05-04377_05-04377': 'A*34:02+A*68:02'
# '06-05177_06-05177': 'A*74:01',  TODO - Understand why there's only one HLA allele here

print (vcf_filename)
