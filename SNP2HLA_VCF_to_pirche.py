import sys
import pandas
import gzip

vcf_filename = "./SNP2HLA_Imputation/Penn.PTI/chr6.dose.vcf.gz"
# vcf_df = pandas.read_csv(vcf_filename,delimiter='\t',compression="gzip")

sample_IDs = []
sample_ID_HLA_geno_A = {} # store the HLA allele-level genotype for the A locus - for key "05-01308_05-01308" assign a value of "A*01:01+A*11:01"
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
            for x in range(9,1163):
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


        if (hla_locus != "A"):
            continue
        # print (hla_locus)
        print (hla_specificity)

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
        genotype_dict = {}
        genotype_GT_dict = {}
        for i in range(9,1163):
            sample_ID_index = i - 9
            sample_ID = sample_IDs[sample_ID_index]
            # print (sample_ID)
            genotype = fields[i]
            (geno_GT,geno_DS,geno_HDS,geno_GP) = genotype.split(',')
            # TBD - determine what the genotype string definition is - look at SNP2HLA documentation
            (allele_presence_geno,float_unknown_1,float_unknown2) = geno_GT.split(":")
            (allele1_presence,allele2_presence) = allele_presence_geno.split("|")
            if (allele_presence_geno != "0|0"):  # only look at genotypes where the allele was present - either position can be positive
                print ("Sample ID: " + sample_ID)
                print ("GT_1: " + allele1_presence)
                if (allele1_presence == "1"):
                    print (hla_locus)
                    if (hla_locus == "A"):
                        if (sample_ID not in sample_ID_HLA_geno_A): # put in allele if no previous allele
                            sample_ID_HLA_geno_A[sample_ID] = hla_specificity
                        else: # else make into a genotype
                            sample_ID_HLA_geno_A[sample_ID] = sample_ID_HLA_geno_A[sample_ID] + "+" + hla_specificity
                        print ("A locus type present 1")
                        print (sample_ID_HLA_geno_A)            
                print ("GT_2: " + allele2_presence)
                if (allele2_presence == "1"):
                    if (hla_locus == "A"):
                        if (sample_ID not in sample_ID_HLA_geno_A): # put in allele if no previous allele
                            sample_ID_HLA_geno_A[sample_ID] = hla_specificity
                        else: # else make into a genotype
                            sample_ID_HLA_geno_A[sample_ID] = sample_ID_HLA_geno_A[sample_ID] + "+" + hla_specificity
                        print ("A locus type present 2")
                        print (sample_ID_HLA_geno_A)
                # print ("GT: " + geno_GT)    # 1|0:1.000:1.000    0|1:0.949:0.000
                # print ("DS: " + geno_DS)
            genotype_dict[sample_ID] = genotype
            genotype_GT_dict[sample_ID] = geno_GT

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
print (sample_ID_HLA_geno_A)  # should have a HLA genotype for each ID for locus A

# TODO - do a multi-dimensional or nested Python dictionary so you don't need a separate dictionary per locus

# '05-04377_05-04377': 'A*34:02+A*68:02'
# '06-05177_06-05177': 'A*74:01',  TODO - Understand why there's only one HLA allele here

print (vcf_filename)