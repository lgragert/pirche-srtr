#!/usr/bin/env python
# Converting SRTR data to PRICHE format

import random
import pandas as pd


# Choose one haplotype pair for patient and one pair for donor independently
# Use weighted choice from https://scaron.info/blog/python-weighted-choice.html
def weighted_choice(seq, weights):
    assert len(weights) == len(seq)
    assert abs(1. - sum(weights)) < 1e-6

    x = random.random()
    for i, elmt in enumerate(seq):
        if x <= weights[i]:
            return elmt
        x -= weights[i]


# Open up file
freqs_dir = 'NEMO-python/ACBDRBXDRB1DQA1DQB1DPA1DPB1_GLOBAL'

pop = ['AFA','API','CAU','HIS','NAM','UNK']

SRTR = pd.DataFrame()
for pops in pop:
    srtr_filename = "./" + freqs_dir + "/impute.srtr." + pops + ".csv.gz"
    
    SRTR_file = pd.read_csv(srtr_filename, delimiter=',',header=None, compression='gzip')
    SRTR = pd.concat([SRTR,SRTR_file], ignore_index=True)

# Subset the columns with respecitive terms
col_names = ["PX_ID", "RANK", "HAPPAIR_1", "HAPPAIR_2", "FREQ"]
SRTR.columns = col_names

# pd.set_option('display.max_columns', None)
# print(SRTR.head())

num_rows = len(SRTR)
print("Number of Rows: ", num_rows)  # 3490955 rows in AFA
# print(SRTR.head())


# Create 2 data frames for all donors and recips using filter criteria
DONOR = SRTR[(SRTR['PX_ID'].str.startswith('D'))]
# DONOR = DONOR.rename(columns={'PX_ID':'DONOR_ID'})

RECIP = SRTR[(SRTR['PX_ID'].str.startswith('R'))]
# RECIP = RECIP.rename(columns={'PX_ID':'RECIP_ID'})

# print(DONOR.head())
# print(RECIP.head())


# Create a dictionary of pair IDs (without the character)
#        Output the pairs
id_index = SRTR['PX_ID'].str.split("[DR]", expand= True)
print('Counted Values for Each ID: ', id_index.value_counts())
id_list = pd.Series(id_index[1]).drop_duplicates().to_list()

id_total = {}
for i in range(len(id_list)):
    id_total[i+1] = id_list[i]
# print("Happair ID Total: ", id_total)


# Compute cumulative genotype frequency totals per subject
# Taken but modified from kidney-outcomes-sfvt/aa_mm_biopython_runmatch_9loc.py
happair_id_total = {}  # this includes both D and R characters
happair_probs = {}
happair_hla = {}
for r in range(len(SRTR)):

    ID = SRTR.iloc[r,0]        # PX_ID column
    FREQ = SRTR.iloc[r,4]      # FREQ column
    happair_freq = float(FREQ)

    if ID not in happair_id_total:
        happair_id_total[ID] = 0
    happair_id_total[ID] += float(happair_freq)


# Renormalization loop
for r in range(len(SRTR)):

    ID = SRTR.iloc[r,0]        # PX_ID column
    FREQ = SRTR.iloc[r,4]      # FREQ column
    HAPPAIR_1 = SRTR.iloc[r,2] # HAPPAIR_1 column
    HAPPAIR_2 = SRTR.iloc[r,3] # HAPPAIR_2 column

    happair_freq = float(FREQ)
    happair = HAPPAIR_1 + "+" + HAPPAIR_2

    if ID not in happair_probs:
        happair_probs[ID] = list()
    happair_probs[ID].append(happair_freq / happair_id_total[ID])
    if ID not in happair_hla:
        happair_hla[ID] = list()
    happair_hla[ID].append(happair)



# Loop through IDs to select the row from the donor table and recip table
happair_selected_donor = {}
happair_selected_recip = {}
for x in id_total:
    id_recip = "R" + id_total[x]
    id_donor = "D" + id_total[x]

    if id_recip not in happair_hla:
        print("Missing ID from Imputation Output: " + id_recip)
        continue
    if id_donor not in happair_hla:
        print("Missing ID from Donor Output: " + id_donor)
        continue

    haplist = happair_hla[id_recip]
    problist = happair_probs[id_recip]

    happair_recip = weighted_choice(haplist,problist)

    haplist = happair_hla[id_donor]
    problist = happair_probs[id_donor]

    happair_donor = weighted_choice(haplist, problist)

    happair_selected_recip[id_recip] = happair_recip
    happair_selected_donor[id_donor] = happair_donor



# Make dictionary into DataFrame
happair_recip_df = pd.DataFrame.from_dict(happair_selected_recip, orient='index', columns=['HAPPAIRS'])
happair_recip_df = happair_recip_df.rename_axis('PX_ID').reset_index()
happair_donor_df = pd.DataFrame.from_dict(happair_selected_donor, orient='index', columns=['HAPPAIRS'])
happair_donor_df = happair_donor_df.rename_axis('PX_ID').reset_index()


# Append recip and donor rows together. print a comma by itself to make PIRCHE
PIRCHE = pd.DataFrame(columns=['PX_ID','HAPPAIRS'])

for i in id_total:
    id_recip = 'R' + id_total[i]
    id_donor = 'D' + id_total[i]

    if id_recip not in happair_hla:
        continue
    if id_donor not in happair_hla:
        continue
    
    blank_data = [['','']]
    comma_df = pd.DataFrame(blank_data, columns=['PX_ID','HAPPAIRS'])

    pirche_recip = happair_recip_df.loc[happair_recip_df['PX_ID'] == id_recip]
    pirche_donor = happair_donor_df.loc[happair_donor_df['PX_ID'] == id_donor]

    PIRCHE_merge = pd.concat([pirche_recip,pirche_donor,comma_df], ignore_index=True)
    PIRCHE = pd.concat([PIRCHE, PIRCHE_merge], ignore_index=True)


# Separate haplotypes
PIRCHE[['HAPPAIR_1', 'HAPPAIR_2']] = PIRCHE.HAPPAIRS.str.split('+', expand=True)
PIRCHE[["A_1", "C_1", "B_1", "DRB3_1", "DRB1_1", "DQA1_1", "DQB1_1", "DPA1_1", "DPB1_1"]] = PIRCHE.HAPPAIR_1.str.split('~',expand = True)
PIRCHE[["A_2", "C_2", "B_2", "DRB3_2", "DRB1_2", "DQA1_2", "DQB1_2", "DPA1_2", "DPB1_2"]] = PIRCHE.HAPPAIR_2.str.split('~',expand = True)


# Drop columns we don't need for PIRCHE
# PIRCHE does not need: DRB3, DQA1, DPA1, and DPB1
# Only columns you need: PX_ID, A_1, C_1, B_1, DRB1_1, DQB1_1, A_2, C_2, B_2, DRB1_2, DQB1_2
PIRCHE = PIRCHE.drop(columns=['HAPPAIRS', 'HAPPAIR_1', 'HAPPAIR_2', 'DRB3_1', 'DQA1_1', 'DPA1_1', 'DPB1_1', 'DRB3_2', 'DQA1_2', 'DPA1_2', 'DPB1_2'])
PIRCHE = PIRCHE[:-1]  # need to drop the last set of commas
print(PIRCHE.head())


# Go from DataFrame to PIRCHE format CSV
PIRCHE.to_csv('impute_srtr_priche.csv', header=False, index=False)

