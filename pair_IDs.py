# Formatting the new CSV file from VCF data into the PIRCHE format
# by Alyssa Paynter

import pandas as pd
import sys

# pd.set_option('display.max_columns', None)
subjHLA_PennPTI = 'Subj_HLA_PennPTI.csv'
subjHLA = pd.read_csv(subjHLA_PennPTI, sep=',', header=None)
# print (subjHLA)

# Get pair IDs list
pairID_PennPTI = './SNP2HLA_Imputation/Penn.PTI/Penn.PTI.DR.pairs.IDs.csv'
pairID = pd.read_csv(pairID_PennPTI, sep=',')

pairIDs = pairID.dropna()           # get rid of recip data only
recipID = pd.DataFrame()
donorID = pd.DataFrame()
recipID[0] = pairIDs['IID'].astype(str) + "_" + pairIDs['IID'].astype(str)     # format in the vcf file had the name twice
donorID[0] = pairIDs['IID_D'].astype(str) + "_" + pairIDs['IID_D'].astype(str)
# print(recipID)
# print(donorID)

recip = pd.merge(recipID, subjHLA)
donor = pd.merge(donorID, subjHLA)

# Create new, unique names for every recipient and donor
nums = list(range(1, len(pairIDs) + 1))
penn_recip_IDs = pd.DataFrame({0: nums})
penn_recip_IDs[0] = 'R_PennPTI_' + penn_recip_IDs[0].astype(str)
recip.insert(1, 'PX_ID', penn_recip_IDs)
recip = recip.drop(recip.columns[0], axis=1)
# print(recip.head())

penn_donor_IDs = pd.DataFrame({0: nums})
penn_donor_IDs[0] = 'D_PennPTI_' + penn_donor_IDs[0].astype(str)
donor.insert(1, 'PX_ID', penn_donor_IDs)
donor = donor.drop(donor.columns[0], axis=1)
# print(donor.head())

# Append recipient and donor rows with a blank space inbetween
PIRCHE = pd.DataFrame()
blank_data = [['','','','','','','','','','','','','','','','','']]
comma_df = pd.DataFrame(blank_data)
for i in range(len(pairIDs)):

    pirche_recip = recip.iloc[[i]]
    pirche_donor = donor.iloc[[i]]

    PIRCHE_merge = pd.concat([pirche_recip,pirche_donor,comma_df], ignore_index=True)
    PIRCHE = pd.concat([PIRCHE, PIRCHE_merge], ignore_index=True)

# Drop last row which is a blank
PIRCHE = PIRCHE[:-1]
# print(PIRCHE.tail())

PIRCHE.to_csv('impute_PennPTI_pirche.csv', header=False, index=False)



print("Work on PMBB files")
# PMBB is slightly different as some of their donors come from PennPTI data
subjHLA_PMBB = 'Subj_HLA_PMBB.csv'
subjHLA1 = pd.read_csv(subjHLA_PMBB, sep=',', header=None)

# Get pair IDs list
pairID_PMBB = './SNP2HLA_Imputation/PMBB/PMBB.DR.pairs.IDs.csv'
pairID_1 = pd.read_csv(pairID_PMBB, sep=',')

pairIDs_1 = pairID_1.dropna()           # get rid of recip data only
recipID_1 = pd.DataFrame()
donorID_1 = pd.DataFrame()
recipID_1[0] = pairIDs_1['IID'].astype(str) + "_" + pairIDs_1['IID'].astype(str)     # format in the vcf file was the name in a str twice
donorID_1[0] = pairIDs_1['IID_D'].astype(str) + "_" + pairIDs_1['IID_D'].astype(str)

recip1 = pd.merge(recipID_1, subjHLA1)
donor1 = pd.merge(donorID_1, subjHLA1)
donor1 = pd.merge(donorID_1, subjHLA)   # some donors are from PennPTI list
# print(donor1)
# print(recip1)

# Create new, unique names for every recipient and donor
nums = list(range(1, len(recipID_1) + 1))
PMBB_recip_IDs = pd.DataFrame({0: nums})
PMBB_recip_IDs[0] = 'R_PMBB_' + PMBB_recip_IDs[0].astype(str)
recip1.insert(1, 'PX_ID', PMBB_recip_IDs)
recip1 = recip1.drop(recip1.columns[0], axis=1)

nums = list(range(1, len(donorID_1) + 1))
PMBB_donor_IDs = pd.DataFrame({0: nums})
PMBB_donor_IDs[0] = 'D_PMBB_' + PMBB_donor_IDs[0].astype(str)
donor1.insert(1, 'PX_ID', PMBB_donor_IDs)
donor1 = donor1.drop(donor1.columns[0], axis=1)

# Append recipient and donor rows with a blank space inbetween
PIRCHE1 = pd.DataFrame()
for i in range(len(pairIDs_1)):

    pirche_recip1 = recip1.iloc[[i]]
    pirche_donor1 = donor1.iloc[[i]]

    PIRCHE_merge1 = pd.concat([pirche_recip1,pirche_donor1,comma_df], ignore_index=True)
    PIRCHE1 = pd.concat([PIRCHE1, PIRCHE_merge1], ignore_index=True)

# Drop last row which is a blank
PIRCHE1 = PIRCHE1[:-1]
# print(PIRCHE1.tail())

PIRCHE1.to_csv('impute_PMBB_pirche.csv', header=False, index=False)

