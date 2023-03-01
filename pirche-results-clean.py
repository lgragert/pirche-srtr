# Clean PIRCHE CSV file

import pandas as pd
import sys

pirche_output_file = sys.argv[1]

# Parse CSV from PIRCHE

pirche_file = pd.read_csv(pirche_output_file, sep=',', header=None)

# Remove header row
pirche_file.drop(index=pirche_file.index[0], axis=0, inplace=True)
new_header = pirche_file.iloc[0]
pirche_file = pirche_file[1:]
pirche_file.columns = new_header

# Remove blank lines
PIRCHE_rows = pirche_file.dropna(how='all')
PIRCHE = PIRCHE_rows.dropna(axis='columns',thresh=2)
PIRCHE_clean = PIRCHE.reset_index(drop=True)

# Output new CSV
PIRCHE_clean.to_csv('clean_pirche_100.csv', index=False)



# Open up clean CSV to add population information
clean_file = pd.read_csv('clean_pirche_100.csv', sep=',')
clean_file.rename(columns={'Patient/Donor_ID': 'PX_ID'}, inplace=True)
PIRCHE_ids = clean_file['PX_ID']
id_index = clean_file['PX_ID'].str.split("[DR]", expand=True)
id_list = pd.Series(id_index[1]).drop_duplicates().to_list()


# Merge Population Data column to clean data file
directory = '../kidney-outcomes-sfvt/'
grffail_filename = directory + 'SRTR_AA_MM_matrix_grffail_1.txt'
grffail = pd.read_csv(grffail_filename, sep='\t')

# Only need CAN_RACE and DON_RACE column with specific PX_ID's
pop_info = grffail[["PX_ID", "CAN_RACE", "DON_RACE"]]
pop_info.PX_ID = pop_info.PX_ID.astype(str)


can_pop = pop_info[(pop_info['PX_ID'].isin(id_list))]
can_pop = can_pop.drop(columns=['DON_RACE'])
can_pop.rename(columns={'CAN_RACE': 'RACE'}, inplace=True)
can_pop['PX_ID'] = "R" + can_pop['PX_ID']
don_pop = pop_info[(pop_info['PX_ID'].isin(id_list))]
don_pop = don_pop.drop(columns=['CAN_RACE'])
don_pop.rename(columns={'DON_RACE': 'RACE'}, inplace=True)
don_pop['PX_ID'] = "D" + don_pop['PX_ID']


merge_pop = pd.concat([can_pop, don_pop], ignore_index=True)
print(merge_pop.head())

PIRCHE_pop = pd.merge(PIRCHE_ids, merge_pop, how='outer', on='PX_ID')
PIRCHE_pop.drop(columns=['PX_ID'])
print(PIRCHE_pop.head())

clean_file.insert(1, 'RACE', PIRCHE_pop.loc[:, 'RACE'])
print(clean_file.head())

clean_file.to_csv('clean_pirche100_pops.csv', index=False)

