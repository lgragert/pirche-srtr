# Clean PIRCHE CSV file

import pandas as pd

# Parse CSV from PIRCHE
pirche_file = pd.read_csv('pirche_result_SRTR_100.csv', sep=',', header=None)

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
PIRCHE_clean.to_csv('clean_priche_100.csv', index=False)



# Open up clean CSV and drop the D and R character from PX_ID
clean_file = pd.read_csv('clean_priche_100.csv', sep=',')
PIRCHE_ids = clean_file['Patient/Donor_ID'].str.split('[DR]', expand=True)
PIRCHE_ids.columns = ['index','PX_ID']
print(PIRCHE_ids.head())
id_list = pd.Series(PIRCHE_ids['PX_ID']).drop_duplicates().to_list()


# Merge Ethnicity Data column to clean data file
directory = '../kidney-outcomes-sfvt/'
grffail_filename = directory + 'SRTR_AA_MM_matrix_grffail_1.txt'
grffail = pd.read_csv(grffail_filename, sep='\t')

# Only need CAN_RACE column with specific PX_ID's
can_race = grffail[["PX_ID", "CAN_RACE"]]
can_race.PX_ID = can_race.PX_ID.astype(str)

can_race = can_race[(can_race['PX_ID'].isin(id_list))]
print(can_race.head())


# Merges the PIRCHE PX_ID with the grffail PX_ID column
pirche_pop = pd.merge(PIRCHE_ids, can_race, how='outer', on='PX_ID')
pirche_pop.drop(columns=['PX_ID','index'])
print(pirche_pop.head())

clean_file.insert(1, 'CAN_RACE', pirche_pop.loc[:, 'CAN_RACE'])
print(clean_file.head())


# Output new CSV, but has population information in it now
clean_file.to_csv('clean_pirche100_pops.csv', index=False)

