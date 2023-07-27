import pandas as pd
import re

cohort = "VIENNA"

subjHLA_filename = "Subj_HLA_" + cohort + ".csv"
subjHLA = pd.read_csv(subjHLA_filename, sep=',', names=["PX_ID", "A_1", "A_2", "C_1", "C_2", "B_1", "B_2",
                                                       "DRB1_1", "DRB1_2", "DQA1_1", "DQA1_2",
                                                       "DQB1_1", "DQB1_2", "DPA1_1", "DPA1_2",
                                                       "DPB1_1", "DPB1_2"])
total = len(subjHLA)
print(total)

HLA_list = ["A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]

HLA_dict = {"Summary of " + cohort: ['Missing Loci Total: ', 'Extra Loci Total: ']}
HLA_frac_dict = {"Summary of " + cohort: ['Missing Loci Percent: ', 'Extra Loci Percent: ']}
extra_happairs = pd.DataFrame()
for HLA in HLA_list:
    
    HLA_1 = HLA + "_1"
    HLA_2 = HLA + "_2"
    
    missing1 = subjHLA[HLA_1].isna().sum()
    missing2 = subjHLA[HLA_2].isna().sum()

    missing_total = missing1 + missing2
    missing_perc = format(100* missing_total / total, '.2f')

    combo1 = subjHLA[HLA_1].str.contains('[+]').sum()
    df = subjHLA[subjHLA[HLA_1].str.contains('[+]', na=False)]
    df = df['PX_ID']
    combo2 = subjHLA[HLA_2].str.contains('[+]').sum()
    df2 = subjHLA[subjHLA[HLA_2].str.contains('[+]', na=False)]
    df2 = df2['PX_ID']

    extra_df = pd.concat([df,df2])
    extra_happairs = pd.concat([extra_happairs, extra_df])

    combo = combo1 + combo2
    combo_perc = format(100* combo / total, '.2f')

    HLA_dict[HLA] = [missing_total, combo]
    HLA_frac_dict[HLA] = [missing_perc, combo_perc]


# List of missing loci by PX_ID
missing_df = subjHLA[subjHLA.isna().any(axis=1)]
missing_id = missing_df['PX_ID']
missing_id = missing_id.reset_index(drop=True)

# List of extra loci by PX_ID
extra_happairs = extra_happairs.reset_index(drop=True)
missing_extra = pd.concat([missing_id, extra_happairs], axis=1)
missing_extra.columns = ['Missing HLA', 'Extra HLA']

# Sum of missing more than one loci
more_missing =str(len(subjHLA.loc[subjHLA.isnull().sum(1)>1])) + "/" + str(total)
more_missing_perc = format(100 * len(subjHLA.loc[subjHLA.isnull().sum(1)>1]) /  total, '.2f')
print(more_missing)
print(more_missing_perc)

missing_extra.to_csv("Missing_HLA_in_" + cohort + ".csv", index=False)

HLA_info = pd.DataFrame.from_dict(HLA_dict)
HLA_frac = pd.DataFrame.from_dict(HLA_frac_dict)
HLA_summary = pd.concat([HLA_info, HLA_frac])
print(HLA_summary)

HLA_summary.to_csv("summary_HLA_missing_" + cohort + ".csv", index=False)


