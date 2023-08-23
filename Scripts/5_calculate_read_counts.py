import pandas as pd

def calculate_read_counts(mapping_file_path, read_count_file):
    # load raw data
    column_name = ['count','sgRNA','bin','operon']
    raw_data = pd.read_csv(mapping_file_path, sep="\t", names = column_name, usecols=[4, 6, 7, 11])
    # sum up the count if the columns promoter,sgRNA,bin are the same
    read_count = raw_data.groupby(['operon', 'sgRNA', 'bin']).agg({'count':'sum'})
    # use bin number in column bin as new columns
    read_count = read_count.pivot_table(index=['operon', 'sgRNA'], columns='bin', values='count')
    # add the non-existing promoter-sgRNA pair rows back with NaN as count for every bin
    read_count = read_count.reindex(pd.MultiIndex.from_product(read_count.index.levels, names=['operon','sgRNA'])).reset_index()
    # replace NaN with 0.0
    read_count = read_count.fillna(0)
    read_count.to_csv(read_count_file, index=False)
    
# first three are M9 glucose replicates, the next two are LB replicates, the final three are M9 glycerol replicates
# Glycerol_rep2_2 is a NGS technical replicates of Glycerol_rep2
folder_names = ['Glucose_rep_1', 'Glucose_rep_2', 'Glucose_rep_3', 'LB_rep1', 'LB_rep2', 'Glycerol_rep1', 'Glycerol_rep2', 'Glycerol_rep2_2']
for folder_name in folder_names:
    mapping_file_path = '../SAM_files/' + folder_name + '/filtered_closest_operon_map_all.bed'
    read_count_file = '../Processed_data/read_count_'+ folder_name + '.csv'
    calculate_read_counts(mapping_file_path, read_count_file)
