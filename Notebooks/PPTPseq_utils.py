import numpy as np
import pandas as pd

def update_tf_name(df):
    sgRNA_name_change = pd.read_csv('../External_data/gene_annotation/tf_name_change.txt', 
                        sep="\t", names = ['new','synonym'])
    name_change_dict = pd.Series(sgRNA_name_change.new.values, index=sgRNA_name_change.synonym).to_dict()
    rows_need_update = df['tf_gene'].isin(sgRNA_name_change.synonym)
    df.loc[rows_need_update,'tf_gene'] = df.loc[rows_need_update,'tf_gene'].map(name_change_dict)
    return df

# Read tf-operon network from regulonDB
def read_tf_operon_network():
    file_path = '../External_data/RegulonDB/network_tf_operon.txt'
    column_name = ['tf', 'operon', 'effect', 'evidence', 'evidence_type']
    tf_operon_net = pd.read_csv(file_path, sep="\t", names = column_name, skiprows = 35)
    # Remove the information in [] in the operon column
    tf_operon_net['operon'] = tf_operon_net['operon'].str.extract(r"(^\S+)\[", expand=False)

    tf_operon_net_1to1_df = tf_operon_net[tf_operon_net['tf'].isin(homo_tf['TF_name'])].reset_index(drop = True)
    tf_operon_net_2to1_df = tf_operon_net[tf_operon_net['tf'].isin(hetero_tf['TF_name'])].reset_index(drop = True)
    tf_operon_net_2to1_df_copy = tf_operon_net_2to1_df.copy()

    # Find TF's gene name
    file_path = '../External_data/gene_annotation/TF_gene_map.txt'
    column_name = ['TF_name', 'tf_gene_name', 'bnumber', 'num_gene']
    tf_gene_map_df = pd.read_csv(file_path, sep="\t", names = column_name)
    homo_tf = tf_gene_map_df[tf_gene_map_df['num_gene']=='single']
    # For hetero TF, a TF-operon interaction is split into two interactions
    hetero_tf = tf_gene_map_df[tf_gene_map_df['num_gene']=='double']
    hetero_tf1 = hetero_tf.drop_duplicates(subset = ['TF_name'], keep = 'first')
    hetero_tf2 = hetero_tf.drop_duplicates(subset = ['TF_name'], keep = 'last')
    tf_gene_map1 = pd.Series(homo_tf.tf_gene_name.values,index=homo_tf.TF_name).to_dict()
    tf_gene_map2 = pd.Series(hetero_tf1.tf_gene_name.values,index=hetero_tf1.TF_name).to_dict()
    tf_gene_map3 = pd.Series(hetero_tf2.tf_gene_name.values,index=hetero_tf2.TF_name).to_dict()

    # Find gene name for every TF. 
    tf_operon_net_1to1_df['tf_gene'] = tf_operon_net_1to1_df['tf'].map(tf_gene_map1)
    tf_operon_net_2to1_df['tf_gene'] = tf_operon_net_2to1_df['tf'].map(tf_gene_map2)
    tf_operon_net_2to1_df_copy['tf_gene'] = tf_operon_net_2to1_df_copy['tf'].map(tf_gene_map3)
    tf_operon_net_full = pd.concat([tf_operon_net_1to1_df,
                                    tf_operon_net_2to1_df,
                                    tf_operon_net_2to1_df_copy]).reset_index(drop = True)

    tf_operon_net_full.loc[tf_operon_net_full.duplicated(subset = ['operon','tf_gene'], keep=False),'effect']='+-'
    tf_operon_net_full = tf_operon_net_full.sort_values(by = 'evidence_type').drop_duplicates(subset = ['operon','tf_gene'], keep='first')
    tf_operon_net_full = update_tf_name(tf_operon_net_full)
    return tf_operon_net_full

# Read gene operon table
def read_operon_table():
    operon_gene_table = pd.read_csv('../External_data/gene_annotation/gene_operon_map.txt', sep = "\t", names = ['ispseudo','start','end','strand','gene_name','type','bnumber','operon_name','operon_id','is_first_gene'])
    return operon_gene_table



# Read COG annotation
def read_COG_anno_ECO():
    #load COG annotation for E. coli operons
    file_path = '../External_data/gene_annotation/operon_cog_map.txt'
    column_name = ['operon', 'COG_cat','count']
    ECO_COG_df = pd.read_csv(file_path, sep="\t", names = column_name)

    #load COG annotation for conversion from category to full name
    file_path = '../External_data/gene_annotation/fun-20.tab.txt'
    column_name = ['COG_cat', 'color','name']
    COG_func_df = pd.read_csv(file_path, sep="\t", names = column_name)
    return ECO_COG_df,COG_func_df

# known binding sites in RegulonDB v10.9
def read_known_binding_sites():
    file_name = '../External_data/RegulonDB/BindingSiteSet_organized.txt'
    Binding_site_df = pd.read_csv(file_name, sep="\t")
    Binding_site_df = update_tf_name(Binding_site_df)
    Binding_site_df['abs_pos_TSS'] = abs(Binding_site_df['position_to_TSS'])
    Binding_site_df['known_TFBS'] = 1
    Binding_site_unique_df = Binding_site_df.sort_values(by = 'abs_pos_TSS').drop_duplicates(subset=['tf_gene','operon'], keep='first')
    return Binding_site_unique_df[['tf_gene', 'operon', 'position_to_TSS', 'gene_expression_effect', 'known_TFBS']]

# Data source: regulonDB
# for i in {365..528}; do curl https://regulondb-datasets.ccg.unam.mx/wdps/RHTECOLIBSD00${i}/authorData/cvs >>gselex_dataset.csv; done;
def read_gSELEX():
    file_name = '../External_data/RegulonDB/gSELEX.tsv'
    SELEX_df = pd.read_csv(file_name, sep="\t", usecols = [0,2,3,4,5])
    SELEX_df = update_tf_name(SELEX_df)
    SELEX_df['gSELEX'] = 1
    SELEX_unique_df = SELEX_df.sort_values(by = 'Peak Intensity Fold Change/Binding intensity (%)').drop_duplicates(subset=['tf_gene','operon'], keep='last')
    return SELEX_unique_df

# Data source: regulonDB
# for i in {336..364}; do curl https://regulondb-datasets.ccg.unam.mx/wdps/RHTECOLIBSD00${i}/authorData/cvs >>chip-seq_dataset.csv; done;
def read_ChIP_seq():
    file_name = '../External_data/RegulonDB/ChIPseq.tsv'
    ChIPseq_df = pd.read_csv(file_name, sep="\t")
    ChIPseq_df = update_tf_name(ChIPseq_df)
    ChIPseq_df = ChIPseq_df.drop_duplicates(subset=['tf_gene','operon']).reset_index(drop = True)
    ChIPseq_df['ChIP_seq'] = 1
    return ChIPseq_df[['tf_gene','operon','ChIP_seq']]

# Data source: regulonDB
# for i in {242..335}; do curl https://regulondb-datasets.ccg.unam.mx/wdps/RHTECOLIBSD00${i}/authorData/cvs >>chip-exo_dataset.csv; done;
def read_ChIP_exo():
    file_name = '../External_data/RegulonDB/ChIPexo.tsv'
    ChIPexo_df = pd.read_csv(file_name, sep="\t")
    ChIPexo_df = update_tf_name(ChIPexo_df)
    ChIPexo_df = ChIPexo_df.drop_duplicates(subset=['tf_gene','operon']).reset_index(drop = True)
    ChIPexo_df['ChIP_exo'] = 1
    return ChIPexo_df[['tf_gene','operon','ChIP_exo']]

# DAP-seq results (Baumgart et al., 2021)
def read_DAP_seq():
    file_name = '../External_data/DAP_operon.bed'
    column_name = ['site_start','site_end','tf_gene','FC_DAP','operon']
    DAP_df = pd.read_csv(file_name, sep="\t", names = column_name, usecols = [1, 2, 3, 6, 12])
    DAP_df['tf_gene'] = DAP_df['tf_gene'].str.extract(r"(^\S+?)_", expand=False)
    DAP_df = update_tf_name(DAP_df)
    DAP_df['DAP_seq'] = 1
    num_binding_map = (DAP_df['tf_gene']+'_'+DAP_df['operon']).value_counts().to_dict()
    DAP_df['num_binding_site'] = (DAP_df['tf_gene']+'_'+DAP_df['operon']).map(num_binding_map)
    DAP_unique_df = DAP_df.sort_values(by = 'FC_DAP').drop_duplicates(subset=['tf_gene','operon'], keep='last')
    return DAP_unique_df

# Transcription start site (TSS) information
def read_TSS_map():
    file_name = '../External_data/gene_annotation/TSS_operon.bed'
    column_name = ['promoter_name','pos_1','operon','strand']
    TSS_df = pd.read_csv(file_name, sep="\t", names = column_name, usecols = [3, 6, 10, 12])
    # using the closest TSS to the start codon
    TSS_unique_df = pd.concat(
        [TSS_df[TSS_df['strand']=='+'].sort_values(by = 'pos_1').drop_duplicates(subset=['operon'], keep='last'),
         TSS_df[TSS_df['strand']=='-'].sort_values(by = 'pos_1').drop_duplicates(subset=['operon'], keep='first')])
    TSS_unique_df = TSS_unique_df.sort_values(by = 'pos_1').reset_index(drop = True)
    operon_to_TSS_map = pd.Series(TSS_unique_df.pos_1.values, index=TSS_unique_df.operon).to_dict()
    operon_to_strand_map = pd.Series(TSS_unique_df.strand.values, index=TSS_unique_df.operon).to_dict()
    return operon_to_TSS_map, operon_to_strand_map

def calculate_relative_distance_to_TSS(df):
    df['dist1'] = df['site_start'] - df['TSS']
    df['dist2'] = df['site_end'] - df['TSS']
    df['rel_pos'] = (df['dist1']+df['dist2'])/2
    # To calculate minimal distance to TSS
    # df['rel_pos'] = 0
    # df.loc[(df['dist1']<0) & (df['dist2']<0),'rel_pos'] = df.loc[(df['dist1']<0) & (df['dist2']<0),'dist2']
    # df.loc[(df['dist1']>0) & (df['dist2']>0),'rel_pos'] = df.loc[(df['dist1']>0) & (df['dist2']>0),'dist1']
    df.loc[df['operon_strand']=='-','rel_pos'] = - df.loc[df['operon_strand']=='-','rel_pos']
    return df['rel_pos']

# TF protein expression level by Ribo-seq (Li et al., 2014)
def read_Ribo_seq():
    file_name = '../External_data/Li2014_minimal.txt'
    column_name = ['gene_name','expression_level','confidence']
    Protein_level_df = pd.read_csv(file_name, names = column_name, sep="\t", skiprows = 1)
    gene_to_expression_map = pd.Series(Protein_level_df.expression_level.values,index=Protein_level_df.gene_name).to_dict()
    return gene_to_expression_map

# Escherichia coli proteome (Schmidt et al., 2016)
def read_proteome_MS():
    file_name = '../External_data/Schmidt2016.tsv'
    Protein_level_df = pd.read_csv(file_name, sep="\t")
    gene_to_expression_map = Protein_level_df.set_index('Gene')[['Glucose','LB','Glycerol']].to_dict()
    return gene_to_expression_map

# p value correction for multiple hypothesis testing.
# q-value method (Storey and Tibshirani, 2003)
from scipy.interpolate import UnivariateSpline
def qvalue(pvals):
    m, pvals = len(pvals), np.asarray(pvals)
    ind = np.argsort(pvals)
    rev_ind = np.argsort(ind)
    pvals = pvals[ind]

    # Estimate proportion of features that are truly null.
    kappa = np.arange(0, 0.96, 0.01)
    pik = [sum(pvals > k) / (m*(1-k)) for k in kappa]
    cs = UnivariateSpline(kappa, pik, k=3, s=None, ext=0)
    pi0 = float(cs(1.))
    print('The estimated proportion of truly null features is %.3f' % pi0)

    if (pi0 < 0 or pi0 > 1):
        pi0 = 1
        print('Smoothing estimator did not converge in [0, 1]')

    # Compute the q-values.
    qvals = np.zeros(np.shape(pvals))
    qvals[-1] = pi0*pvals[-1]
    for i in np.arange(m-2, -1, -1):
        qvals[i] = min(pi0*m*pvals[i]/float(i+1), qvals[i+1])

    # Order the q-values according to the original order of the p-values.
    qvals = qvals[rev_ind]
    return qvals

# Benjamini-Hochberg method
def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]

# using IQR method to remove outliers in negative control samples
def remove_outlier_IQR(df):
    n_col = df.shape[1]
    df_filtered = df.copy()
    Q1 = df.quantile(0.25, axis = 1)
    Q3 = df.quantile(0.75, axis = 1)
    IQR = Q3 - Q1
    for i_col in range(n_col):
        df_filtered.iloc[:,i_col]=df.iloc[:,i_col].where((df.iloc[:,i_col]>(Q1-1.5*IQR))&(df.iloc[:,i_col]<(Q3+1.5*IQR)))
    return df_filtered

def avoid_zero(x):
    SMALL_NUMBER = 1e-10
    n_bins = 16
    y = (1 - n_bins * SMALL_NUMBER) * x + SMALL_NUMBER
    return y

# normalize scale for different replicates
from sklearn import linear_model
def rescale_linear(input_log, x_log, y_log):
    #rescale input from linear scale of x to linear scale of y
    x_lin = np.exp(x_log)
    y_lin = np.exp(y_log)
    notna_idx = (x_lin.notna())&(y_lin.notna())
    #robust_linear_regression
    model = linear_model.HuberRegressor()
    model.fit(x_lin[notna_idx].to_numpy().reshape(-1, 1), y_lin[notna_idx])
    output_log = np.log(model.coef_[0] * np.exp(input_log) + model.intercept_)
    return output_log

CONSISTENCE_CUTOFF = 0.7