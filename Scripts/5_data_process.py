import numpy as np
import pandas as pd
from scipy.optimize import minimize
from tqdm.notebook import trange,tqdm
import multiprocessing
from math import erf

def calculate_read_counts(mapping_file_path):
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
    return read_count

# log normal distribution is used for fitting
# x axis data in sort-seq is in log-scale, so gaussian cdf is used here
def gaussian_cdf(x,*param): 
    mu,sigma = param
    if sigma == 0:
        return np.nan
    else:
        return (1 + erf((x - mu) / sigma / np.sqrt(2))) / 2

def avoid_zero(x):
    SMALL_NUMBER = 1e-10
    n_bins = 16
    y = (1 - n_bins * SMALL_NUMBER) * x + SMALL_NUMBER
    return y
    
# minimizing forward KL divergence is equivalent maximum likelihood estimation (MLE)
def forward_KL_Divergence_loss(param, px):
    # px is the estimator of cell fractions 
    mu = param[0]
    sigma = param[1]
    n_bins = len(px)
    # bin-width is equal in log scale in the sort-seq setting
    # the scale of x range is converted later.
    # the boundaries are log10 values
    BOUNDARIES = [1.656616775, 1.840292852, 2.069421667, 2.292622975,
                    2.515824283, 2.739025591, 2.9622269, 3.185428208,
                    3.408629516, 3.631830824, 3.855032132, 4.07823344,
                    4.301434748, 4.524636057, 4.747837365]
    x = np.array(BOUNDARIES)
#     avoid 0 in calculation
    px = avoid_zero(px)
    qx = np.zeros(16)
    for i in range(n_bins): # from 0 to 15
        if i==0: # the first bin index 0
            qx[i] = avoid_zero(gaussian_cdf(x[0],mu,sigma) - 0)
        elif i==n_bins-1: # the last bin index 15
            qx[i] = avoid_zero(1 - gaussian_cdf(x[i-1],mu,sigma))
        else:
            qx[i] = avoid_zero(gaussian_cdf(x[i],mu,sigma) - gaussian_cdf(x[i-1],mu,sigma))
    KLD = sum(px * np.log(px/qx))
    return KLD

def denoise_cell_distribution(cell_fraction_df):
    # assuming 5% sorting errors
    # p_i, the observed fraction of cells in bin i in FACS-seq
    # after denoise, the adjusted cell fraction in bin i should be
    # p_adj_i = (p_i - 0.05)/sum(p_i-0.05), i is the bin index, and we have 16 bins here
    # so p_adj_i = 5 * p_i - 0.25
    cell_fraction_df = cell_fraction_df * 5 - 0.25
    # if p_i < 0.05, p_adj_i = 0
    cell_fraction_df = cell_fraction_df.mask(cell_fraction_df<0,0)
    # redo normalization by setting as p_adj_new_i = p_adj_i/sum(p_adj_i)
    cell_fraction_df = cell_fraction_df.div(cell_fraction_df.sum(axis = 1), axis=0)
    return cell_fraction_df

def get_fit_para(r):
    if r.sum()>0:
        n_bin = len(r)
        x = np.array([1.54730281, 1.735373418, 1.964796054, 2.187997362,
            2.41119867, 2.634399978, 2.857601286, 3.080802595,
            3.304003903, 3.527205211, 3.750406519, 3.973607827,
            4.196809135, 4.420010443, 4.643211752, 5.085286478])
        # the initial guess
        mean = x.dot(r)
        std = np.sqrt(np.sum((x - mean)**2 * r))
        guess = [mean, std]
        bounds = [(0.0,0.1),(5,0.5)]
        # fit Gaussian distribution using maximum likelihood estimation
        para_forward = minimize(forward_KL_Divergence_loss, x0 = guess, args = r, bounds = bounds, method='Nelder-Mead')
        fKLD = forward_KL_Divergence_loss(para_forward['x'],r)
#         para_forward = minimize(reverse_KL_Divergence_loss, x0 = guess, args = r, bounds = bounds, method='Nelder-Mead')
#         fKLD = reverse_KL_Divergence_loss(para_forward['x'],r)
        return para_forward['x'], fKLD, para_forward['success']
    else:
        return np.array([0.0,0.0]), -1, False

def multiprocess_fitting(estimated_cell_count):
    estimated_cell_fraction = estimated_cell_count.div(estimated_cell_count.sum(axis = 1), axis=0)
    estimated_cell_fraction = denoise_cell_distribution(estimated_cell_fraction)
    cell_fraction = estimated_cell_fraction.to_numpy() #calculation in numpy format is faster 
    cf_iter = [cell_fraction[i] for i in range(cell_fraction.shape[0])]
    fit_result = list(tqdm(multiprocessing.get_context('fork').Pool().imap(get_fit_para, cf_iter), total=cell_fraction.shape[0]))
    fit_stat = post_fitting_process(fit_result, estimated_cell_count)
    return fit_stat    
    
def post_fitting_process(fit_result, estimated_cell_count):
    logc = np.log(10)
    fit_stat_df = pd.DataFrame(fit_result, index = estimated_cell_count.index, columns = ['fitted_parameters','KLD', 'success'])
    fit_stat_df[['mu_log10','sigma_log10']] = pd.DataFrame(fit_stat_df['fitted_parameters'].tolist(), index = estimated_cell_count.index)
    fit_stat_df['mu_scaled'] = logc * fit_stat_df['mu_log10'] # convert log10 scale to natural log scale
    fit_stat_df['sigma_scaled'] = logc * fit_stat_df['sigma_log10'] # convert log10 scale to natural log scale
    fit_stat_df['log_mean'] = fit_stat_df['mu_scaled'] + fit_stat_df['sigma_scaled']**2 / 2    
    fit_stat_df['n_cells'] = estimated_cell_count.sum(1)
    return fit_stat_df[['mu_scaled','sigma_scaled','log_mean','KLD','n_cells']]

def mask_variants(df, metric, cell_cutoff):
    mask_df = df.copy()
    conditions = ((mask_df['n_cells']>cell_cutoff)  
                & (mask_df['log_mean']<5 * np.log(10))  
                & (mask_df['log_mean']>1.5 * np.log(10)) 
                & (mask_df['KLD']<1))
    mask_df[metric] = mask_df[metric].where(conditions)
    return mask_df
    
# construct a dataframe of read counts
# data from different replicates can be accessed by changing the path
# replicate #1
mapping_file_path_rep1 = '../SAM_files/20210601/filtered_closest_operon_map_all.bed'
# replicate #2
mapping_file_path_rep2 = '../SAM_files/20210823/filtered_closest_operon_map_all.bed'
# replicate #3
mapping_file_path_rep3 = '../SAM_files/20211123/filtered_closest_operon_map_all.bed'

read_count_rep1 = calculate_read_counts(mapping_file_path_rep1)
read_count_rep2 = calculate_read_counts(mapping_file_path_rep2)
read_count_rep3 = calculate_read_counts(mapping_file_path_rep3)
# save the read count data
read_count_rep1.to_csv('../Processed_data/read_count_rep_1.csv', index=False)
read_count_rep2.to_csv('../Processed_data/read_count_rep_2.csv', index=False)
read_count_rep3.to_csv('../Processed_data/read_count_rep_3.csv', index=False)

# total read count constants
TOTAL_READ_COUNT1 = np.array([526274, 1194008, 2265607, 9058702, 20542654, 40805271, 42182137, 53164051, 52619159, 21398173, 11956414, 7061493, 5887774, 1834363, 255317, 77753])
TOTAL_READ_COUNT2 = np.array([9149788, 17773547, 22946265, 23785948, 43608136, 69176247, 72663308, 103891017, 55248775, 35914221, 34111528, 18752959, 10675065, 6528797, 1455479, 303507])
TOTAL_READ_COUNT3 = np.array([6456869, 10039262, 16268575, 20649209, 34314468, 23931989, 39719819, 75408669, 57853559, 37099995, 29507816, 16034625, 8341981, 6904891, 2148394, 408292])
# total cell count constants
TOTAL_CELL_COUNT1 = np.array([50770, 112627, 355751, 1049891, 2437232, 4199762, 4431353, 5001263, 4803172, 3102923, 1690719, 913290, 529518, 270436, 88071, 11023])
TOTAL_CELL_COUNT2 = np.array([339464, 439283, 808952, 1326194, 1719727, 2439567, 2548201, 3087259, 2781019, 2013786, 1240092, 739759, 449023, 240175, 92382, 16416])
TOTAL_CELL_COUNT3 = np.array([1337144, 1397684, 2577177, 3879439, 5540160, 8538119, 10358916, 13482740, 13935055, 10962977, 5710509, 3427873, 1865614, 970827, 460179, 83299])
# c_ij = r_ij * C_j / R_j
read_count_rep1.set_index(['operon', 'sgRNA'], inplace=True)
read_count_rep2.set_index(['operon', 'sgRNA'], inplace=True)
read_count_rep3.set_index(['operon', 'sgRNA'], inplace=True)
estimated_cell_count_rep1 = read_count_rep1.multiply(TOTAL_CELL_COUNT1).div(TOTAL_READ_COUNT1)
estimated_cell_count_rep2 = read_count_rep2.multiply(TOTAL_CELL_COUNT2).div(TOTAL_READ_COUNT2)
estimated_cell_count_rep3 = read_count_rep3.multiply(TOTAL_CELL_COUNT3).div(TOTAL_READ_COUNT3)

# fitting to log-normal distribution
log_GFP1 = multiprocess_fitting(estimated_cell_count_rep1)
log_GFP2 = multiprocess_fitting(estimated_cell_count_rep2)
log_GFP3 = multiprocess_fitting(estimated_cell_count_rep3)
# save fitting results
log_GFP1.to_csv('../Processed_data/fit_denoised_rep1.csv')
log_GFP2.to_csv('../Processed_data/fit_denoised_rep2.csv')
log_GFP3.to_csv('../Processed_data/fit_denoised_rep3.csv')

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

# filter low quality fitting data and merge 
cell_cutoff = 1
metric = 'log_mean'
masked_log_GFP1 = mask_variants(log_GFP1, metric, cell_cutoff)
masked_log_GFP2 = mask_variants(log_GFP2, metric, cell_cutoff)
masked_log_GFP3 = mask_variants(log_GFP3, metric, cell_cutoff)
# calculate median promoter activity across all sgRNA conditions
# this is the reference for rescale
masked_pivot1 = masked_log_GFP1.pivot(index = 'operon', columns='sgRNA', values='log_mean')
masked_pivot2 = masked_log_GFP2.pivot(index = 'operon', columns='sgRNA', values='log_mean')
masked_pivot3 = masked_log_GFP3.pivot(index = 'operon', columns='sgRNA', values='log_mean')
median_log_GFP = pd.concat([masked_pivot1.median(1).rename('rep1'),
                            masked_pivot2.median(1).rename('rep2'),
                            masked_pivot3.median(1).rename('rep3')], axis = 1)
log_GFP = pd.concat([masked_log_GFP1.set_index(['operon', 'sgRNA'])['log_mean'].rename('rep1'),
                     masked_log_GFP2.set_index(['operon', 'sgRNA'])['log_mean'].rename('rep2'),
                     masked_log_GFP3.set_index(['operon', 'sgRNA'])['log_mean'].rename('rep3')], 
                     axis = 1)
log_GFP = log_GFP.reindex(pd.MultiIndex.from_product(log_GFP.index.levels, names=['operon','sgRNA']))
# normalization based on robust linear regression
log_GFP_normalized = log_GFP.copy()
log_GFP_normalized['rep1'] = rescale_linear(log_GFP_normalized['rep1'], median_log_GFP['rep1'], median_log_GFP['rep3'])
log_GFP_normalized['rep2'] = rescale_linear(log_GFP_normalized['rep2'], median_log_GFP['rep2'], median_log_GFP['rep3'])
log_GFP_normalized.to_csv('../Processed_data/normalized_rep_mean_value.csv')

log_GFP_normalized.reset_index(inplace = True)
log_GFP_normalized_nc = log_GFP_normalized[log_GFP_normalized['sgRNA'].isin(['NC_35','NC_82','NC_84','NC_89'])].set_index(['operon', 'sgRNA']).unstack()
log_GFP_normalized_exp = log_GFP_normalized[~log_GFP_normalized['sgRNA'].isin(['NC_31', 'NC_35','NC_82','NC_84','NC_89'])].set_index(['operon', 'sgRNA'])

stat_exp = pd.concat([
    log_GFP_normalized_exp.mean(axis=1).rename('mean'),
    log_GFP_normalized_exp.std(axis=1).rename('std'),
    np.exp(log_GFP_normalized_exp).std(axis=1).rename('std_linear'),
    log_GFP_normalized_exp.notna().sum(axis=1).rename('n_rep')],axis=1)

stat_nc = pd.concat([
    log_GFP_normalized_nc.mean(axis=1).rename('mean'),
    log_GFP_normalized_nc.std(axis=1).rename('std'),
    np.exp(log_GFP_normalized_nc).std(axis=1).rename('std_linear'),
    log_GFP_normalized_nc.notna().sum(axis=1).rename('n_rep')],axis=1)

stat_exp.to_csv('../Processed_data/normalized_stat_exp.csv')
stat_nc.to_csv('../Processed_data/normalized_stat_nc.csv')