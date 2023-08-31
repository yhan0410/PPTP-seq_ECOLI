import numpy as np
import pandas as pd
from scipy.optimize import minimize
from tqdm.notebook import trange,tqdm
import multiprocessing
from math import erf
from functools import partial
from PPTPseq_utils import avoid_zero, remove_outlier_IQR, rescale_linear

# log normal distribution is used for fitting
# x axis data in sort-seq is in log-scale, so gaussian cdf is used here
def gaussian_cdf(x,*param): 
    mu,sigma = param
    if sigma == 0:
        return np.nan
    else:
        return (1 + erf((x - mu) / sigma / np.sqrt(2))) / 2

# minimizing forward KL divergence is equivalent maximum likelihood estimation (MLE)
def forward_KL_Divergence_loss(param, px, FACS_setup):
    # px is the estimator of cell fractions
    mu = param[0]
    sigma = param[1]
    n_bins = len(px)
    # bin-width is equal in log scale in the sort-seq setting
    # the scale of x range is converted later.
    # the boundaries are log10 values

    # This set up is in BD FACSAria II
    if FACS_setup == 'Aria':
        BOUNDARIES = [1.656616775, 1.840292852, 2.069421667, 2.292622975,
                    2.515824283, 2.739025591, 2.9622269, 3.185428208,
                    3.408629516, 3.631830824, 3.855032132, 4.07823344,
                    4.301434748, 4.524636057, 4.747837365]
    # This set up is in BD FACSMelody
    elif FACS_setup == 'Melody':
        BOUNDARIES = [1.75583333, 2.        , 2.24416667, 2.48833333, 
                    2.7325    , 2.97666667, 3.22083333, 3.465     , 
                    3.70916667, 3.95333333, 4.1975    , 4.44166667, 
                    4.68583333, 4.93      , 5.17416667]
    else:
        print('We do not have this setup')
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

def get_fit_para(r, FACS_setup):
    if r.sum()>0:
        n_bin = len(r)
        if FACS_setup == 'Aria':
            x = np.array(
                [1.54730281, 1.735373418, 1.964796054, 2.187997362,
                2.41119867, 2.634399978, 2.857601286, 3.080802595,
                3.304003903, 3.527205211, 3.750406519, 3.973607827,
                4.196809135, 4.420010443, 4.643211752, 5.085286478])
            bounds = [(0.0,0.1),(5,0.5)]
        elif FACS_setup == 'Melody':
            x = np.array(
                [1.63375   , 1.87791667, 2.12208334, 2.36625   ,
                2.61041667, 2.85458334, 3.09875   , 3.34291667,
                3.58708334, 3.83125   , 4.07541667, 4.31958334,
                4.56375   , 4.80791667, 5.05208334, 5.29625])
            bounds = [(0.0,0.1),(5.5,0.7)]
        else:
            print('We do not have this setup')
        # the initial guess
        mean = x.dot(r)
        std = np.sqrt(np.sum((x - mean)**2 * r))
        guess = [mean, std]
        # fit Gaussian distribution using maximum likelihood estimation
        para_forward = minimize(partial(forward_KL_Divergence_loss, FACS_setup = FACS_setup), x0 = guess, args = r, bounds = bounds, method='Nelder-Mead')
        fKLD = forward_KL_Divergence_loss(para_forward['x'],r, FACS_setup)
        return para_forward['x'], fKLD, para_forward['success']
    else:
        return np.array([0.0,0.0]), -1, False

def multiprocess_fitting(estimated_cell_count, FACS_setup):
    estimated_cell_fraction = estimated_cell_count.div(estimated_cell_count.sum(axis = 1), axis=0)
    estimated_cell_fraction = denoise_cell_distribution(estimated_cell_fraction)
    cell_fraction = estimated_cell_fraction.values #calculation in numpy format is faster 
    cf_iter = [cell_fraction[i] for i in range(cell_fraction.shape[0])]
    fit_result = list(tqdm(multiprocessing.get_context('fork').Pool().imap(partial(get_fit_para, FACS_setup = FACS_setup), cf_iter), total=cell_fraction.shape[0]))
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

def exclude_outliners(df, metric, cell_cutoff, FACS_setup):
    mask_df = df.copy()
    if FACS_setup == 'Aria':
        EXP_BOUNDARY = [1.5, 5]
    elif FACS_setup == 'Melody':
        EXP_BOUNDARY = [1.5, 5.5]
    else:
        print('We do not have this setup')
    conditions = ((mask_df['n_cells']>cell_cutoff)  
                & (mask_df['log_mean']<EXP_BOUNDARY[1] * np.log(10))  
                & (mask_df['log_mean']>EXP_BOUNDARY[0] * np.log(10)) 
                & (mask_df['KLD']<1))
    mask_df[metric] = mask_df[metric].where(conditions)
    return mask_df

def get_median_promoter_activity(df):
    activity_matrix = df.reset_index().pivot(index = 'operon', columns='sgRNA', values='log_mean')
    return activity_matrix.median(1)

def get_replicate_stats(df, condition):
    log_GFP_normalized = df.copy()
    log_GFP_normalized.reset_index(inplace = True)
    log_GFP_normalized_nc = log_GFP_normalized[log_GFP_normalized['sgRNA'].isin(['NC_35','NC_82','NC_84','NC_89'])].set_index(['operon', 'sgRNA']).unstack()
    log_GFP_normalized_nc = remove_outlier_IQR(log_GFP_normalized_nc)
    log_GFP_normalized_exp = log_GFP_normalized[~log_GFP_normalized['sgRNA'].isin(['NC_31', 'NC_35','NC_82','NC_84','NC_89'])].set_index(['operon', 'sgRNA'])

    stat_exp = pd.concat([
        log_GFP_normalized_exp.mean(axis=1).rename('mean'),
        log_GFP_normalized_exp.std(axis=1).rename('std'),
        np.exp(log_GFP_normalized_exp).std(axis=1).rename('std_linear'),
        log_GFP_normalized_exp.notna().sum(axis=1).rename('n_rep')], axis=1)
    stat_nc = pd.concat([
        log_GFP_normalized_nc.mean(axis=1).rename('mean'),
        log_GFP_normalized_nc.std(axis=1).rename('std'),
        np.exp(log_GFP_normalized_nc).std(axis=1).rename('std_linear'),
        log_GFP_normalized_nc.notna().sum(axis=1).rename('n_rep')], axis=1)
    stat_exp.to_csv('../Processed_data/normalized_stat_exp_'+condition+'.csv')
    stat_nc.to_csv('../Processed_data/normalized_stat_nc_'+condition+'.csv')
    
replicates = ['Glucose_rep1', 'Glucose_rep2', 'Glucose_rep3', 'LB_rep1', 'LB_rep2', 'Glycerol_rep1', 'Glycerol_rep2_1', 'Glycerol_rep2_2']
bio_replicates = ['Glucose_rep1', 'Glucose_rep2', 'Glucose_rep3', 'LB_rep1', 'LB_rep2', 'Glycerol_rep1', 'Glycerol_rep2']
read_count_dfs = []
for replicate in replicates:
    read_count_file = '../Processed_data/read_count_'+ replicate + '.csv'
    read_count_dfs.append(pd.read_csv(read_count_file, index_col = ['operon', 'sgRNA']))
# merge the reads from two tech replicates for M9 glycerol biological replicate #2
read_count_dfs[6] = read_count_dfs[6].add(read_count_dfs[7], fill_value=0)
read_count_dfs.pop()
# the number of biological replicates in total
num_bio_reps = len(read_count_dfs)

TOTAL_READ_COUNT = pd.read_csv('../Processed_data/total_read_count_summary.csv')
TOTAL_CELL_COUNT = pd.read_csv('../Processed_data/total_cell_count_summary.csv')
# set cell count for bin1 in M9 glycerol to 0 due to potential unwanted mutations.
TOTAL_CELL_COUNT.loc[[5,6],'1'] = 0

for i in range(num_bio_reps):
    estimated_cell_count = read_count_dfs[i].multiply(TOTAL_CELL_COUNT.iloc[i].values).div(TOTAL_READ_COUNT.iloc[i].values)
    if i < 3:
        FACS_setup = 'Aria'
    else:
        FACS_setup = 'Melody'
    raw_fitting_result = multiprocess_fitting(estimated_cell_count, FACS_setup = FACS_setup)
    # filter low quality fitting data
    fitting_result_df = exclude_outliners(raw_fitting_result, metric = 'log_mean', cell_cutoff = 1, FACS_setup = FACS_setup)
    mean_log_GFP = fitting_result_df['log_mean'].rename(bio_replicates[i])
    median_log_GFP_per_promoter = get_median_promoter_activity(fitting_result_df).rename(bio_replicates[i])
    if i == 0:
        mean_log_GFP_df = mean_log_GFP
        median_log_GFP_per_promoter_df = median_log_GFP_per_promoter
    else:
        mean_log_GFP_df = pd.concat([mean_log_GFP_df, mean_log_GFP], axis = 1)
        median_log_GFP_per_promoter_df = pd.concat([median_log_GFP_per_promoter_df, median_log_GFP_per_promoter], axis = 1)
# include all operon-sgRNA combinations
# mean_log_GFP_df.reindex(pd.MultiIndex.from_product(mean_log_GFP_df.index.levels, names=['operon','sgRNA']))

# rescale
bio_replicate_to_be_scaled = ['Glucose_rep1', 'Glucose_rep2', 'LB_rep1', 'Glycerol_rep1']
bio_replicate_references = ['Glucose_rep3', 'Glucose_rep3', 'LB_rep2', 'Glycerol_rep2']
log_GFP_normalized_df = mean_log_GFP_df.copy()
for i in range(len(bio_replicate_to_be_scaled)):
    input_log = mean_log_GFP_df[bio_replicate_to_be_scaled[i]]
    x_log = median_log_GFP_per_promoter_df[bio_replicate_to_be_scaled[i]]
    y_log = median_log_GFP_per_promoter_df[bio_replicate_references[i]]
    log_GFP_normalized_df[bio_replicate_to_be_scaled[i]] = rescale_linear(input_log, x_log, y_log)
    
get_replicate_stats(log_GFP_normalized_df[['Glucose_rep1', 'Glucose_rep2', 'Glucose_rep3']], 'Glu')
get_replicate_stats(log_GFP_normalized_df[['LB_rep1', 'LB_rep2']], 'LB')
get_replicate_stats(log_GFP_normalized_df[['Glycerol_rep1', 'Glycerol_rep2']], 'Gly')
