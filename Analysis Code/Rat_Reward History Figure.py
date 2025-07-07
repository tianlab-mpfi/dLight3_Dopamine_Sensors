# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 15:16:38 2024

@author: tanne
"""

# %% imports

import utils
import plot_utils
import fp_utils
import numpy as np
import beh_analysis_helpers as bah
import fp_analysis_helpers as fpah
from fp_analysis_helpers import Alignment as Align
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MultipleLocator
import pandas as pd
import pickle
import seaborn as sb

import statsmodels.api as sm
import statsmodels.stats.multitest as smm

# %% Load Data

with open('beh_data.pkl', 'rb') as f:
    sess_data = pickle.load(f)

with open('fp_data.pkl', 'rb') as f:
    fp_data = pickle.load(f)

subj_ids = np.unique(sess_data['subjid'])
sess_ids = sess_data.groupby('subjid')['sessid'].to_dict()
implant_info = fp_data['implant_info']
fp_data = fpah.process_fp_data(fp_data['fp_data'], sess_ids)

# get aggregated session information
sess_info = sess_data.groupby('subjid').agg(n_sess=('sessid', 'nunique'), n_trials=('sessid', 'size')).reset_index()
sess_info = pd.concat([sess_info, 
                       pd.DataFrame({'subjid': ['Avg'],
                                     'n_sess': [sess_info['n_sess'].mean()],
                                     'n_trials': [sess_info['n_trials'].mean()]})], ignore_index=True)

# %% Build signal matrices aligned to alignment points

signal_types = ['dff_iso', 'z_dff_iso']
alignments = [Align.cue, Align.reward]
regions = ['DMS', 'PL']
xlims = {Align.cue: {'DMS': [-1,2], 'PL': [-3,5]},
         Align.reward: {'DMS': [-1,2], 'PL': [-3,20]}}

aligned_signals = {subjid: {sessid: {sig_type: {align: {region: [] for region in regions} 
                                                for align in alignments}
                                     for sig_type in signal_types}
                            for sessid in sess_ids[subjid]} 
                   for subjid in subj_ids}

for subj_id in subj_ids:
    for sess_id in sess_ids[subj_id]:
        if sess_id in fpah.__sess_ignore:
            continue

        sub_fp_data = fp_data[subj_id][sess_id]

        trial_data = sess_data[sess_data['sessid'] == sess_id]

        ts = sub_fp_data['time']
        trial_start_ts = sub_fp_data['trial_start_ts'][:-1]
        cue_ts = trial_start_ts + trial_data['response_cue_time']
        cpoke_out_ts = trial_start_ts + trial_data['cpoke_out_time']
        resp_ts = trial_start_ts + trial_data['response_time']
        reward_ts = trial_start_ts + trial_data['reward_time']

        for signal_type in signal_types:
            for align in alignments:
                match align:
                    case Align.cue:
                        align_ts = cue_ts
                        # mask out reward
                        reward_ts_mask = reward_ts.to_numpy(copy=True)
                        reward_ts_mask[~trial_data['rewarded']] = np.nan
                        if np.isnan(reward_ts_mask[-1]):
                            reward_ts_mask[-1] = np.inf
                        reward_ts_mask = pd.Series(reward_ts_mask).bfill().to_numpy()
                        
                        mask_lims = np.hstack((np.full_like(align_ts, 0)[:, None], reward_ts_mask[:, None]))
                        
                    case Align.reward:
                        align_ts = reward_ts
                        # mask out next reward
                        next_reward_ts = reward_ts[1:].to_numpy(copy=True)
                        next_reward_ts[~trial_data['rewarded'][1:]] = np.nan
                        next_reward_ts = pd.Series(np.append(next_reward_ts, np.inf))
                        next_reward_ts = next_reward_ts.bfill().to_numpy()

                        mask_lims = np.hstack((np.zeros_like(align_ts)[:, None], next_reward_ts[:, None]))
                
                for region in sub_fp_data['processed_signals'].keys():
                    if region in regions:
                        signal = sub_fp_data['processed_signals'][region][signal_type]

                        lims = xlims[align][region]
                        
                        mat, t = fp_utils.build_signal_matrix(signal, ts, align_ts, -lims[0], lims[1], mask_lims=mask_lims)
                        aligned_signals[subj_id][sess_id][signal_type][align][region] = mat

aligned_signals['t'] = {align: {region: [] for region in regions} for align in alignments}
dt = fp_data['dec_info']['decimated_dt']
for align in alignments:
    for region in regions:
        aligned_signals['t'][align][region] = np.arange(xlims[align][region][0], xlims[align][region][1]+dt, dt)

# %% Analyze aligned signals

rew_n_back = 10
rew_rate_n_back = 3
bah.calc_rew_rate_hist(sess_data, n_back=rew_rate_n_back, kernel='uniform')

# get bins output by pandas for indexing
# make sure 0 is included in the first bin, intervals are one-sided
n_rew_hist_bins = 4
rew_hist_bin_edges = np.linspace(-0.001, 1.001, n_rew_hist_bins+1)
rew_hist_bins = pd.IntervalIndex.from_breaks(rew_hist_bin_edges)
rew_hist_bin_strs = {b:'{:.0f}-{:.0f}%'.format(abs(b.left)*100, b.right*100) for b in rew_hist_bins}

alignments = [Align.cue, Align.reward] #  
signal_types = ['z_dff_iso', 'dff_iso'] # 

norm_baseline = True
analyze_peaks = True

filter_props = {Align.cue: {'DMS': {'filter': True, 'use_filt_signal_props': False, 'cutoff_f': 8},
                            'PL': {'filter': True, 'use_filt_signal_props': True, 'cutoff_f': 1}},
                Align.reward: {'DMS': {'filter': True, 'use_filt_signal_props': False, 'cutoff_f': 8},
                            'PL': {'filter': True, 'use_filt_signal_props': True, 'cutoff_f': 1}}}

peak_find_props = {Align.cue: {'DMS': {'min_dist': 0.05, 'peak_tmax': 0.45, 'peak_edge_buffer': 0.08, 'lim_peak_width_to_edges': True},
                               'PL': {'min_dist': 0.2, 'peak_tmax': 1.5, 'peak_edge_buffer': 0.2, 'lim_peak_width_to_edges': True}},
                   Align.reward: {'DMS': {'min_dist': 0.05, 'peak_tmax': 0.45, 'peak_edge_buffer': 0.08, 'lim_peak_width_to_edges': True},
                                  'PL': {'min_dist': 0.5, 'peak_tmax': 3.5, 'peak_edge_buffer': 0.2, 'lim_peak_width_to_edges': False}}}

sides = ['contra', 'ipsi']

t = aligned_signals['t']
stacked_signals = {s: {a: {r: {} for r in regions} 
                       for a in alignments} 
                   for s in signal_types}

if analyze_peaks:
    peak_metrics = []

reward_times = {r: {} for r in regions}

def stack_mat(stacked_mats, key, mat):
    if not key in stacked_mats:
        stacked_mats[key] = np.zeros((0, mat.shape[1]))
    else:
        stacked_mats[key] = np.vstack((stacked_mats[key], mat))

for subj_id in subj_ids:
    print('Analyzing peaks for subj {}'.format(subj_id))
    for sess_id in sess_ids[subj_id]:
        
        trial_data = sess_data[sess_data['sessid'] == sess_id]
        rewarded = trial_data['rewarded'].to_numpy()
        responded = ~np.isnan(trial_data['response_time']).to_numpy()
        rew_hist = pd.cut(trial_data['rew_rate_hist_all'], rew_hist_bins)
        choice = trial_data['choice']
        reward_time = trial_data['reward_time'].to_numpy()[:,None]
        stays = choice[:-1].to_numpy() == choice[1:].to_numpy()
        switches = np.insert(~stays, 0, False)
        stays = np.insert(stays, 0, False)
        prev_rewarded = np.insert(rewarded[:-1], 0, False)
        prev_unrewarded = np.insert(~rewarded[:-1], 0, False)
        
        resp_rewarded = rewarded[responded]
        
        for region in regions:

            region_side = implant_info[subj_id][region]['side']
            choice_side = choice.apply(lambda x: fpah.get_implant_rel_side(x, region_side) if not x == 'none' else 'none').to_numpy()
            
            # save reward times
            for rew_bin in rew_hist_bins:
                rew_sel = rew_hist == rew_bin
                bin_str = rew_hist_bin_strs[rew_bin]
                
                stack_mat(reward_times[region], 'rew_hist_'+bin_str, reward_time[rew_sel & responded])
                stack_mat(reward_times[region], 'rew_hist_'+bin_str+'_rewarded', reward_time[rew_sel & responded & rewarded])
                stack_mat(reward_times[region], 'rew_hist_'+bin_str+'_unrewarded', reward_time[rew_sel & responded & ~rewarded])
                
                for side in sides:
                    side_sel = choice_side == side
                    stack_mat(reward_times[region], 'rew_hist_'+bin_str+'_'+side, reward_time[rew_sel & responded & side_sel,:])
                    stack_mat(reward_times[region], 'rew_hist_'+bin_str+'_rewarded_'+side, reward_time[rew_sel & responded & rewarded & side_sel,:])
                    stack_mat(reward_times[region], 'rew_hist_'+bin_str+'_unrewarded_'+side, reward_time[rew_sel & responded & ~rewarded & side_sel,:])
            

            for signal_type in signal_types:
                if not signal_type in aligned_signals[subj_id][sess_id]:
                    continue
                for align in alignments:
                    if not align in aligned_signals[subj_id][sess_id][signal_type]:
                        continue

                    t_r = t[align][region]
                    mat = aligned_signals[subj_id][sess_id][signal_type][align][region]

                    # calculate peak properties on a trial-by-trial basis
                    if analyze_peaks:
                        contra_choices = choice_side == 'contra'
                        contra_choices = contra_choices[responded]
                        
                        resp_idx = 0
                        for i in range(mat.shape[0]):
                            if responded[i]:
                                metrics = fpah.calc_peak_properties(mat[i,:], t_r, 
                                                                    filter_params=filter_props[align][region],
                                                                    peak_find_params=peak_find_props[align][region],
                                                                    fit_decay=False)
                                
                                if resp_idx < rew_n_back:
                                    buffer = np.full(rew_n_back-resp_idx, False)
                                    rew_hist_vec = np.flip(np.concatenate((buffer, resp_rewarded[:resp_idx])))
                                    contra_hist_vec = np.flip(np.concatenate((buffer, contra_choices[:resp_idx])))
                                    ipsi_hist_vec = np.flip(np.concatenate((buffer, ~contra_choices[:resp_idx])))
                                else:
                                    rew_hist_vec = np.flip(resp_rewarded[resp_idx-rew_n_back:resp_idx])
                                    contra_hist_vec = np.flip(contra_choices[resp_idx-rew_n_back:resp_idx])
                                    ipsi_hist_vec = np.flip(~contra_choices[resp_idx-rew_n_back:resp_idx])
    
                                peak_metrics.append(dict([('subj_id', subj_id), ('sess_id', sess_id), ('signal_type', signal_type), 
                                                         ('align', align.name), ('region', region), ('trial', i),
                                                         ('rewarded', rewarded[i]), ('rew_hist_bin', rew_hist.iloc[i]), 
                                                         ('side', choice_side[i]), ('rew_hist', rew_hist_vec),
                                                         ('contra_hist', contra_hist_vec), ('ipsi_hist', ipsi_hist_vec),
                                                         ('reward_time', reward_time[i]), ('RT', trial_data['RT'].iloc[i]),
                                                         ('cpoke_out_latency', trial_data['cpoke_out_latency'].iloc[i]), *metrics.items()]))
                                
                                resp_idx += 1
                    
                    if norm_baseline:
                        # normalize all grouped matrices to the pre-event signal of the lowest reward rate
                        baseline_mat = mat[(rew_hist == rew_hist_bins[0]) & responded,:]
                        if baseline_mat.shape[0] > 0:
                            baseline_sel = (t_r >= -0.1) & (t_r < 0)
                            baseline = np.nanmean(baseline_mat[:,baseline_sel])
                        else:
                            baseline = 0
                            
                        mat = mat - baseline
                    
                    # group trials together and stack across sessions
                    match align:
                        case Align.cue:
                            for side in sides:
                                side_sel = choice_side == side
                                stack_mat(stacked_signals[signal_type][align][region], 'stay_prev_reward_'+side, mat[stays & prev_rewarded & side_sel,:])
                                stack_mat(stacked_signals[signal_type][align][region], 'stay_prev_unreward_'+side, mat[stays & prev_unrewarded & side_sel,:])
                                stack_mat(stacked_signals[signal_type][align][region], 'switch_prev_reward_'+side, mat[switches & prev_rewarded & side_sel,:])
                                stack_mat(stacked_signals[signal_type][align][region], 'switch_prev_unreward_'+side, mat[switches & prev_unrewarded & side_sel,:])
                        case Align.reward:
                            for side in sides:
                                side_sel = choice_side == side
                                stack_mat(stacked_signals[signal_type][align][region], 'stay_rewarded_prev_reward_'+side, mat[stays & rewarded & prev_rewarded & side_sel,:])
                                stack_mat(stacked_signals[signal_type][align][region], 'stay_rewarded_prev_unreward_'+side, mat[stays & rewarded & prev_unrewarded & side_sel,:])
                                stack_mat(stacked_signals[signal_type][align][region], 'switch_rewarded_prev_reward_'+side, mat[switches & rewarded & prev_rewarded & side_sel,:])
                                stack_mat(stacked_signals[signal_type][align][region], 'switch_rewarded_prev_unreward_'+side, mat[switches & rewarded & prev_unrewarded & side_sel,:])
                                
                                stack_mat(stacked_signals[signal_type][align][region], 'stay_unrewarded_prev_reward_'+side, mat[stays & ~rewarded & prev_rewarded & side_sel,:])
                                stack_mat(stacked_signals[signal_type][align][region], 'stay_unrewarded_prev_unreward_'+side, mat[stays & ~rewarded & prev_unrewarded & side_sel,:])
                                stack_mat(stacked_signals[signal_type][align][region], 'switch_unrewarded_prev_reward_'+side, mat[switches & ~rewarded & prev_rewarded & side_sel,:])
                                stack_mat(stacked_signals[signal_type][align][region], 'switch_unrewarded_prev_unreward_'+side, mat[switches & ~rewarded & prev_unrewarded & side_sel,:])
                    
                    for rew_bin in rew_hist_bins:
                        rew_sel = rew_hist == rew_bin
                        bin_str = rew_hist_bin_strs[rew_bin]

                        stack_mat(stacked_signals[signal_type][align][region], 'rew_hist_'+bin_str+'_rewarded', mat[rew_sel & responded & rewarded,:])
                        stack_mat(stacked_signals[signal_type][align][region], 'rew_hist_'+bin_str+'_unrewarded', mat[rew_sel & responded & ~rewarded,:])
                        
                        for side in sides:
                            side_sel = choice_side == side
                            stack_mat(stacked_signals[signal_type][align][region], 'rew_hist_'+bin_str+'_rewarded_'+side, mat[rew_sel & responded & rewarded & side_sel,:])
                            stack_mat(stacked_signals[signal_type][align][region], 'rew_hist_'+bin_str+'_unrewarded_'+side, mat[rew_sel & responded & ~rewarded & side_sel,:])

                        match align:
                            case Align.cue:

                                stack_mat(stacked_signals[signal_type][align][region], 'rew_hist_'+bin_str, mat[rew_sel & responded,:])

                                for side in sides:
                                    side_sel = choice_side == side
                                    stack_mat(stacked_signals[signal_type][align][region], 'rew_hist_'+bin_str+'_'+side, mat[rew_sel & responded & side_sel,:])

                            #case Align.reward:
if analyze_peaks:
    peak_metrics = pd.DataFrame(peak_metrics)
    # drop unused columns
    peak_metrics.drop(['decay_tau', 'decay_params', 'decay_form'], axis=1, inplace=True)

# %% declare common plotting stuff

def calc_error(mat, use_se):
    if use_se:
        return utils.stderr(mat, axis=0)
    else:
        return np.nanstd(mat, axis=0, ddof=1)
    
side_labels = {'ipsi': 'Ipsi', 'contra': 'Contra'}
rew_labels = {'rewarded': 'rew', 'unrewarded': 'unrew'}
#bin_labels = {b:'{:.0f}-{:.0f}'.format(np.abs(np.ceil(b.left*rew_rate_n_back)), np.floor(b.right*rew_rate_n_back)) for b in rew_hist_bins}
bin_labels = {b:'{:.0f}'.format(np.mean([np.abs(np.ceil(b.left*rew_rate_n_back)), np.floor(b.right*rew_rate_n_back)])) for b in rew_hist_bins}
group_labels_dict = {'rew_hist_{}_{}'.format(rew_hist_bin_strs[rew_bin], rk): '{} {}'.format(bin_labels[rew_bin], rv)
                for rew_bin in rew_hist_bins
                for rk, rv in rew_labels.items()}
group_labels_dict.update({'rew_hist_{}'.format(rew_hist_bin_strs[rew_bin]): bin_labels[rew_bin]
                for rew_bin in rew_hist_bins})
group_labels_dict.update({'rew_hist_{}_{}'.format(rew_hist_bin_strs[rew_bin], side_type): '{} {}'.format(bin_labels[rew_bin], side_label)
                for rew_bin in rew_hist_bins 
                for side_type, side_label in side_labels.items()})
group_labels_dict.update({'rew_hist_{}_{}_{}'.format(rew_hist_bin_strs[rew_bin], rk, side_type): '{} {} {}'.format(bin_labels[rew_bin], side_label, rv)
                for rew_bin in rew_hist_bins 
                for side_type, side_label in side_labels.items()
                for rk, rv in rew_labels.items()})

parameter_titles = {'peak_time': 'Time to Peak', 'peak_height': 'Peak Amplitude',
                    'peak_width': 'Peak Width', 'decay_tau': 'Decay τ'}
parameter_labels = {'peak_time': 'Time to Peak (s)', 'peak_height': 'Peak Amplitude ({})',
                    'peak_width': 'Peak FWHM (s)', 'decay_tau': 'Decay τ (s)'}
align_labels = {'cue': 'Response Cue', 'reward': 'Reward Delivery'}

# %% Prepare peak metrics

# make subject ids categories
peak_metrics['subj_id'] = peak_metrics['subj_id'].astype('category')
peak_metrics['rew_hist_bin_label'] = peak_metrics['rew_hist_bin'].apply(lambda x: bin_labels[x])
peak_metrics['align_label'] = peak_metrics['align'].apply(lambda x: align_labels[x])

ignore_outliers = True
ignore_any_outliers = True
outlier_thresh = 10

t_min = 0.02
t_max = {a: {r: peak_find_props[a][r]['peak_tmax'] - t_min for r in regions} for a in alignments} 

parameters = ['peak_time', 'peak_height', 'peak_width'] #, 'decay_tau'

filt_peak_metrics = peak_metrics.copy()

# remove outliers on a per-subject basis:
if ignore_outliers:
    
    # first get rid of peaks with times too close to the edges of the peak window (10ms from each edge)
    peak_sel = np.full(len(peak_metrics), False)
    for align in alignments:    
        for region in regions:
            align_region_sel = (peak_metrics['align'] == align) & (peak_metrics['region'] == region)
            sub_peak_metrics = peak_metrics[align_region_sel]
            peak_sel[align_region_sel] = ((sub_peak_metrics['peak_height'] > 0) & 
                                          (sub_peak_metrics['peak_time'] > t_min) &
                                          (sub_peak_metrics['peak_time'] < t_max[align][region]))

    filt_peak_metrics = filt_peak_metrics[peak_sel]
    
    # first add iqr multiple columns
    for param in parameters:
        filt_peak_metrics['iqr_mult_'+param] = np.nan
    
    # calculate iqr multiple for potential outliers
    outlier_grouping = ['subj_id', 'sess_id', 'signal_type']
    
    # compute IQR on different groups of trials based on the alignment and region
    for align in alignments:
        # separate peaks by outcome at time of reward
        if align == Align.reward:
            align_outlier_grouping = outlier_grouping+['rewarded']
        else:
            align_outlier_grouping = outlier_grouping
            
        for region in regions:
            # separate peaks by side for DMS since very sensitive to choice side
            if region == 'DMS':
                region_outlier_grouping = align_outlier_grouping+['side']
            else:
                region_outlier_grouping = align_outlier_grouping
                
            align_region_sel = (filt_peak_metrics['align'] == align) & (filt_peak_metrics['region'] == region)
            
            filt_peak_metrics.loc[align_region_sel, :] = fpah.calc_iqr_multiple(filt_peak_metrics[align_region_sel], region_outlier_grouping, parameters)
    
    # then remove outlier values
    if ignore_any_outliers:

        outlier_sel = np.full(len(filt_peak_metrics), False)
        for param in parameters:
            outlier_sel = outlier_sel | (np.abs(filt_peak_metrics['iqr_mult_'+param]) >= outlier_thresh)
            
        filt_peak_metrics.loc[outlier_sel, parameters] = np.nan
        
    else:
        for param in parameters:
            outlier_sel = np.abs(filt_peak_metrics['iqr_mult_'+param]) >= outlier_thresh

            filt_peak_metrics.loc[outlier_sel, param] = np.nan


# %% Plot average traces across all groups (Fig S9-E)

plot_regions = ['DMS', 'PL'] # 'DMS', 'PL'
plot_aligns = [Align.cue, Align.reward]
plot_signals = ['z_dff_iso']

# plot formatting
plot_dec = {'DMS': 1, 'PL': 2}
x_inc = {'DMS': 0.3, 'PL': 3}
y_inc = {'DMS': 0.2, 'PL': 0.2}

all_color = '#08AB36'
rew_color = '#BC141A'
unrew_color = '#1764AB'
dms_color = '#53C43B'
pl_color = '#BB6ED8'

separate_outcome = False

if separate_outcome:
    gen_groups = {Align.cue: {'rewarded': 'rew_hist_{}_rewarded', 'unrewarded': 'rew_hist_{}_unrewarded'},
                  Align.reward: {'rewarded': 'rew_hist_{}_rewarded', 'unrewarded': 'rew_hist_{}_unrewarded'}}
    group_labels = {'rewarded': 'Rewarded', 'unrewarded': 'Unrewarded'}

    region_colors = {'DMS': {'rewarded': rew_color, 'unrewarded': unrew_color},
                     'PL': {'rewarded': rew_color, 'unrewarded': unrew_color}}
else:
    gen_groups = {Align.cue: {'all': 'rew_hist_{}'},
                  Align.reward: {'all': 'rew_hist_{}_rewarded'}}
    group_labels = {'all': '_'}

    region_colors = {'DMS': {'all': dms_color}, 'PL': {'all': pl_color}}
    
groups = {a: {k: [v.format(rew_hist_bin_strs[rew_bin]) for rew_bin in rew_hist_bins] for k, v in gen_groups[a].items()} for a in plot_aligns}

cue_to_reward = np.nanmedian(sess_data['reward_time'] - sess_data['response_cue_time'])

plot_lims = {Align.cue: {'DMS': [-0.1,0.8], 'PL': [-1,1]},
             Align.reward: {'DMS': [-0.1,1.2], 'PL': [-1,12]}}

width_ratios = [np.diff(plot_lims[align][plot_regions[0]])[0] for align in plot_aligns]
#width_ratios = [2,10.5]
#width_ratios = [0.7,0.9]

n_rows = len(plot_regions)
n_cols = len(plot_aligns)
t = aligned_signals['t']
x_label = 'Time (s)'

for signal_type in plot_signals:
    signal_title, y_label = fpah.get_signal_type_labels(signal_type)
    
    fig, axs = plt.subplots(n_rows, n_cols, layout='constrained', figsize=(5.5, 3*n_rows+0.1), sharey='row', width_ratios=width_ratios)
    
    if n_rows == 1 and n_cols == 1:
        axs = np.array(axs)

    axs = axs.reshape((n_rows, n_cols))
    
    fig.suptitle(signal_title)
    
    plot_name = '{}_all_trials_outcome_{}_{}'.format('_'.join(plot_aligns), signal_type, '_'.join(plot_regions))

    for i, region in enumerate(plot_regions):
        colors = region_colors[region]
        for j, align in enumerate(plot_aligns):
            match align:
                case Align.cue:
                    title = 'Response Cue'
                    
                case Align.reward:
                    title = 'Reward Delivery'

            ax = axs[i,j]
            
            region_signals = stacked_signals[signal_type][align][region]
            t_r = t[align][region]
            
            # normalize to pre-cue levels on a trial-by-trial basis
            cue_signals = stacked_signals[signal_type][Align.cue][region]
            baseline_sel = (t_r >= -0.1) & (t_r < 0)
            
            for key, stack_groups in groups[align].items():
                stacked_signal = np.zeros((0,len(t_r)))
                
                for group in stack_groups:
                    baseline = np.nanmean(cue_signals[group][:,baseline_sel], axis=1)[:,None]
                    stacked_signal = np.vstack((stacked_signal, region_signals[group] - baseline))

                t_sel = (t_r > plot_lims[align][region][0]) & (t_r < plot_lims[align][region][1])
                error = calc_error(stacked_signal, True)
                
                plot_utils.plot_psth(t_r[t_sel][::plot_dec[region]], np.nanmean(stacked_signal, axis=0)[t_sel][::plot_dec[region]], error[t_sel][::plot_dec[region]], 
                                     ax, label=group_labels[key], color=colors[key], plot_x0=False)
    
            ax.set_title('{} {}'.format(region, title))
            plot_utils.plot_dashlines(0, ax=ax)
    
            if j == 0:
                ax.set_ylabel(y_label)
            #else:
                # ax.yaxis.set_tick_params(which='both', labelleft=True)
            ax.legend(loc='upper right')

            ax.set_xlabel(x_label)

            ax.xaxis.set_major_locator(MultipleLocator(x_inc[region]))
            ax.yaxis.set_major_locator(MultipleLocator(y_inc[region]))
            
        fpah.save_fig(fig, fpah.get_figure_save_path('Two-armed Bandit', 'Reward History', plot_name), format='pdf')

# %% Plot average traces by reward history (Fig 5J)

plot_regions = ['DMS', 'PL'] # 
plot_aligns = [Align.reward] # Align.cue, 
plot_signal_types = ['z_dff_iso']

# plot formatting
plot_dec = {'DMS': 1, 'PL': 2}
x_inc = {'DMS': 0.3, 'PL': 3}
y_inc = {'DMS': 0.3, 'PL': 0.3}

#gen_groups = {Align.cue: ['rew_hist_{}_rewarded', 'rew_hist_{}_unrewarded'], Align.reward: ['rew_hist_{}_rewarded', 'rew_hist_{}_unrewarded']}
gen_groups = {Align.cue: ['rew_hist_{}'], Align.reward: ['rew_hist_{}_rewarded', 'rew_hist_{}_unrewarded']}
groups = {a: [group.format(rew_hist_bin_strs[rew_bin]) for group in gen_groups[a] for rew_bin in rew_hist_bins] for a in plot_aligns}

rew_hist_all_colors = plt.cm.Greens(np.linspace(0.4,1,len(rew_hist_bins)))
rew_hist_rew_colors = plt.cm.Reds(np.linspace(0.4,1,len(rew_hist_bins)))
rew_hist_unrew_colors = plt.cm.Blues(np.linspace(0.4,1,len(rew_hist_bins)))

#colors = {Align.cue: np.vstack((rew_hist_rew_colors, rew_hist_unrew_colors)), Align.reward: np.vstack((rew_hist_rew_colors, rew_hist_unrew_colors))}
colors = {Align.cue: rew_hist_all_colors, Align.reward: np.vstack((rew_hist_rew_colors, rew_hist_unrew_colors))}

plot_lims = {Align.cue: {'DMS': [-0.1,0.5], 'PL': [-1,1]},
             Align.reward: {'DMS': [-0.2,1.2], 'PL': [-2,12]}}

#width_ratios = [1,2]
#width_ratios = [0.7,0.9]
if len(plot_aligns) == 1:
    width_ratios = [1]
else:
    width_ratios = [np.diff(plot_lims[Align.cue][plot_regions[0]])[0], np.diff(plot_lims[Align.reward][plot_regions[0]])[0]]

n_rows = len(plot_regions)
n_cols = len(plot_aligns)
t = aligned_signals['t']
x_label = 'Time (s)'

plot_data = {}

for signal_type in plot_signal_types:
    signal_title, y_label = fpah.get_signal_type_labels(signal_type)
    
    plot_data[signal_type] = {}
    
    fig, axs = plt.subplots(n_rows, n_cols, layout='constrained', figsize=(5.5, 3*n_rows+0.1), sharey='row', width_ratios=width_ratios)
    
    if n_rows == 1 and n_cols == 1:
        axs = np.array(axs)

    axs = axs.reshape((n_rows, n_cols))
    
    fig.suptitle(signal_title)
    
    plot_name = 'reward_hist_{}_back_{}_{}_{}'.format(rew_rate_n_back, signal_type, '_'.join(plot_aligns), '_'.join(plot_regions))

    for i, region in enumerate(plot_regions):
        plot_data[signal_type][region] = {}
        
        for j, align in enumerate(plot_aligns):
            match align:
                case Align.cue:
                    title = 'Response Cue'
                    legend_cols = 1
                    
                case Align.reward:
                    title = 'Reward Delivery'
                    legend_cols = 2

            ax = axs[i,j]
            
            plot_data[signal_type][region][align] = {}
            
            region_signals = stacked_signals[signal_type][align][region]
            
            t_r = t[align][region]
            t_sel = (t_r > plot_lims[align][region][0]) & (t_r < plot_lims[align][region][1])
            t_r = t_r[t_sel][::plot_dec[region]]
            
            plot_data[signal_type][region][align]['time'] = t_r
    
            for group, color in zip(groups[align], colors[align]):
                act = region_signals[group]
                error = calc_error(act, True)
                
                avg_act = np.nanmean(act, axis=0)[t_sel][::plot_dec[region]]
                avg_err = error[t_sel][::plot_dec[region]]
                
                plot_data[signal_type][region][align][group_labels_dict[group]+'_avg'] = avg_act
                plot_data[signal_type][region][align][group_labels_dict[group]+'_err'] = avg_err
                
                plot_utils.plot_psth(t_r, avg_act, avg_err, ax, label=group_labels_dict[group], color=color, plot_x0=False)
    
            ax.set_title('{} {}'.format(region, title))
            plot_utils.plot_dashlines(0, ax=ax)
    
            if j == 0:
                ax.set_ylabel(y_label)
            #else:
                # ax.yaxis.set_tick_params(which='both', labelleft=True)
            ax.legend(ncols=legend_cols, loc='upper right', title='# Rewards in last {} Trials'.format(rew_rate_n_back))

            ax.set_xlabel(x_label)
            
            ax.xaxis.set_major_locator(MultipleLocator(x_inc[region]))
            ax.yaxis.set_major_locator(MultipleLocator(y_inc[region]))
            
            plot_data[signal_type][region][align] = pd.DataFrame.from_dict(plot_data[signal_type][region][align])
            
        fpah.save_fig(fig, fpah.get_figure_save_path('Two-armed Bandit', 'Reward History', plot_name), format='pdf')
        

# %% make peak property comparison figures comparing regions at cue and reward (Fig S9-F)

parameters = ['peak_height', 'peak_width'] # 'peak_time', 'decay_tau']

regions = ['DMS', 'PL']
region_colors = ['#53C43B', '#BB6ED8']
plot_aligns = ['Response Cue', 'Reward Delivery']
subj_ids = np.unique(filt_peak_metrics['subj_id'])
plot_signals = ['dff_iso']

# have same jitter for each subject
noise = 0.075
n_neg = int(np.floor(len(subj_ids)/2))
n_pos = int(np.ceil(len(subj_ids)/2))
jitters = np.concatenate([np.random.uniform(-1, -0.1, n_neg), np.random.uniform(0.1, 1, n_pos)]) * noise

for signal_type in plot_signals:
    signal_label, y_label = fpah.get_signal_type_labels(signal_type)
    
    peak_sel = filt_peak_metrics['signal_type'] == signal_type
    # only show rewarded peaks for reward delivery
    peak_sel = peak_sel & ((filt_peak_metrics['align'] == Align.cue) | ((filt_peak_metrics['align'] == Align.reward) & filt_peak_metrics['rewarded']))
    
    sub_peak_metrics = filt_peak_metrics[peak_sel]
    
    # Compare responses across regions grouped by alignment
    fig, axs = plt.subplots(len(parameters), 1, figsize=(4, 4*len(parameters)), layout='constrained')
    axs = np.array(axs).reshape((len(parameters)))
    
    for i, param in enumerate(parameters):
        ax = axs[i]
        ax.set_title(parameter_titles[param])

        subj_avgs = sub_peak_metrics.groupby(['subj_id', 'region', 'align_label']).agg({param: np.nanmean}).reset_index()
    
        # plot sensor averages in boxplots
        sb.boxplot(data=sub_peak_metrics, x='align_label', y=param, hue='region', palette=region_colors,
                   order=plot_aligns, hue_order=regions, ax=ax, showfliers=False, whis=(5,95), log_scale=True)
        if param == 'peak_height':
            ax.set_ylabel(parameter_labels[param].format(y_label))
        else:
            ax.set_ylabel(parameter_labels[param])
        ax.set_xlabel('')
    
        # add subject averages for each alignment with lines connecting them
        dodge = 0.2

        for j, align in enumerate(plot_aligns):
            x = np.array([j - dodge, j + dodge])
    
            for subj_id, jitter in zip(subj_ids, jitters):
                subj_avg = subj_avgs[(subj_avgs['subj_id'] == subj_id) & (subj_avgs['align_label'] == align)]
    
                y = [subj_avg.loc[subj_avg['region'] == r, param] for r in regions]
    
                ax.plot(x+jitter, y, color='black', marker='o', linestyle='solid', alpha=0.75, linewidth=1, markersize=5)
                
        # set y max to be a multiple of 10
        y_min, y_max = ax.get_ylim()
        ax.set_ylim(y_min, utils.convert_to_multiple(y_max, 10, direction='up'))
                
    plot_name = 'cue_reward_peak_comp_{}_{}'.format('_'.join(parameters), signal_type)
    fpah.save_fig(fig, fpah.get_figure_save_path('Two-armed Bandit', 'Reward History', plot_name), format='pdf')

# %% make reward history comparison figures per peak property (Fig 5K)

parameters = ['peak_height'] # 'peak_time', 'peak_height', 'peak_width', 'decay_tau'

rew_hist_rew_colors = plt.cm.Reds(np.linspace(0.4,1,len(rew_hist_bins)))
rew_hist_palette = sb.color_palette(rew_hist_rew_colors) 

plot_signals = ['dff_iso']
plot_aligns = ['reward'] #'cue', 
plot_regions = ['DMS', 'PL']

split_by_side = False
sides = ['contra', 'ipsi']
n_regions = len(regions)

n_aligns = len(plot_aligns)
align_labels = {'cue': 'Response Cue', 'reward': 'Reward Delivery'}
hatch_order = ['//\\\\', '']
line_order = ['dashed', 'solid']

plot_data = {}

for signal_type in plot_signals:
    signal_label, y_label = fpah.get_signal_type_labels(signal_type)
    
    plot_data[signal_type] = {}
    
    for param in parameters:
        
        plot_data[signal_type][param] = {}
        
        plot_name = '{}_reward_hist_{}_{}_back_{}'.format('_'.join(plot_aligns), param, rew_rate_n_back, signal_type)
    
        # Compare responses across alignments grouped by region
        fig, axs = plt.subplots(n_regions, n_aligns, figsize=(4*n_aligns, 4*n_regions), layout='constrained', sharey='row')
        fig.suptitle('{}, {}'.format(parameter_titles[param], signal_label))
        
        axs = np.array(axs).reshape((n_regions, n_aligns))
    
        for i, region in enumerate(plot_regions):
            
            plot_data[signal_type][param][region] = {}
            
            for j, align in enumerate(plot_aligns):
                peak_sel = (filt_peak_metrics['align'] == align) & (filt_peak_metrics['signal_type'] == signal_type) & (filt_peak_metrics['region'] == region)
    
                match align:
                    #case Align.cue:

                    case Align.reward:
                        peak_sel = peak_sel & filt_peak_metrics['rewarded']

                region_metrics = filt_peak_metrics[peak_sel]
                plot_data[signal_type][param][region][align] = {'raw': region_metrics}
                
                ax = axs[i,j]
                ax.set_title('{} {}'.format(region, align_labels[align]))
                
                if split_by_side:
                    # plot reward history group averages in boxplots
                    sb.boxplot(data=region_metrics, x='rew_hist_bin_label', y=param, hue='side', hue_order=sides,
                               whis=(5,95), ax=ax, showfliers=False, legend=False)
                    
                    # update colors and fills of boxes
                    for k, patch in enumerate(ax.patches):
                        # Left boxes first, then right boxes
                        if k < len(rew_hist_bins):
                            # add hatch to cues
                            patch.set_hatch(hatch_order[int(k/len(rew_hist_bins))])

                        patch.set_facecolor(rew_hist_palette[k % len(rew_hist_bins)])

                    # add subject averages for each alignment with lines connecting them
                    subj_avgs = region_metrics.groupby(['subj_id', 'rew_hist_bin_label', 'side']).agg({param: np.nanmean}).reset_index()
                    plot_data[signal_type][param][region][align]['avg'] = subj_avgs
                    
                    group_labels = np.unique(region_metrics['rew_hist_bin_label'])
                    region_subj_ids = np.unique(region_metrics['subj_id'])
                    dodge = [-0.2, 0.2]
                    for k, (side, d) in enumerate(zip(sides, dodge)):
                        x = np.arange(len(group_labels)) + d
                        for subj_id in region_subj_ids:
                            subj_avg = subj_avgs[(subj_avgs['subj_id'] == subj_id) & (subj_avgs['side'] == side)]
                            y = [subj_avg.loc[subj_avg['rew_hist_bin_label'] == g, param] for g in group_labels]
                            ax.plot(x, y, color='black', marker='o', linestyle=line_order[k], alpha=0.7)

                    # Add the custom legend to the figure (or to one of the subplots)
                    legend_patches = [Patch(facecolor='none', edgecolor='black', hatch=hatch_order[k], label=side_labels[s]) for k, s in enumerate(sides)]
                    ax.legend(handles=legend_patches, frameon=False)
                else:
                    # plot reward history group averages in boxplots
                    sb.boxplot(data=region_metrics, x='rew_hist_bin_label', y=param,
                               hue='rew_hist_bin_label', palette=rew_hist_palette, whis=(5,95),
                               ax=ax, showfliers=False, legend=False, saturation=0.7)
    
                    # add subject averages for each alignment with lines connecting them
                    subj_avgs = region_metrics.groupby(['subj_id', 'rew_hist_bin_label']).agg({param: np.nanmean}).reset_index()
                    plot_data[signal_type][param][region][align]['avg'] = subj_avgs
            
                    group_labels = np.unique(region_metrics['rew_hist_bin_label'])
                    region_subj_ids = np.unique(region_metrics['subj_id'])
                    for subj_id in region_subj_ids:
                        subj_avg = subj_avgs[subj_avgs['subj_id'] == subj_id]
                        y = [subj_avg.loc[subj_avg['rew_hist_bin_label'] == g, param] for g in group_labels]
        
                        ax.plot(np.arange(len(group_labels)), y, color='black', marker='o', linestyle='solid', alpha=0.7, linewidth=1, markersize=5)
                        
                if param == 'peak_height':
                    ax.set_ylabel(parameter_labels[param].format(y_label))
                else:
                    ax.set_ylabel(parameter_labels[param])
                ax.set_xlabel('# Rewards in last {} Trials'.format(rew_rate_n_back))
                
        fpah.save_fig(fig, fpah.get_figure_save_path('Two-armed Bandit', 'Reward History', plot_name), format='pdf')

# %% Get mean and SE for peak properties
# by region/align
signal_type = 'dff_iso'
parameters = ['peak_height', 'peak_width'] # 'peak_time', 'decay_tau']
plot_regions = ['DMS', 'PL']
plot_aligns = ['cue', 'reward'] #

for param in parameters:
    for region in plot_regions:
        for align in plot_aligns:
            peak_sel = (filt_peak_metrics['align'] == align) & (filt_peak_metrics['signal_type'] == signal_type) & (filt_peak_metrics['region'] == region)
    
            if align == Align.reward:
                    peak_sel = peak_sel & filt_peak_metrics['rewarded']
    
            region_metrics = filt_peak_metrics[peak_sel]
            
            print('{} {}, {}-aligned {}: {:.3f} +/- {:.3f}\n'.format(
                   signal_type, region, align, param, np.nanmean(region_metrics[param]), utils.stderr(region_metrics[param])))
            
# %% Perform statistical tests on the peak properties

# look at each region and property separately
parameters = ['peak_height'] # ['peak_time', 'peak_height', 'peak_width', 'decay_tau'] #
regions = ['DMS', 'PL']
aligns = [Align.reward] #Align.cue, 
signals = ['dff_iso']
include_side = False
use_bins_as_cats = True # whether to use the bin labels as categories or regressors
print_fit_results = False

for signal_type in signals:
    signal_label, y_label = fpah.get_signal_type_labels(signal_type)
    
    for param in parameters:
        for region in regions:
            for align in aligns:
                
                match align:
                    case Align.cue:
                        region_metrics = filt_peak_metrics[(filt_peak_metrics['region'] == region) & (filt_peak_metrics['align'] == align) 
                                                      & (filt_peak_metrics['signal_type'] == signal_type)]
                        #re_form = '1 + side'
                        #vc = {'subj_id': '0 + C(subj_id)', 'side': '0 + C(side)'}
                    case Align.reward:
                        # only show rewarded peaks heights in reward alignment
                        region_metrics = filt_peak_metrics[(filt_peak_metrics['region'] == region) & (filt_peak_metrics['align'] == align) 
                                                      & (filt_peak_metrics['signal_type'] == signal_type) & filt_peak_metrics['rewarded']]
                        #re_form = '1 + side'
                        #vc = {'subj_id': '0 + C(subj_id)'}
                        
                # drop nans
                region_metrics = region_metrics[~np.isnan(region_metrics[param])]
                
                #NOTE: mixed effects models fit a common variance for a zero-mean gaussian distribution
                # for all random effects where the individual random effect estimates are drawn from that distribution
                # thus adding both a subject and side term to the variance means that separate variances will be fit for those random effects across all groups
                if include_side:
                    # if subj_id is not specified, will fit two means per subject, one for each side which does slightly worse than fitting the subject average in addition to any variability based on side
                    vc_form={'subj_id': '1', 'side': '0 + C(side)'} 
                else:
                    vc_form={'subj_id': '1'} #same as 0 + C(subj_id)
        
                # group_lm = ols(param+' ~ C(group_label)', data=region_metrics).fit()
                # subj_lm = ols(param+' ~ C(subj_id)', data=region_metrics).fit()
                # group_subj_lm = ols(param+' ~ C(group_label)*C(subj_id)', data=region_metrics).fit()
                # subj_group_lm = ols(param+' ~ C(subj_id)*C(group_label)', data=region_metrics).fit()
                
                # # print('{} {}-aligned variant model fit:\n {}\n'.format(region, align, variant_lm.summary()))
                # # print('{} {}-aligned variant & subject model fit:\n {}\n'.format(region, align, variant_subj_lm.summary()))
                
                # print('{}: {} {}, {}-aligned reward history model ANOVA:\n {}\n'.format(param, signal_type, region, align, anova_lm(group_lm)))
                # print('{}: {} {}, {}-aligned subject model ANOVA:\n {}\n'.format(param, signal_type, region, align, anova_lm(subj_lm)))
                # print('{}: {} {}, {}-aligned reward history & subject model ANOVA:\n {}\n'.format(param, signal_type, region, align, anova_lm(group_subj_lm)))
                # print('{}: {} {}, {}-aligned subject & reward history model ANOVA:\n {}\n'.format(param, signal_type, region, align, anova_lm(subj_group_lm)))
                
                # print('{}: {} {}, {}-aligned comparison between reward history & reward history/subject models:\n {}\n'.format(param, signal_type, region, align, anova_lm(group_lm, group_subj_lm)))
                # print('{}: {} {}, {}-aligned comparison between subject & subject/reward history models:\n {}\n'.format(param, signal_type, region, align, anova_lm(subj_lm, subj_group_lm)))
                
                        
                # mem = sm.MixedLM.from_formula(param+' ~ C(group_label)', groups='subj_id', data=region_metrics, missing='drop')
                # print('{}: {} {}, {}-aligned fixed variants, random subjects:\n {}\n'.format(param, signal_type, region, align, mem.fit().summary()))
                
                # Create N mixed effects models where the first reward history group is rotated to compare it with all other groups
                rew_hist_groups = np.unique(region_metrics['rew_hist_bin_label'])
                param_vals = region_metrics[param].to_numpy()
                subj_id_vals = region_metrics['subj_id'].to_numpy()
                side_vals = region_metrics['side'].to_numpy()
                
                if use_bins_as_cats:
                    p_vals = []
                    for i, first_group in enumerate(rew_hist_groups):
                        group_mapping = {first_group: '0'}
                        for other_group in rew_hist_groups[~np.isin(rew_hist_groups, first_group)]:
                            group_mapping.update({other_group: 'bin_'+other_group})
                            
                        rew_hist_vals = region_metrics['rew_hist_bin_label'].apply(lambda x: group_mapping[x])
                        
                        model_data = pd.DataFrame.from_dict({param: param_vals, 'subj_id': subj_id_vals, 'group': rew_hist_vals, 'side': side_vals})
    
                        mem = sm.MixedLM.from_formula(param+' ~ C(group)', groups='subj_id', data=model_data, missing='drop', vc_formula=vc_form).fit() #vc_formula=vc_form, 
                        
                        # get the group comparison p-values
                        p_val_sel = mem.pvalues.index.str.contains('C(group)', regex=False)
                        group_p_vals = mem.pvalues[p_val_sel]
                        group_p_vals.index = first_group+' -> '+group_p_vals.index
                        # only save the unique comparisons
                        p_vals.append(group_p_vals.iloc[i:])

                        if print_fit_results:
                            print('{}: {} {}, {}-aligned Mixed-effects Model, \'{}\' group compared against other groups:\n {}\n'.format(
                                   param, signal_type, region, align, first_group, mem.summary()))
                            print('Random Effects:\n{}\n'.format(mem.random_effects))
                            
                    p_vals = pd.concat(p_vals)
                    reject, _, _, corrected_alpha = smm.multipletests(p_vals, alpha=0.05, method='bonferroni')
                    p_vals = pd.DataFrame(p_vals).rename(columns={0: 'p_val'})
                    p_vals['reject null'] = reject
                    p_vals['corrected alpha'] = corrected_alpha
                    
                    print('{}: {} {}, {}-aligned pairwise group comparison p-values:\n {}\n'.format(
                           param, signal_type, region, align, p_vals))
                else:
                    group_mapping = {g: i for i, g in enumerate(rew_hist_groups)}
                    rew_hist_vals = region_metrics['rew_hist_bin_label'].apply(lambda x: group_mapping[x])
                    model_data = pd.DataFrame.from_dict({param: param_vals, 'subj_id': subj_id_vals, 'group': rew_hist_vals, 'side': side_vals})

                    mem = sm.MixedLM.from_formula(param+' ~ group', groups='subj_id', data=model_data, missing='drop', vc_formula=vc_form).fit() #vc_formula=vc_form, 
                    
                    if print_fit_results:
                        print('{}: {} {}, {}-aligned Mixed-effects Model, linear regression on rew history bin #:\n {}\n'.format(
                               param, signal_type, region, align, mem.summary()))
                        print('Random Effects:\n{}\n'.format(mem.random_effects))
                        
                    print('{}: {} {}, {}-aligned slope regression p-values:\n {}\n'.format(
                           param, signal_type, region, align, mem.pvalues))
                    
                
# %% declare common methods for performing regression over time

# define method to perform the regression
def regress_over_time(signals, predictors):
    params = []
    ci_lower = []
    ci_upper = []
    se = []
    p_vals = []
    rmse = []
    for i in range(signals.shape[1]):
        t_sig = signals[:,i]
        # remove trials with nans at this time step
        rem_nan_sel = ~np.isnan(t_sig)
        
        if np.sum(rem_nan_sel) == 0:
            nan_vals = {col: np.nan for col in predictors.columns}
            params.append(nan_vals)
            ci_lower.append(nan_vals)
            ci_upper.append(nan_vals)
            se.append(nan_vals)
            p_vals.append(nan_vals)
            rmse.append(np.nan)
        else:
            lm = sm.OLS(t_sig[rem_nan_sel], predictors[rem_nan_sel])
            lm = lm.fit()
            #lm = lm.fit_regularized(alpha=0.1, L1_wt=0.5)
            params.append(lm.params.to_dict())
            cis = lm.conf_int(0.05)
            ci_lower.append(cis[0].to_dict())
            ci_upper.append(cis[1].to_dict())
            se.append(lm.bse.to_dict())
            p_vals.append(lm.pvalues.to_dict())
            rmse.append(np.sqrt(np.mean(lm.resid**2)))
        
    return {'params': pd.DataFrame(params), 
            'ci_lower': pd.DataFrame(ci_lower), 
            'ci_upper': pd.DataFrame(ci_upper),
            'se': pd.DataFrame(se),
            'p_vals': pd.DataFrame(p_vals),
            'rmse': np.array(rmse)}

# define common plotting routine
def plot_regress_over_time(params, t, plot_cols, ax, ci_lower=None, ci_upper=None, error=None, sig=None, p_vals=None, group_labels=None,
                           t_sel=None, colors=None, plot_y0=True, plot_x0=True, dec=1, x_inc=None, y_inc=None, filt_outliers=False):
    if len(plot_cols) == 0:
        return
    
    if group_labels is None:
        group_labels = {}
    
    plot_data = {}
    
    sig_y_dist = 0.03
    if t_sel is None:
        t_sel = np.full_like(t, True)
    else:
        t_sel = t_sel.copy()
        
    plot_t = t[t_sel][::dec]
    plot_data['time'] = plot_t
    
    line_colors = []
    for i, col in enumerate(plot_cols):
        vals = params[col].to_numpy()
        plot_vals = vals[t_sel][::dec]
        
        if filt_outliers:
            t_sel = t_sel & (np.abs(vals) < 1e3)
            
        if colors is None:
            color = None
        else:
            color = colors[i]
            
        col_label = group_labels.get(col, col)

        if not ci_lower is None and not ci_upper is None:
            error = np.abs(np.vstack((ci_lower[col], ci_upper[col])) - vals[None,:])
            plot_err = error[:,t_sel][::dec]
            line, _ = plot_utils.plot_psth(plot_t, plot_vals, plot_err, ax=ax, label=col_label, plot_x0=False, color=color)
        elif not error is None:
            plot_err = error[col][t_sel][::dec]
            line, _ = plot_utils.plot_psth(plot_t, plot_vals, plot_err, ax=ax, label=col_label, plot_x0=False, color=color)
        else:
            plot_err = None
            line, _ = plot_utils.plot_psth(plot_t, plot_vals, ax=ax, label=col_label, plot_x0=False, color=color)
            
        plot_data[col_label+'_coeffs'] = plot_vals
        plot_data[col_label+'_err'] = plot_err
        
        line_colors.append(line.get_color())
            
    if plot_x0:
        plot_utils.plot_dashlines(0, dir='v', ax=ax)
    if plot_y0:
        plot_utils.plot_dashlines(0, dir='h', ax=ax)
    
    # plot significance from 0    
    if plot_sig and not sig is None:
        y_min, y_max = ax.get_ylim()
        y_offset = (y_max-y_min)*sig_y_dist

        for i, col in enumerate(plot_cols):
            # # perform correction
            # reject, corrected_pvals, _, _  = smm.multipletests(p_vals[col][t_sel], alpha=0.05, method='fdr_bh')
            sig_t = t[sig[col] & t_sel]
            ax.scatter(sig_t, np.full_like(sig_t, y_max+i*y_offset), color=line_colors[i], marker='.', s=10)
            
            col_label = group_labels.get(col, col)
            plot_data[col_label+'_sig'] = sig[col][t_sel][::dec]
            if not p_vals is None:
                plot_data[col_label+'_p'] = p_vals[col][t_sel][::dec]

    ax.set_xlabel('Time (s)')
    ax.legend(loc='best')
    
    if not x_inc is None:
        ax.xaxis.set_major_locator(MultipleLocator(x_inc))
    if not y_inc is None:
        ax.yaxis.set_major_locator(MultipleLocator(y_inc))
        
    return plot_data
    
    
# %% Prep for reward history regression over time

regions = ['DMS', 'PL'] # 'DMS', 'PL'
aligns = [Align.cue, Align.reward] # Align.cue, Align.reward
signals = ['z_dff_iso']
included_subjs = np.array(subj_ids)

n_back = 3
normalize = True
exclude_n_back = True # whether to exclude trials less than the n_back history

#bah.calc_rew_rate_hist(sess_data, n_back=n_back, kernel='uniform')

# first build predictor and response matrices
t = aligned_signals['t']
    
subj_stacked_signals = {subj_id: {s: {a: {r: np.zeros((0, len(t[a][r]))) for r in regions} 
                                      for a in aligns} 
                                  for s in signals}
                        for subj_id in included_subjs}

subj_predictors = {subj_id: {r: [] for r in regions} for subj_id in included_subjs}

for subj_id in included_subjs:
    for sess_id in sess_ids[subj_id]:
        
        # build predictor matrix by trial
        trial_data = sess_data[sess_data['sessid'] == sess_id]
        responded = ~np.isnan(trial_data['response_time']).to_numpy()
        rewarded = trial_data['rewarded'].to_numpy()[responded].astype(int)
        unrewarded = (~trial_data['rewarded'].to_numpy())[responded].astype(int)
        choice = trial_data['choice'].to_numpy()[responded]
        left_choice = (choice == 'left').astype(int)
        right_choice = (choice == 'right').astype(int)
        switches = np.concatenate(([False], choice[:-1] != choice[1:])).astype(int)
        stays = np.concatenate(([False], choice[:-1] == choice[1:])).astype(int)
        next_stays = np.concatenate((choice[:-1] == choice[1:], [False])).astype(int)
        rew_hist = pd.cut(trial_data['rew_rate_hist_all'], rew_hist_bins)
        rew_time = trial_data['reward_time'].to_numpy()[responded]

        # make buffered predictors to be able to go n back
        # Note: need to add 0 at the end for generic indexing logic (i:-n_back+i-1) - have to index to at least -1, can't do 0
        if not exclude_n_back:
            buff_reward = np.concatenate((np.full((n_back), 0), rewarded, [0]))
            buff_unreward = np.concatenate((np.full((n_back), 0), unrewarded, [0]))
            buff_left_choice = np.concatenate((np.full((n_back), 0), left_choice, [0]))
            buff_right_choice = np.concatenate((np.full((n_back), 0), right_choice, [0]))
        else:
            buff_reward = np.concatenate((rewarded, [0]))
            buff_unreward = np.concatenate((unrewarded, [0]))
            buff_left_choice = np.concatenate((left_choice, [0]))
            buff_right_choice = np.concatenate((right_choice, [0]))
            choice = choice[n_back:]
            switches = switches[n_back:]
            stays = stays[n_back:]
            next_stays = next_stays[n_back:]
            rew_time = rew_time[n_back:]
            
        for region in regions:
            
            # build predictors by region
            region_side = implant_info[subj_id][region]['side']
            choice_side = [fpah.get_implant_rel_side(x, region_side) for x in choice]
            
            preds = {'reward ({})'.format(i-n_back): buff_reward[i:-n_back+i-1] for i in range(n_back, -1, -1)}
            preds.update({'unreward ({})'.format(i-n_back): buff_unreward[i:-n_back+i-1] for i in range(n_back, -1, -1)})
            if region_side == 'left':
                preds.update({'ipsi choice ({})'.format(i-n_back): buff_left_choice[i:-n_back+i-1] for i in range(n_back-1, -1, -1)})
                preds.update({'contra choice ({})'.format(i-n_back): buff_right_choice[i:-n_back+i-1] for i in range(n_back-1, -1, -1)})
            else:
                preds.update({'ipsi choice ({})'.format(i-n_back): buff_right_choice[i:-n_back+i-1] for i in range(n_back-1, -1, -1)})
                preds.update({'contra choice ({})'.format(i-n_back): buff_left_choice[i:-n_back+i-1] for i in range(n_back-1, -1, -1)})
            preds.update({'choice': choice_side})
            preds.update({'switch': switches})
            preds.update({'stay': switches})
            preds.update({'next_stay': switches})
            preds.update({'reward_time': rew_time})
                
            subj_predictors[subj_id][region].append(pd.DataFrame(preds))

            for signal_type in signals:
                if not signal_type in aligned_signals[subj_id][sess_id]:
                    continue
                for align in aligns:
                    if not align in aligned_signals[subj_id][sess_id][signal_type]:
                        continue
            
                    t_r = t[align][region]
                    mat = aligned_signals[subj_id][sess_id][signal_type][align][region]
                    
                    # normalize all grouped matrices to the average pre-event signal of the lowest reward rate
                    if normalize:
                        baseline_mat = mat[(rew_hist == rew_hist_bins[0]) & responded,:]
                        #baseline_mat = mat[responded,:]
                        if baseline_mat.shape[0] > 0:
                            baseline_sel = (t_r >= -0.1) & (t_r < 0)
                            baseline = np.nanmean(baseline_mat[:,baseline_sel])
                        else:
                            baseline = 0
                            
                        mat = mat - baseline
                        
                    mat = mat[responded,:]
                    if exclude_n_back:
                        mat = mat[n_back:,:]
                    
                    subj_stacked_signals[subj_id][signal_type][align][region] = np.vstack((subj_stacked_signals[subj_id][signal_type][align][region], mat))
    
    for region in regions:
        subj_predictors[subj_id][region] = pd.concat(subj_predictors[subj_id][region])
    
# for analyzing meta subject, create new subject entry with everything stacked
subj_predictors['all'] = {r: [] for r in regions}
subj_stacked_signals['all'] = {s: {a: {r: [] for r in regions} 
                                   for a in aligns} 
                               for s in signals}
for region in regions:
    subj_predictors['all'][region] = pd.concat([subj_predictors[subj_id][region] for subj_id in included_subjs])
    
    for signal_type in signals:
        for align in aligns:
            subj_stacked_signals['all'][signal_type][align][region] = np.vstack([subj_stacked_signals[subj_id][signal_type][align][region] for subj_id in included_subjs])



# %% perform regression with various options

limit_rewarded_trials = True
include_current_reward = False # only relevant if not including outcome interaction, mostly for cue-related alignments
include_outcome_reward_interaction = True # Separate regressions based on outcome
fit_ind_subj = False
fit_meta_subj = True

analyzed_subjs = []
if fit_ind_subj:
    analyzed_subjs.extend(included_subjs.tolist())

if fit_meta_subj:
    analyzed_subjs.append('all')
    
reg_params = {subj_id: {s: {r: {a: {} for a in aligns} 
                            for r in regions} 
                        for s in signals}
              for subj_id in analyzed_subjs}

for subj_id in analyzed_subjs:
    for region in regions:
        # build the predictor matrix based on the current options
        region_preds = subj_predictors[subj_id][region]
        pred_mat = {}

        contra_choice = region_preds['choice'] == 'contra'
        ipsi_choice = region_preds['choice'] == 'ipsi'
        full_rewarded = region_preds['reward (0)'].astype(bool)
        rewarded = full_rewarded
        unrewarded = region_preds['unreward (0)'].astype(bool)
        
        if limit_rewarded_trials:
            region_preds = region_preds[rewarded]
            contra_choice = contra_choice[rewarded]
            rewarded = rewarded[rewarded]

        # determine how to model the intercept and current reward
        elif include_outcome_reward_interaction:
            pred_mat['rewarded'] = rewarded.astype(int)
            if not limit_rewarded_trials:
                pred_mat['unrewarded'] = unrewarded.astype(int)
        else:
            pred_mat['intercept'] = 1
            if include_current_reward and not limit_rewarded_trials:
                pred_mat['rewarded'] = rewarded.astype(int)

        # add in reward history, don't go all the way to the current trial
        for i in range(n_back-1, -1, -1):
            rew_str = 'reward ({})'.format(i-n_back)
            rew_preds = region_preds[rew_str]
            
            if include_outcome_reward_interaction:
                pred_mat['rewarded, '+rew_str] = rew_preds * rewarded
                if not limit_rewarded_trials:
                    pred_mat['unrewarded, '+rew_str] = rew_preds * unrewarded
            else:
                pred_mat[rew_str] = rew_preds

        pred_mat = pd.DataFrame(pred_mat)

        for signal_type in signals:
            for align in aligns:
                
                signal_mat = subj_stacked_signals[subj_id][signal_type][align][region]
                if limit_rewarded_trials:
                    signal_mat = signal_mat[full_rewarded,:]
                 
                reg_params[subj_id][signal_type][region][align] = regress_over_time(signal_mat, pred_mat)
                
                
# %% plot regression coefficients over time (Fig 5L)

# Create the reward history colormap
cmap = LinearSegmentedColormap.from_list('red_to_blue', [plt.cm.Reds(0.7), plt.cm.Blues(0.7)])

plot_signals = ['z_dff_iso']
plot_regions = ['DMS', 'PL'] # 'DMS', 'PL'
plot_aligns = [Align.reward] #  Align.cue, Align.reward

# plot formatting
plot_dec = {'DMS': 1, 'PL': 2}
x_inc = {'DMS': 0.3, 'PL': 3}
y_inc = {'DMS': 0.3, 'PL': 0.3}
plot_lims = {Align.cue: {'DMS': [-0.1,0.6], 'PL': [-1,2]},
             Align.reward: {'DMS': [-0.2,1.2], 'PL': [-2,12]}}

plot_n_back = n_back # 

plot_ind_subj = False
plot_subj_average = False
plot_meta_subj = True
plot_current_reward_separate = False
plot_rmse = False
use_ci_errors = False
plot_sig = True

# Process p-values for significance
sig_lvl = 0.05
method = 'bonferroni' # 'bonferroni' 'holm'

if plot_sig:
    for subj_id in analyzed_subjs:
        for region in plot_regions:
            for signal_type in plot_signals:
                for align in plot_aligns:
    
                    p_vals = reg_params[subj_id][signal_type][region][align]['p_vals']
                    sig = p_vals.apply(lambda x: smm.multipletests(x, alpha=sig_lvl, method=method)[0], axis=1, result_type='broadcast').astype(bool)
                    reg_params[subj_id][signal_type][region][align]['sig'] = sig

# Build plot groups based on regression options

# Simplify side/reward interaction labels only if not plotting them separately
if not plot_current_reward_separate:
    group_labels = {'contra, rewarded': 'rewarded', 'contra, unrewarded': 'unrewarded', 'ipsi, rewarded': 'rewarded', 'ipsi, unrewarded': 'unrewarded'}
else:
    group_labels = {}
    
rew_hist_label = 'reward ({})'
group_labels.update({g.format(rew_hist_label.format(i-plot_n_back)): rew_hist_label.format(i-plot_n_back) for i in range(plot_n_back-1, -1, -1) for g in 
                     ['{}', 'contra, {}', 'ipsi, {}', 'rewarded, {}', 'unrewarded, {}', 'contra, rewarded, {}', 'contra, unrewarded, {}', 'ipsi, rewarded, {}', 'ipsi, unrewarded, {}']})

# reward history groups
    
if include_outcome_reward_interaction:
    plot_group_labels = ['Reward History before Rewarded Choices']
    hist_groups = ['rewarded, reward ({})']
    
    if not limit_rewarded_trials:
        plot_group_labels.append('Reward History before Unrewarded Choices')
        hist_groups.append('unrewarded, reward ({})')
    
    if plot_current_reward_separate:
        plot_groups = [[g.format(i-plot_n_back) for i in range(plot_n_back-1, -1, -1)] for g in hist_groups]
    else:
        rew_groups = [['rewarded']]
        if not limit_rewarded_trials:
            rew_groups.append(['unrewarded'])
        plot_groups = [r+[g.format(i-plot_n_back) for i in range(plot_n_back-1, -1, -1)] for r, g in zip(rew_groups, hist_groups)]
else:
    plot_group_labels = ['Reward History']
    if plot_current_reward_separate:
        plot_groups = [['reward ({})'.format(i-plot_n_back) for i in range(plot_n_back-1, -1, -1)]]
    else:
        if include_current_reward:
            plot_groups = [['rewarded']+['reward ({})'.format(i-plot_n_back) for i in range(plot_n_back-1, -1, -1)]]
        else:
            plot_groups = [['reward ({})'.format(i-plot_n_back) for i in range(plot_n_back-1, -1, -1)]]
    
group_colors = [cmap(np.linspace(0,1,len(g))) for g in plot_groups]

# group all other parameters into one plot    
plot_group_labels.append('Other Parameters')
other_group = []

if include_outcome_reward_interaction:
    if plot_current_reward_separate:
        other_group.append('rewarded')
        if not limit_rewarded_trials:
            other_group.append('unrewarded')
else:
    other_group.append('intercept')
    if plot_current_reward_separate and include_current_reward:
        other_group.append('rewarded')
    
if any(other_group):    
    plot_groups.append(other_group)
    group_colors.append(['C{}'.format(i) for i, _ in enumerate(other_group)])

width_ratios = [np.diff(plot_lims[align]['DMS'])[0] for align in plot_aligns]
    
n_rows = len(plot_regions)
n_cols = len(plot_aligns)

plot_data = {}

for signal_type in plot_signals:
    signal_label, y_label = fpah.get_signal_type_labels(signal_type)
    
    plot_data[signal_type] = {}
    
    if plot_ind_subj:
        for plot_group, group_label, colors in zip(plot_groups, plot_group_labels, group_colors):
            for subj_id in included_subjs:
                fig, axs = plt.subplots(n_rows, n_cols, layout='constrained', figsize=(5.5, 3*n_rows+0.1), width_ratios=width_ratios, sharey='row')
                axs = np.array(axs).reshape((n_rows, n_cols))
                
                fig.suptitle('{} Regression, Subj {}'.format(group_label, subj_id))
                
                for i, region in enumerate(plot_regions):
                    for j, align in enumerate(plot_aligns):
                        ax = axs[i,j]
                        
                        t_r = t[align][region]
                        subj_params = reg_params[subj_id][signal_type][region][align]
                        t_sel = (t_r > plot_lims[align][region][0]) & (t_r < plot_lims[align][region][1])
                        
                        if use_ci_errors:
                            plot_regress_over_time(subj_params['params'], t_r, plot_group, ax, region,
                                                   ci_lower=subj_params['ci_lower'], ci_upper=subj_params['ci_upper'], 
                                                   sig=subj_params['sig'], t_sel=t_sel, colors=colors,
                                                   dec=plot_dec[region], x_inc=x_inc[region], y_inc=y_inc[region])
                        else:
                            plot_regress_over_time(subj_params['params'], t_r, plot_group, ax, region,
                                                   error=subj_params['se'], sig=subj_params['sig'], 
                                                   t_sel=t_sel, colors=colors, dec=plot_dec[region], x_inc=x_inc[region], y_inc=y_inc[region])   
                        
                        ax.set_title('{}, {}-aligned'.format(region, align))
                        
                        if j == 0:
                            ax.set_ylabel('Coefficient ({})'.format(y_label))
                            
        if plot_rmse:
            for subj_id in included_subjs:
                fig, axs = plt.subplots(n_rows, n_cols, layout='constrained', figsize=(5.5, 3*n_rows+0.1), width_ratios=width_ratios, sharey='row')
                axs = np.array(axs).reshape((n_rows, n_cols))
                fig.suptitle('Regression RMSE, Subj {}'.format(subj_id))
                
                for i, region in enumerate(plot_regions):
                    for j, align in enumerate(plot_aligns):
                        ax = axs[i,j]
                        t_r = t[align][region]
                        
                        rmse = reg_params[subj_id][signal_type][region][align]['rmse']
                        t_sel = (t_r > plot_lims[align][region][0]) & (t_r < plot_lims[align][region][1])
                        plot_regress_over_time(pd.DataFrame({'rmse': rmse}), t_r, ['rmse'], ax, region, t_sel=t_sel, plot_y0=False)
                        
                        ax.set_title('{}, {}-aligned'.format(region, align))
                        
                        if j == 0:
                            ax.set_ylabel('RMSE ({})'.format(y_label))
                        
    if plot_subj_average:
        for plot_group, group_label, colors in zip(plot_groups, plot_group_labels, group_colors):
            fig, axs = plt.subplots(n_rows, n_cols, layout='constrained', figsize=(5.5, 3*n_rows+0.1), width_ratios=width_ratios, sharey='row')
            axs = np.array(axs).reshape((n_rows, n_cols))
            fig.suptitle('{} Regression, Subject Avg'.format(group_label))
            
            for i, region in enumerate(plot_regions):
                for j, align in enumerate(plot_aligns):
                    ax = axs[i,j]
                    t_r = t[align][region]
                    
                    # average coefficients across subjects
                    all_params = pd.concat([reg_params[subj_id][signal_type][region][align]['params'] for subj_id in included_subjs])
                    param_avg = all_params.groupby(level=0).mean()
                    param_se = all_params.groupby(level=0).std() / np.sqrt(len(included_subjs))
                    
                    t_sel = (t_r > plot_lims[align][region][0]) & (t_r < plot_lims[align][region][1])
            
                    plot_regress_over_time(param_avg, t_r, plot_group, ax, region,
                                           error=param_se, t_sel=t_sel, colors=colors,
                                           dec=plot_dec[region], x_inc=x_inc[region], y_inc=y_inc[region])
                    
                    ax.set_title('{}, {}-aligned'.format(region, align))
                    
                    if j == 0:
                        ax.set_ylabel('Coefficient ({})'.format(y_label))
                        
    if plot_meta_subj:
        for plot_group, group_label, colors in zip(plot_groups, plot_group_labels, group_colors):
            subj_id = 'all'
            fig, axs = plt.subplots(n_rows, n_cols, layout='constrained', figsize=(5.5, 3*n_rows+0.1), width_ratios=width_ratios, sharey='row')
            axs = np.array(axs).reshape((n_rows, n_cols))
            
            fig.suptitle('{} Regression, Subj {}'.format(group_label, subj_id))
            
            for i, region in enumerate(plot_regions):
                if not region in plot_data[signal_type]:
                    plot_data[signal_type][region] = {}
                    
                for j, align in enumerate(plot_aligns):
                    ax = axs[i,j]
                    t_r = t[align][region]
                    
                    subj_params = reg_params[subj_id][signal_type][region][align]
                    t_sel = (t_r > plot_lims[align][region][0]) & (t_r < plot_lims[align][region][1])
                    
                    if not align in plot_data[signal_type][region]:
                        plot_data[signal_type][region][align] = {}
                    
                    if use_ci_errors:
                        group_plot_data = plot_regress_over_time(subj_params['params'], t_r, plot_group, ax, region,
                                               ci_lower=subj_params['ci_lower'], ci_upper=subj_params['ci_upper'], group_labels=group_labels,
                                               sig=subj_params['sig'], p_vals=subj_params['p_vals'], t_sel=t_sel, colors=colors,
                                               dec=plot_dec[region], x_inc=x_inc[region], y_inc=y_inc[region])
                    else:
                        group_plot_data = plot_regress_over_time(subj_params['params'], t_r, plot_group, ax, region,
                                               error=subj_params['se'], sig=subj_params['sig'], p_vals=subj_params['p_vals'], group_labels=group_labels,
                                               t_sel=t_sel, colors=colors, dec=plot_dec[region], x_inc=x_inc[region], y_inc=y_inc[region])
                    
                    plot_data[signal_type][region][align].update(group_plot_data)
                    
                    ax.set_title('{}, {}-aligned'.format(region, align))
                    
                    if j == 0:
                        ax.set_ylabel('Coefficient ({})'.format(y_label))
                        
            plot_name = '{}_{}_time_regression_{}_{}'.format('_'.join(plot_regions), '_'.join(plot_aligns), group_label, signal_type)
            fpah.save_fig(fig, fpah.get_figure_save_path('Two-armed Bandit', 'Reward History', plot_name), format='pdf')
        
        for region in plot_regions:
            for align in plot_aligns:
                plot_data[signal_type][region][align] = pd.DataFrame.from_dict(plot_data[signal_type][region][align])
        
        if plot_rmse:
            fig, axs = plt.subplots(n_rows, n_cols, layout='constrained', figsize=(5.5, 3*n_rows+0.1), width_ratios=width_ratios, sharey='row')
            axs = np.array(axs).reshape((n_rows, n_cols))
            
            fig.suptitle('Regression RMSE, Subj {}'.format(subj_id))
            
            for i, region in enumerate(plot_regions):
                for j, align in enumerate(plot_aligns):
                    ax = axs[i,j]
                    t_r = t[align][region]
                    
                    rmse = reg_params[subj_id][signal_type][region][align]['rmse']
                    t_sel = (t_r > plot_lims[align][region][0]) & (t_r < plot_lims[align][region][1])
                    plot_regress_over_time(pd.DataFrame({'rmse': rmse}), t_r, ['rmse'], ax, region, t_sel=t_sel, plot_y0=False)
                    
                    ax.set_title('{}, {}-aligned'.format(region, align))
                    
                    if j == 0:
                        ax.set_ylabel('RMSE ({})'.format(y_label))
            
            