# -*- coding: utf-8 -*-
"""
Helper functions

@author: tanner stevenson
"""

import numpy as np
import os
import os.path as path
import matplotlib.colors as mc
import colorsys
from scipy.stats import binomtest
import warnings
from math import gcd


# %% Number utils
def is_scalar(x):
    return np.ndim(x) == 0

def convert_to_multiple(value, factor, direction='nearest'):

    match direction:
        case 'nearest':
            return factor * np.round(value/factor)
        case 'up':
            return factor * np.ceil(value/factor)
        case 'down':
            return factor * np.floor(value/factor)


def to_int(x):
    return np.round(x).astype(int)

def nancount(x, axis=0):
    return np.sum(~np.isnan(x), axis=axis)

def stderr(x, axis=0, ignore_nan=True):
    ''' Calculates standard error '''

    x = np.array(x)
    # ignore warnings from all nans in x
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)

        if ignore_nan:
            std = np.nanstd(x, axis=axis, ddof=1)
            n = nancount(x, axis=axis)
        else:
            std = np.std(x, axis=axis, ddof=1)
            n = np.shape(x)[axis]

        se = std/np.sqrt(n)

        # handle cases where n = 0
        if is_list(se):
            se[np.isinf(se)] = np.nan
        elif np.isinf(se):
            se = np.nan

    return se

def z_score(x):
    ''' Calculates z_score over all values in x '''

    # Flatten input to calculate z-score over all values in x
    x = np.array(x)
    flat_x = x.flatten()
    # ignore warnings from all nans in x
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', category=RuntimeWarning)

        flat_z = (flat_x - np.nanmean(flat_x))/np.nanstd(flat_x)

    # reshape the z scores to match original x
    return np.reshape(flat_z, x.shape)

def rescale(x, new_min, new_max, axis=None):
    # if an axis is specified, need to create the axis array to compute the min/max over
    if axis is None:
        ax_arr = None
    else:
        ax_arr = np.arange(len(x.shape)).astype(int)
        ax_arr = tuple(ax_arr[ax_arr != axis])
        
    old_min = np.broadcast_to(np.min(x, axis=ax_arr, keepdims=True), x.shape)
    old_max = np.broadcast_to(np.max(x, axis=ax_arr, keepdims=True), x.shape)
    return new_min + (x - old_min) * (new_max - new_min) / (old_max - old_min)

def binom_cis(k, n, level=0.95):
    # handle non-integer values
    if k % 1 != 0 or n % 1 != 0:
        raise ValueError('k and n must be integers. k={0}, n={1}'.format(k, n))
    else:
        k = int(k)
        n = int(n)
    if n > 0:
        return np.array(binomtest(k, n).proportion_ci(confidence_level=level)[:])
    else:
        return np.array([np.nan, np.nan])

def get_sequence_lengths(x):
    ''' Get a list of lengths of repeated values in array indexed by value '''
    unique_vals = np.unique(x)
    lengths = {v: [] for v in unique_vals}

    seq_len = 1
    for i in range(len(x)-1):
        if x[i] == x[i+1]:
            seq_len += 1

        if x[i] != x[i+1] or i == len(x)-2:
            lengths[x[i]].append(seq_len)
            seq_len = 1

    return lengths

def lcm(x, y):
    # get the least common multiple of two numbers
    return x * y // gcd(x, y)

# %% OS utils
def check_make_dir(full_path):
    path_dir = os.path.dirname(full_path)
    if not os.path.exists(path_dir):
        os.makedirs(path_dir)

def get_user_home():
    return path.expanduser('~')

# %% Object utils

def flatten(lists, recursive=True):
    ''' flattens a structure of lists '''
    if recursive:
        if is_dict(lists):
            return flatten([i for v in lists.values() for i in v], recursive)
        elif is_list(lists) and len(lists) > 0 and is_list(lists[0]):
            return flatten([i for v in lists for i in v], recursive)
        else:
            return lists
    else:
        if is_dict(lists):
            return [i for v in lists.values() for i in v]
        elif is_list(lists) and is_list(lists[0]):
            return [i for v in lists for i in v]
        else:
            return lists


def is_dict(value):
    return isinstance(value, dict)

def is_list(value):
    return isinstance(value, list) or isinstance(value, tuple) or isinstance(value, np.ndarray)

def get_all_combinations(*items):
    ''' Gets all possible combinations of the elements in the provided 1 dimensional item lists '''
    
    # Create a meshgrid from the input arrays
    mesh = np.meshgrid(*items, indexing='ij')
    
    # Stack and reshape to get a 2D array of combinations
    combos = np.stack(mesh, axis=-1).reshape(-1, len(items))
    
    return combos

# %% Plot utils

def change_color_lightness(color, mult):
    """
    Lightens or darkens the given color by changing the the current lightness by the given multiplier (from -1 to 1).
    To lighten the color, the multiplier should be positive, to darken it should be negative. 
    Input can be matplotlib color string, hex string, or RGB tuple.

    From: https://gist.github.com/technic/8bf3932ad7539b762a05da11c0093ed5

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """

    try:
        c = mc.cnames[color]
    except:
        c = color
    c = np.array(colorsys.rgb_to_hls(*mc.to_rgb(c)))
    
    if mult > 0: # lighten, send closer to 1
        new_l = c[1] + mult*(1-c[1])
    else: # darken, send closer to 0
        new_l = c[1] + mult*c[1]
        
    return colorsys.hls_to_rgb(c[0], new_l, c[2])
