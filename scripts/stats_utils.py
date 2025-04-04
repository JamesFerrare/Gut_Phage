
import numpy
from scipy import stats

import gzip
import config
import os
from math import radians, cos, sin, asin, sqrt


def make_survival_dist(data, range_, probability=True):

    '''
    (Slowly) converts an array into a survival distribution over range range_
    '''

    data = data[numpy.isfinite(data)]
    survival_array = [sum(data>=i) for i in range_]
    survival_array = numpy.asarray(survival_array)

    if probability == True:
        survival_array = survival_array/len(data)

    return survival_array



def permutation_two_sample_ks_test(array_1, array_2, n=1000):

    '''
    Two-sample KS test via permutation
    '''

    statistic, pvalue = stats.ks_2samp(array_1, array_2, alternative='two-sided')

    array_merged = numpy.concatenate((array_1, array_2), axis=None)
    
    statistic_null = []
    for n_i in range(n):

        numpy.random.shuffle(array_merged)
        statistic_null.append(stats.ks_2samp(array_merged[:len(array_1)], array_merged[len(array_1):], alternative='two-sided')[0])

    statistic_null = numpy.asarray(statistic_null)

    p_value = sum(statistic_null > statistic)/n

    return statistic, p_value


def haversine(lon1, lat1, lon2, lat2):
   
    """
    Calculate the great circle distance in kilometers between two points 
    on the earth (specified in decimal degrees)

    Useful for calculating distance betwen two spatial points separated by large distances (e.g., Italy vs. China)
    """
    
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles. Determines return value units.
    
    return c * r



def computed_poisson_thinning(diffs, opportunities):
    
    # apply this to calculation of all dN/dS
    # specifically when calculating dS
    thinned_diffs_1 = numpy.random.binomial(diffs, 0.5)
    thinned_diffs_2 = diffs - thinned_diffs_1
    d1 = thinned_diffs_1 / (opportunities.astype(float) / 2)
    d2 = thinned_diffs_2 / (opportunities.astype(float) / 2)
   
    return d1, d2



