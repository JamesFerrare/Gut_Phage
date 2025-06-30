import copy

import os
import json
import config
import data_utils
import numpy
import sys
import pickle
from collections import Counter


import random
#from collections import Counter
from itertools import combinations
import matplotlib.pyplot as plt
from scipy import special


n_allele_prob_dist_dict_path = config.data_directory + 'n_allele_prob_dist_dict.pkl'


def make_n_allele_prob_dist_dict(max_fraction_nan = 0.05):

    votu_all = data_utils.get_single_votus()
    n_alleles_count = {}
    n_alleles_count['probability'] = {}
    
    for i in range(1,5):
        n_alleles_count['probability'][i] = []

    votu_all_final = []
    for votu in votu_all:
            
        if votu in data_utils.votu_to_skip:
            continue
        
        allele_counts_map_path_ = data_utils.allele_counts_map_path % votu
        
        allele_counts_map = pickle.load(open(allele_counts_map_path_, "rb"))
        
        sites =  list(allele_counts_map['aligned_sites'].keys())
        
        genomes = allele_counts_map['genomes']
        genome_pairs_idx = list(combinations(range(len(genomes)), 2))
        n_genomes = len(genomes)
        
        # get sites with not too many NaNs
        n_alleles_all = []
        for s in sites:
            
            alleles_s = allele_counts_map['aligned_sites'][s]['alleles']
            fraction_nan = alleles_s.count('-')/n_genomes
            
            # insufficient number of informative sites
            if fraction_nan > max_fraction_nan:
                continue 
        
            allele_count_dict = dict(Counter(alleles_s))
            # ignore sites with > 2 alleles        
            nucleotide_intersect = set(allele_count_dict.keys()) & set(data_utils.nucleotides)
            
            n_alleles_all.append(len(nucleotide_intersect))
            
            
        allele_count_dict = dict(Counter(n_alleles_all))
        
        keys = list(allele_count_dict.keys())
        keys.sort()
        
        total_n = sum(allele_count_dict.values())
                
        for k in keys:
            n_alleles_count['probability'][k].append(allele_count_dict[k]/total_n)
        
        
        votu_all_final.append(votu)
        
        
    n_alleles_count['votu_all'] = votu_all_final
    
    with open(n_allele_prob_dist_dict_path, 'wb') as f:
        pickle.dump(n_alleles_count, f)



def plot_fourfold_dist():
    
    fourfold_dist_dict = pickle.load(open(n_allele_prob_dist_dict_path, "rb"))

    fig, ax = plt.subplots(figsize=(4,4))
    
    
    x = []
    y = []
    for k in fourfold_dist_dict['prob_fourfold'].keys():
        
        x.extend([k] * len(fourfold_dist_dict['prob_fourfold'][k]))
        y.extend(fourfold_dist_dict['prob_fourfold'][k])
        
    
    ax.scatter(x, y, s=10, alpha=0.3)
    ax.set_yscale('log', base=10)
    
    ax.set_xlabel('Number alleles', fontsize=8)
    ax.set_ylabel('Probability density', fontsize=8)


    fig.subplots_adjust(hspace=0.15, wspace=0.15)
    fig_name = "%sfourfold_dist.png" % (config.analysis_directory)
    fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
    plt.close()
        
    
    




plot_fourfold_dist()
