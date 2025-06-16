import numpy
import sys
import os

import pickle
from itertools import combinations
from collections import Counter
import data_utils
import config


ld_counts_dict_path = config.data_directory + 'ld_counts_dict_all/%s.pkl'


def build_ld_counts_dict(votu, max_fraction_nan=0.05, max_d=1e3):
    
    allele_counts_map_path_ = data_utils.allele_counts_map_path % votu
    allele_counts_map = pickle.load(open(allele_counts_map_path_, "rb"))
    
    sites =  list(allele_counts_map['aligned_sites'].keys())
    
    genomes = allele_counts_map['genomes']
    genome_pairs_idx = list(combinations(range(len(genomes)), 2))
    n_genomes = len(genomes)
    
    # get sites with not too many NaNs
    sites_final = []
    for s in sites:
        
        #alleles = numpy.asarray(allele_counts_map['aligned_sites'][s]['alleles'])
        #fourfold_status = numpy.asarray(allele_counts_map['aligned_sites'][s]['fourfold_status'])
        #fraction_nan = sum(alleles=='-')/len(alleles)
        alleles_s = allele_counts_map['aligned_sites'][s]['alleles']
        fraction_nan = alleles_s.count('-')/n_genomes
        
        # insufficient number of informative sites
        if fraction_nan > max_fraction_nan:
            continue 
       
        allele_count_dict = dict(Counter(alleles_s))
        # ignore sites with > 2 alleles        
        nucleotide_intersect = set(allele_count_dict.keys()) & set(data_utils.nucleotides)
        if len(nucleotide_intersect) > 2:
            continue
                
        # define major allele
        major_allele = max(list(nucleotide_intersect), key=lambda k: allele_count_dict[k])
        #major_allele = list(nucleotide_intersect - set(minor_allele))[0]
        allele_counts_map['aligned_sites'][s]['major_allele'] = major_allele
        allele_counts_map['aligned_sites'][s]['n_alleles'] = len(nucleotide_intersect)
        allele_counts_map['aligned_sites'][s]['n_obs_no_nan'] = sum(allele_count_dict.values())
        
        # make numpy arrays
        no_nan_bool_idx = numpy.asarray([x != '-' for x in  alleles_s])
        # True if = minor allele or '-'
        allele_bool_idx = numpy.asarray([x != major_allele for x in  alleles_s])
        
        allele_counts_map['aligned_sites'][s]['no_nan_bool_idx'] = no_nan_bool_idx
        allele_counts_map['aligned_sites'][s]['allele_bool_idx'] = allele_bool_idx
        
        fourfold_status = allele_counts_map['aligned_sites'][s]['fourfold_status']
        allele_counts_map['aligned_sites'][s]['fourfold_status'] = numpy.asarray(fourfold_status)
        
        sites_final.append(s)
    
    
    
    ld_count_dict = {}
    ld_count_dict['data'] = {}
    ld_count_dict['genomes'] = genomes
    for variant_type in data_utils.variant_types + ['all']:
        ld_count_dict['data'][variant_type] = {}
        ld_count_dict['data'][variant_type]['site_pairs'] = []
        #ld_count_dict['data'][variant_type]['ns'] = []
        ld_count_dict['data'][variant_type]['n11s'] = []
        ld_count_dict['data'][variant_type]['n10s'] = []
        ld_count_dict['data'][variant_type]['n01s'] = []
        ld_count_dict['data'][variant_type]['n00s'] = []
        
        
    #sites_final = sites_final[:20]
    #sites_final_pairs = list(combinations(sites_final, 2))
    print(len(sites_final))
    
    #for site_pair_idx, site_pair in enumerate(sites_final_pairs):
    n_pairs_processed = 0
    for site_1_idx in range(len(sites_final)):
        
        for site_2_idx in range(site_1_idx):
                
            #s_1 = site_pair[0]
            #s_2 = site_pair[1]
            
            s_1 = sites_final[site_1_idx]
            s_2 = sites_final[site_2_idx]
            
            site_pair = (s_1, s_2)
            
            if abs(s_1-s_2) >= max_d:
                continue
            
            if n_pairs_processed % 1000000 == 0:                
                sys.stderr.write("%d site pairs processed...\n" % n_pairs_processed)  
            
            # genomes with nucleotides in both sites
            no_nan_bool_idx_1 = allele_counts_map['aligned_sites'][s_1]['no_nan_bool_idx']
            no_nan_bool_idx_2 = allele_counts_map['aligned_sites'][s_2]['no_nan_bool_idx']
            no_nan_bool_idx_inter = no_nan_bool_idx_1 * no_nan_bool_idx_2
            
            #n = sum(no_nan_bool_idx_inter)
            
            allele_bool_idx_1 = allele_counts_map['aligned_sites'][s_1]['allele_bool_idx']
            allele_bool_idx_2 = allele_counts_map['aligned_sites'][s_2]['allele_bool_idx']

            allele_bool_idx_final_1 = allele_bool_idx_1[no_nan_bool_idx_inter]
            allele_bool_idx_final_2 = allele_bool_idx_2[no_nan_bool_idx_inter]
            
            n11 = sum(allele_bool_idx_final_1*allele_bool_idx_final_2)
            n10 = sum(allele_bool_idx_final_1*(~allele_bool_idx_final_2))
            n01 = sum((~allele_bool_idx_final_1)*allele_bool_idx_final_2)
            n00 = sum((~allele_bool_idx_final_1)*(~allele_bool_idx_final_2))
            
            
            ld_count_dict['data']['all']['site_pairs'].append(site_pair)
            #ld_count_dict['data']['all']['ns'].append(n)
            ld_count_dict['data']['all']['n11s'].append(n11)
            ld_count_dict['data']['all']['n10s'].append(n10)
            ld_count_dict['data']['all']['n01s'].append(n01)
            ld_count_dict['data']['all']['n00s'].append(n00)
        
            fourfold_status_1 = allele_counts_map['aligned_sites'][s_1]['fourfold_status']
            fourfold_status_2 = allele_counts_map['aligned_sites'][s_2]['fourfold_status']
            
            for variant_type in data_utils.variant_types:
                # check whether both sites have trhe same fourfold status in each genome
                variant_type_idx = (fourfold_status_1 == variant_type) & (fourfold_status_2 == variant_type)
                
                # no genomes that have the same variant type in both sites
                if sum(variant_type_idx) == 0:
                    continue
                
                #n_var = sum(no_nan_bool_idx_inter*variant_type_idx)
                allele_bool_var_idx_1 = allele_bool_idx_1[no_nan_bool_idx_inter*variant_type_idx]
                allele_bool_var_idx_2 = allele_bool_idx_2[no_nan_bool_idx_inter*variant_type_idx]
                
                n11_var = sum(allele_bool_var_idx_1*allele_bool_var_idx_2)
                n10_var = sum(allele_bool_var_idx_1*(~allele_bool_var_idx_2))
                n01_var = sum((~allele_bool_var_idx_1)*allele_bool_var_idx_2)
                n00_var = sum((~allele_bool_var_idx_1)*(~allele_bool_var_idx_2))

                ld_count_dict['data'][variant_type]['site_pairs'].append(site_pair)
                #ld_count_dict['data'][variant_type]['ns'].append(n_var)
                ld_count_dict['data'][variant_type]['n11s'].append(n11_var)
                ld_count_dict['data'][variant_type]['n10s'].append(n10_var)
                ld_count_dict['data'][variant_type]['n01s'].append(n01_var)
                ld_count_dict['data'][variant_type]['n00s'].append(n00_var)
                
                
            n_pairs_processed += 1  

    
    ld_counts_dict_path_ = ld_counts_dict_path % votu
    with open(ld_counts_dict_path_, 'wb') as f:
        pickle.dump(ld_count_dict, f)

    





if __name__ == "__main__":

    votu = 'vOTU-000010'
    #build_allele_counts_map(votu)
    
    build_ld_counts_dict(votu)
    
    #annotation_dict = pickle.load(open(ld_counts_dict_path % votu, "rb"))
    
