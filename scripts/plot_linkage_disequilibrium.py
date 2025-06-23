
import pickle
import sys
import os

import numpy
import scipy.stats as stats
import scipy.spatial as spatial

import collections
import data_utils
import config

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors


import calculate_linkage_disequilibria
import diversity_utils
import plot_divergence
import plot_utils



def plot_linkage_disequilibrium_decay_alignment(votu, min_n_site_pairs=80, reference_distance=9):
    
    
    cumulative_n_syn_all, cumulative_n_nonsyn_all, cumulative_block_len_syn_all, cumulative_block_len_nonsyn_all, genome_pairs_clean_all = plot_divergence.calculate_syn_div_and_nonsyn_ratio_alignment(votu, min_n_muts=50, min_n_sites=1e3, check_metadata=True)

    ds_all = cumulative_n_syn_all/cumulative_block_len_syn_all
    pi_4d = numpy.median(ds_all)
    
    ld_counts_dict_path = calculate_linkage_disequilibria.ld_counts_dict_path % votu
    
    if os.path.exists(ld_counts_dict_path) == False:
        
        sys.stderr.write("No LD file for %s! Skipping... \n" % votu)
        
        
    else:
        
        sys.stderr.write("Plotting LD for %s ....\n" % votu)
    
        ld_counts_dict = pickle.load(open(ld_counts_dict_path, "rb"))
        
        fig, ax = plt.subplots(figsize=(4,4))

        for fourfold_status in ['all', 0, 3]:
            
            distances = list(ld_counts_dict['data'][fourfold_status].keys())
            distances.sort()

            rsquareds = [ld_counts_dict['data'][fourfold_status][d]['rsquared_numerators']/ld_counts_dict['data'][fourfold_status][d]['rsquared_denominators'] for d in distances]
            n_site_pairs = [ld_counts_dict['data'][fourfold_status][d]['n_site_pairs'] for d in distances]
            
            #print(n_site_pairs)
            distances = numpy.asarray(distances)
            rsquareds = numpy.asarray(rsquareds)
            n_site_pairs = numpy.asarray(n_site_pairs)
                    
            to_keep_idx = (rsquareds>0) & (~numpy.isnan(rsquareds)) & (n_site_pairs >= min_n_site_pairs)
            
            #if sum(to_keep_idx) < 20:
            #    continue

            ax.plot(distances[to_keep_idx], rsquareds[to_keep_idx], lw=1, alpha=0.7, ls='-', c=plot_utils.fourfold_color_dict[fourfold_status], label=plot_utils.fourfold_label_dict[fourfold_status])

            #predict = diversity_utils.predict_qle_ld(0.01, distances, 1.4)
            
            
            #if fourfold_status == 3:
                
            #    if reference_distance not in distances[to_keep_idx]:
                    
            #        reference_distance_idx = numpy.abs(distances[to_keep_idx] - reference_distance).argmin()
            #        reference_distance = distances[to_keep_idx][reference_distance_idx]

            #    rbymu_2, rbymu_4 = diversity_utils.predict_ld_rbymu(distances[to_keep_idx], rsquareds[to_keep_idx], pi_4d, reference_bp=reference_distance)

            #    print(rbymu_4)
            
            

        #theory_ls = numpy.logspace(0, numpy.log10(max(distances)), 100, base=10)
        #theory_NRs = theory_ls/200.0
        #theory_NRs = theory_ls * pi * rbymu_4
        
        #theory_rsquareds = (10+2*theory_NRs)/(22+26*theory_NRs+4*theory_NRs*theory_NRs)
        #ax.plot(theory_ls, theory_rsquareds, lw=2, ls='-', c='k', label='Eq. S6; Garud, Good et al.')


        #ax.set_xlim([min(distances[to_keep_idx]), max(distances[to_keep_idx])])
        ax.set_xlim([1, 1000])
        ax.set_ylim([0.01, 1])


        ax.set_xscale('log', base=10)
        ax.set_yscale('log', base=10)


        ax.set_xlabel('Distance between SNVs, ' + r'$\ell$', fontsize=12)
        ax.set_ylabel('Linkage disequilibrium, ' + r'$\sigma^{2}_{d}$', fontsize=12)

        ax.set_title(votu, fontsize=14)

        ax.legend(loc='lower left')


        fig.subplots_adjust(hspace=0.15, wspace=0.15)
        fig_name = "%slinkage_disequilibrium_decay_alignment/%s.png" % (config.analysis_directory, votu)
        fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        plt.close()





if __name__ == "__main__":

    votu_all = data_utils.get_single_votus()
    #syn_div_nonsyn_ratio_dip_test()
    #start_idx = votu_all.index('vOTU-000005') + 1 
    #votu_all = votu_all[start_idx:]
    
    #votu = 'vOTU-000010'
    #cumulative_n_syn_all, cumulative_n_nonsyn_all, cumulative_block_len_syn_all, cumulative_block_len_nonsyn_all = calculate_syn_div_and_nonsyn_ratio_alignment(votu)
    
    #calculate_divergence_alignment(votu)
    #plot_ds_vs_dnds_dist_axis(votu)
    
    for votu in votu_all:
        
        plot_linkage_disequilibrium_decay_alignment(votu)
        

