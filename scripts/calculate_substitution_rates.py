import config
import os.path
import sys

#import diversity_utils
#import gene_diversity_utils

import stats_utils
import data_utils
from math import log10,ceil
import numpy

#import core_gene_utils
import gzip
import os
import pickle
import glob


substitution_rate_directory = '%ssubstitution_rates/' % (config.data_directory)
intermediate_filename_template = '%s%s.txt.gz'  

#min_coverage = config.min_median_coverage
min_sample_size = 10


if __name__=='__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
    parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
    parser.add_argument("votu", help="Name of specific vOTU to run code on")
    args = parser.parse_args()
    
    debug = args.debug
    chunk_size = args.chunk_size
    species_name=args.votu
    good_species_list = [species_name]

    os.system('mkdir -p %s' % substitution_rate_directory)
    
    
    sys.stderr.write("Loading genome annotation map....\n")
    annotation_dict = pickle.load(open(data_utils.annotation_dict_path, "rb"))
    
    for species_name in good_species_list:

        #sys.stderr.write("Loading haploid samples...\n")
        #sys.stderr.write("Proceeding with %d haploid samples!\n" % len(snp_samples))
        
        sys.stderr.write("Loading aligned FASTA for vOTU %s...\n" % species_name)
        
        fasta_file_path_all = glob.glob('%sSingle_vOTUs/Core_Alignments/*%s*.fna' % (config.data_directory, species_name))
        
        if len(fasta_file_path_all) == 0:
            continue
        else:
            fasta_file_path = fasta_file_path_all[0]
            
        
        fasta_all_genomes = data_utils.classFASTA(fasta_file_path).readFASTA()
        fasta_genome_dict = {x[0]:x[1] for x in fasta_all_genomes}
        
        
        
        
        
        
        
        
        
        