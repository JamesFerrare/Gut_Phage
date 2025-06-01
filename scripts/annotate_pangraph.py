import copy
import os
import json
import config
import numpy
import sys
import pickle
import gzip

import random
#from collections import Counter
from itertools import combinations

random.seed(123456789)

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import pypangraph
import data_utils




# get fasta fules
#votu = 'vOTU-000007'

#data_utils.parse_annotation_table()

votu_all = data_utils.get_single_votus()
start_idx = votu_all.index('vOTU-007481')

for votu in votu_all:

    if votu in data_utils.votu_to_skip:
        continue
    
    #data_utils.build_votu_fasta(votu, build_votu_fasta=True)
    #data_utils.make_syn_sites_votu_dict_from_pangraph(votu)


#    continue
#    






#uhgv_votu_metadata_dict, uhgv_genome_metadata_dict = data_utils.read_uhgv_metadata(checkv_quality=checkv_quality, checkv_quality_cumulative=True)



