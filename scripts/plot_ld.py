
import pickle
import sys

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



votu = 'vOTU-000010'

type_ = 'all'

ld_counts_dict = pickle.load(open(calculate_linkage_disequilibria.ld_counts_dict_path % votu, "rb"))

distances = list(ld_counts_dict['data'][type_].keys())
distances.sort()

ld = [ld_counts_dict['data'][type_][d]['rsquared_numerators']/ld_counts_dict['data'][type_][d]['rsquared_denominators'] for d in distances]

distances = numpy.asarray(distances)
ld = numpy.asarray(ld)

to_keep_idx = (ld>0) & (~numpy.isnan(ld))



fig, ax = plt.subplots(figsize=(4,4))


ax.plot(distances[to_keep_idx], ld[to_keep_idx], lw=2, ls='-', c='dodgerblue', label='Data')

predict = diversity_utils.predict_qle_ld(0.01, distances, 1.4)

ax.plot(distances,predict, lw=2, ls='-', c='k', label='Eq. S6; Garud, Good et al.')



ax.set_xlim([min(distances[to_keep_idx]), max(distances[to_keep_idx])])
ax.set_ylim([0.01, 1])


ax.set_xscale('log', base=10)
ax.set_yscale('log', base=10)


ax.set_xlabel('Pairwise distance, (bp)')
ax.set_ylabel('LD')

ax.set_title(votu, fontsize=14)

ax.legend(loc='lower left')


fig.subplots_adjust(hspace=0.15, wspace=0.15)
fig_name = "%sld_test.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()

#print(distances[to_keep_idx])