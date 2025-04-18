
import numpy
import matplotlib.pyplot as plt
import scipy.stats as stats
import collections
import data_utils
import config
import plot_utils
import stats_utils


fig = plt.figure(figsize = (4, 4))
fig.subplots_adjust(bottom= 0.15)

ax = plt.subplot2grid((1, 1), (0, 0), colspan=1)

rgb_red = plot_utils.make_colormap(len(data_utils.checkv_quality_all))


for q_idx, q in enumerate(data_utils.checkv_quality_all[::-1]):

    votu_dict, vgenome_dict = data_utils.read_uhgv_metadata(checkv_quality=q)

    # convert to counts
    votus_all = []
    for key in votu_dict.keys():

        #if len(phage_metadata_dict[key]['sra']) != len(set(phage_metadata_dict[key]['sra'])):
        #    print(key)

        votus_all.extend([key]*len(set(votu_dict[key]['sra_sample'])))


    votus_counts_all = numpy.asarray(list(collections.Counter(votus_all).values()))

    votus_counts_range = range(1, max(votus_counts_all)+1, 5)
    survival_counts = stats_utils.make_survival_dist(votus_counts_all, votus_counts_range, probability=False)

    
    ax.plot(votus_counts_range, survival_counts, lw=2, ls='-', c=rgb_red(q_idx), label='Quality = ' + q)



ax.legend(loc="upper right", fontsize=8)

ax.set_xscale('log', base=10)
ax.set_yscale('log', base=10)

ax.set_xlabel('Min. sample size', fontsize=12)
ax.set_ylabel('# vOTUs', fontsize=12)


fig.subplots_adjust(hspace=0.3, wspace=0.25)
fig_name = "%ssample_size_survival.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.3, dpi = 600)
plt.close()