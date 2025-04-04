import numpy


from matplotlib import colors
from matplotlib import cm
import matplotlib as mpl

lifestyle_color_dict = {'both':'k', 'temperate':'#87CEEB', 'lytic':'#FF6347'}



def rand_jitter(arr):

    '''
    Adds jitter to datapoints for boxplots, etc
    '''

    stdev = 0.01 * (max(arr) - min(arr))
    return arr + numpy.random.randn(len(arr)) * stdev



def get_latex_pvalue(p_value):

    if p_value < 0.05:
        label = r'$P < 0.05$'
    else:
        label = r'$P \nleq 0.05$'

    return label


def make_colormap(n_entries):
    # some matplotlib tools
    # taxonomic hierarchy colors
    cmap_offset = int(0.2*16)
    # +cmap_offset
    rgb_red_ = cm.Blues(numpy.linspace(0,1,n_entries+5))
    rgb_red_ = mpl.colors.ListedColormap(rgb_red_[cmap_offset:,:-1])

    return rgb_red_
