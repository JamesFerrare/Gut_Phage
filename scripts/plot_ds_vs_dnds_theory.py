
import copy

import os
import json
import config
import data_utils
import numpy
import sys
import pickle

import random
#from collections import Counter
from itertools import combinations
import matplotlib.pyplot as plt
from scipy import special


s = 0.001
f_d = 0.99

def ds_vs_dnds(ds, s, f_d, mu=10**-8):

    x = (s*ds)/(2*mu)

    return (1-f_d) + f_d * ((1 - numpy.exp(-1*x) )/x)



def ds_vs_dnds_approx(ds, s, f_d, mu=10**-8):
    
    x = (s*ds)/(2*mu)

    

    return (1-f_d) + f_d * (1 - x/special.factorial(2) + (x**2)/special.factorial(3) - (x**3)/special.factorial(4))




fig, ax = plt.subplots(figsize=(4,4))



ds = numpy.logspace(-2.6, -1, base=10, num=1000)


#dnds = ds_vs_dnds(ds, s, f_d)
#dnds_approx = ds_vs_dnds_approx(ds, s, f_d)


#ax.plot(ds, ds_vs_dnds(ds, 0.001, f_d), c='b')
#ax.plot(ds, ds_vs_dnds(ds, 0.00001, f_d), c='b')
#ax.plot(ds, ds_vs_dnds(ds, 0.0000001, 0.1), c='b')
#ax.plot(ds, ds_vs_dnds(ds, 0.0000001, 0.5), c='b')
ax.plot(ds, ds_vs_dnds(ds, 0.00001, 0.99), c='b')

#ax.plot(ds, dnds_approx, c='r')

ax.set_xscale('log', base=10)
ax.set_yscale('log', base=10)


fig.subplots_adjust(hspace=0.15, wspace=0.15)
fig_name = "%sds_vs_dnds_theory.png" % (config.analysis_directory)
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
