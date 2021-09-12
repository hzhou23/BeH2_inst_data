#!/usr/bin/python3

import pandas as pd
import numpy as np
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import uncertainties
from uncertainties import ufloat
from uncertainties import unumpy
from uncertainties import ufloat_fromstr
from scipy.optimize import curve_fit
import re
import string

styles = {
'singlenoT' :{'label':' ccECP-reg/AE single-ref ', 'color':'#ff9933', 'linestyle':'--', 'dashes': (3,1), 'marker':'s', 'markersize':6, 'markeredgecolor':'#000000', 'markeredgewidth':0.5},
'multinoT'   :{'label':'AE minimal active space multi-ref',    'color':'#009933', 'linestyle':'--', 'dashes': (6,1), 'marker':'o', 'markersize':6, 'markeredgecolor':'#000000', 'markeredgewidth':0.5},
}

def init():
    font   = {'family' : 'serif', 'size': 18}
    lines  = {'linewidth': 2.0}
    axes   = {'linewidth': 2.0}    # border width
    tick   = {'major.size': 2, 'major.width':2}
    legend = {'frameon':False, 'fontsize':17.0, 'handlelength':2.30, 'labelspacing':0.40, 'handletextpad':0.4, 'loc':'best'}

    mpl.rc('font',**font)
    mpl.rc('lines',**lines)
    mpl.rc('axes',**axes)
    mpl.rc('xtick',**tick)
    mpl.rc('ytick',**tick)
    mpl.rc('legend',**legend)

    mpl.rcParams['text.usetex'] = True
    mpl.rcParams.update({'figure.autolayout':True})
    fig = plt.figure()
    fig.set_size_inches(7.04, 5.28)   # Default 6.4, 4.8
    ax1 = fig.add_subplot(111) # row, column, nth plot
    #gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
    #ax1 = plt.subplot(gs[0])
    #ax2 = plt.subplot(gs[1], sharex = ax1)
    return fig, ax1

###################################################################

dmc_data = pd.read_csv('dmc.csv',sep = ' ',names = ['energy','errorbar'])
timestep = [0.001,0.0015,0.002,0.0025,0.003,0.0035,0.004,0.0045,0.005,0.0055,0.006,0.0065,0.007,0.0075,0.008,0.0085,0.009,0.0095,0.01,0.02,0.035,0.045,0.05,0.055,0.07,0.1]
energy, errorbar = dmc_data['energy'][2:], dmc_data['errorbar'][2:]
fig, ax1 = init()
ax1.tick_params(direction='in', length=8, width=2.0, which='major', pad=10)
ax1.tick_params(direction='in', length=4, width=1.1, which='minor', pad=10)
ax1.set_ylabel('Energy [Ha]')
ax1.set_xlabel('Timestep [Ha$^{-1}$]')
ax1.set_xscale("log")

ax1.errorbar(timestep,energy,yerr = errorbar, **styles['singlenoT'])
#labels = [item.get_text() for item in ax1.get_yticklabels()]
#labels = ['0','-15.9022','-15.9020','-15.9018','-15.9016','-15.9014','-15.9012']
#ax1.set_yticklabels(labels)
ax1.legend(ncol=1)
ax1.grid(b=None, which='major', axis='both', alpha=0.1)
plt.savefig('aesrnoT.pdf')
plt.show()

