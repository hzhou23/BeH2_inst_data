#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import uncertainties
from uncertainties import ufloat
from uncertainties import unumpy
from uncertainties import ufloat_fromstr
from scipy.optimize import curve_fit
from textwrap import wrap
import re
import string

styles = {
        '2eSRnoTmoves' :{'label':'ccECP/[He] single-ref no T-moves', 'color':'g', 'ecolor':'g', 'linestyle':'--', 'dashes': (3,1), 'marker':'s', 'markersize':5, 'markeredgecolor':'#000000', 'markeredgewidth':0.5},
        '2eSRTmoves'   :{'label':'ccECP/[He] single-ref T-moves',    'color':'tab:orange', 'ecolor':'tab:orange', 'linestyle':'--', 'dashes': (3,1), 'marker':'o', 'markersize':5, 'markeredgecolor':'#000000', 'markeredgewidth':0.5},
        '2eCASnoTmoves' :{'label':'ccECP/[He] restr-CI no T-moves', 'color':'b', 'ecolor':'b', 'linestyle':'--', 'dashes': (3,1), 'marker':'^', 'markersize':5, 'markeredgecolor':'#000000', 'markeredgewidth':0.5},
        '2eCASTmoves' :{'label':'ccECP/[He] restr-CI T-moves', 'color':'r', 'ecolor':'r', 'linestyle':'--', 'dashes': (3,1), 'marker':'^', 'markersize':5, 'markeredgecolor':'#000000', 'markeredgewidth':0.5},
}

def init():
    font   = {'family' : 'serif', 'size': 16}
    lines  = {'linewidth': 2.0}
    axes   = {'linewidth': 2.0}    # border width
    tick   = {'major.size': 2.0, 'major.width':2}
    legend = {'frameon':False, 'fontsize':15, 'handlelength':2.30, 'labelspacing':0.40, 'handletextpad':0.4, 'loc':'best'}

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

###  ~~~~~~~ Functions ~~~~~~~~~

class ShorthandFormatter(string.Formatter):
    def format_field(self, value, format_spec):
        if isinstance(value, uncertainties.UFloat):
            return value.format(format_spec+'S')  # Shorthand option added
        # Special formatting for other types can be added here (floats, etc.)
        else:
            # Usual formatting:
            return super(ShorthandFormatter, self).format_field(
                value, format_spec)
frmtr = ShorthandFormatter()

def f2s(x):   # float to string
    try:
        y = frmtr.format("{0:.1u}", x)
    except:
        y = x
    return y

def s2f(x):   # string to float
    try:
        y = ufloat_fromstr(x)
    except:
        y = x
    return y

def pd_s2f(df):   # df to df with uncertainties
    udf = df.copy()
    for column in df.columns:
        try:
            udf[column] = list(map(s2f,list(df[column])))
        except:
            udf[column] = "pd_s2f failed"
    return udf

def pd_f2s(udf):   # df to df with uncertainties
    df = udf.copy()
    for column in udf.columns:
        try:
            df[column] = list(map(f2s,list(udf[column])))
        except:
            df[column] = "pd_f2s failed"
    return df

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ccECP_2e_TM = pd.read_csv('2eccECPTQMC_new.csv')
ccECP_2e_TM = pd_s2f(ccECP_2e_TM)
print(ccECP_2e_TM)

fig, ax1 = init()
ax1.tick_params(direction='in', length=8, width=2.0, which='major', pad=10)
ax1.tick_params(direction='in', length=4, width=1.1, which='minor', pad=10)
ax1.set(ylabel= 'Energy [Ha]',xlabel='Number of steps')

timestep = ccECP_2e_TM['walltime']
energy = unumpy.nominal_values(np.array(ccECP_2e_TM['energy'].values))
errors = unumpy.std_devs(np.array(ccECP_2e_TM['energy'].values))

ax1.errorbar(timestep, energy, yerr=errors, **styles['2eSRTmoves'])

ccECP_2e_TM = pd.read_csv('2eccECPnoTQMC_new.csv')
ccECP_2e_TM = pd_s2f(ccECP_2e_TM)
print(ccECP_2e_TM)

timestep = ccECP_2e_TM['walltime']
energy = unumpy.nominal_values(np.array(ccECP_2e_TM['energy'].values))
errors = unumpy.std_devs(np.array(ccECP_2e_TM['energy'].values))

ax1.errorbar(timestep, energy, yerr=errors, **styles['2eSRnoTmoves'])

ccECP_2e_TM = pd.read_csv('2eccECPCASQMC_new.csv')
ccECP_2e_TM = pd_s2f(ccECP_2e_TM)
print(ccECP_2e_TM)

timestep = ccECP_2e_TM['walltime']
energy = unumpy.nominal_values(np.array(ccECP_2e_TM['energy'].values))
errors = unumpy.std_devs(np.array(ccECP_2e_TM['energy'].values))

ax1.errorbar(timestep, energy, yerr=errors, **styles['2eCASnoTmoves'])

ccECP_2e_TM = pd.read_csv('2eccECPCASTQMC_new.csv')
ccECP_2e_TM = pd_s2f(ccECP_2e_TM)
print(ccECP_2e_TM)

timestep = ccECP_2e_TM['walltime']
energy = unumpy.nominal_values(np.array(ccECP_2e_TM['energy'].values))
errors = unumpy.std_devs(np.array(ccECP_2e_TM['energy'].values))

ax1.errorbar(timestep, energy, yerr=errors, **styles['2eCASTmoves'])



ax1.legend(ncol=1)
ax1.grid(b=None, which='major', axis='both', alpha=0.1)



plt.savefig('2eccECPTQMCall_new.pdf')
plt.show()
