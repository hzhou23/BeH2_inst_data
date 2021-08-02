#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import uncertainties
from uncertainties import ufloat
import string
from scipy.optimize import curve_fit

##~~~~~~~~ Input ~~~~~~~~~~~~
pd.options.display.float_format = '{:,.6f}'.format
folders = ['uccsd-t']	# Methods in increasing accuracy
basis = np.array([2,3,4,5,6])
card1=3		# First cardinal to use in extrapolation
card2=6		# Second cardinal to use in extrapolation
best_method='uccsd-t'	# Best method to use in extrapolation
##~~~~~~~~~~~~~~~~~~~~~~~~~~~

def scf_cbs(n, e_cbs, c, d):	#SCF
        y = e_cbs + c*np.exp(-d*n)
        return y

def corr_cbs(n, e_cbs, c, d):	#Correlation
        y = e_cbs + c/(n+3.0/8.0)**3.0 + d/(n+3.0/8.0)**5.0
        return y

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

scf = pd.DataFrame(columns=basis, index=folders)
corr = pd.DataFrame(columns=basis, index=folders)
extr = pd.DataFrame(columns=["CBS"], index=folders)

best_corr=0.0
for count, i in enumerate(folders):
	for j in basis:
		mypath = os.path.join(i, str(j) + ".csv")
		try:
			df = pd.read_csv(mypath, sep='\s*,\s*', engine='python')	# \s*,\s* gets rid of empty spaces in column names
			scf[j].loc[i]=df['SCF'].iloc[0]
			corr[j].loc[i]=df['CCSD'].iloc[0] - df['SCF'].iloc[0]
		except:
			if i == best_method:
				print('Estimating energy at {} method and {} basis:'.format(i, j))
				scf[j].loc[i] = scf[j].loc[folders[count-1]]
				eps = corr[j-1].loc[i] - corr[j-1].loc[folders[count-1]]
				#print("EPS: ", eps)
				corr[j].loc[i] = corr[j].loc[folders[count-1]] + eps
			else:
				scf[j].loc[i]=np.nan
				corr[j].loc[i]=np.nan

	corr_x = corr.loc[:, card1:card2]
	try:
		#print("Method:",i)
		ydata = list(corr_x.loc[i].values)
		ydata = np.array(ydata, dtype='float64')
		xdata = list(corr_x.columns)
		#print(xdata,ydata)
		
		initial=[min(ydata), 1.0, 1.0]
		limit=( [-np.inf, -np.inf, -np.inf], [min(ydata), np.inf, np.inf] )
		popt_corr, pcov_corr = curve_fit(corr_cbs, xdata, ydata, p0=initial, bounds=limit)
		std_corr = np.sqrt(np.diag(pcov_corr))[0]
		#print("Fit params:", popt_corr, std_corr, "\n")
		if str(std_corr) == "inf":
			my_cbs = popt_corr[0]
			extr['CBS'].loc[i] = my_cbs
		else:
			my_cbs = ufloat(popt_corr[0], np.sqrt(np.diag(pcov_corr))[0])
			extr['CBS'].loc[i] = frmtr.format("{0:.2u}", my_cbs)
		if i == best_method:
			best_corr = my_cbs

		#x=np.linspace(xdata[0],xdata[-1],50)
		#fig, ax = plt.subplots()
		#plt.plot(xdata, ydata, 'o')
		#plt.plot(x, corr_cbs(x, *popt_corr), '--', lw=1, label="Fitted Function")
		#plt.axhline(y=popt_corr[0], ls='dashed', lw=1, color='red', label='CBS limit')
		#ax.set(title=' %s Extrapolation Fit' % i,
		#xlabel='Cardinal n',
		#ylabel='E(hartree)')
		#legend = ax.legend(loc='best', shadow=False)
		###plt.savefig('extrap.png', format='pdf')
		#plt.show()

	except:
		print("Unexpected error at %s:" % i, sys.exc_info()[0])
		#raise

###~~~~~~~~~~~SCF~~~~~~~~~~~~~~~~~~~~~~~~~
scf_x = scf.loc[:, card1:card2]
ydata = list(scf_x.iloc[0].values)
ydata = np.array(ydata,dtype='float64')
xdata = list(scf_x.columns)

initial=[min(ydata), 0.01, 1.0]
limit=( [-np.inf, 0.00, 0.00], [min(ydata), np.inf, np.inf])
popt_scf, pcov_scf = curve_fit(scf_cbs, xdata, ydata, p0=initial, bounds=limit)
std_scf = np.sqrt(np.diag(pcov_corr))[0]
#print("SCF", popt_scf, std_scf)
if str(std_scf) == "inf":
	my_scf = popt_scf[0]
	scf['CBS'] = my_scf
else:
	my_scf = ufloat(popt_scf[0], std_scf)
	scf['CBS'] = frmtr.format("{0:.2u}", my_scf)

#x=np.linspace(xdata[0],xdata[-1],50)
#fig, ax = plt.subplots()
#plt.plot(xdata, ydata, 'o')
#plt.plot(x, scf_cbs(x, *popt_scf), '--', lw=1, label="Fitted Function")
#plt.axhline(y=popt_scf[0], ls='dashed', lw=1, color='red', label='CBS limit')
#ax.set(title='SCF Extrapolation Fit',
#xlabel='Cardinal n',
#ylabel='E(hartree)')
#legend = ax.legend(loc='best', shadow=False)
###plt.savefig('extrap_scf.png', format='pdf')
#plt.show()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#print(scf.to_latex(na_rep=""))
#print(corr.to_latex(na_rep=""))

scf = scf.iloc[[0]]
scf = scf.rename(index={"{}".format(folders[0]):"ROHF"})

corr=pd.concat([corr, extr], axis=1)
corr=pd.concat([corr, scf], axis=0)
my_total = my_scf + best_corr
if isinstance(my_total, float) == True:
	corr.loc['Total','CBS'] = my_total
else:
	corr.loc['Total','CBS'] = frmtr.format("{0:.2u}", my_total)

scf  =  scf.rename(index={"rcisd":"CISD","rccsd-t":"RCCSD(T)","uccsd-t":"UCCSD(T)", "ccsdt-q":"CCSDT(Q)", "fci":"FCI"}, columns={2:"DZ", 3:"TZ", 4:"QZ", 5:"5Z", 6:"6Z"})
corr = corr.rename(index={"rcisd":"CISD","rccsd-t":"RCCSD(T)","uccsd-t":"UCCSD(T)", "ccsdt-q":"CCSDT(Q)", "fci":"FCI"}, columns={2:"DZ", 3:"TZ", 4:"QZ", 5:"5Z", 6:"6Z"})

#scf = scf.transpose()
#corr = corr.transpose()

print(corr.to_latex(na_rep=""))

