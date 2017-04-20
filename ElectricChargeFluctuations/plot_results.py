#!/usr/bin/env python

from numpy import *
from numpy.core import numeric as _nx
from pylab import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import griddata
from matplotlib.patches import Ellipse
import sys, os

mpl.rcParams['pdf.fonttype'] = 42
labelsize = 20
mpl.rcParams['xtick.labelsize'] = labelsize
mpl.rcParams['ytick.labelsize'] = labelsize

#################################################################
# Parameters to run script with
#################################################################
#file information
#filename = sys.argv[1]

#grid information
nDY = 51

panelLabels = ['(a)', '(b)']
panelCounter = 0

chosenParticle = 'pion'

hbarC = 0.197327053

#################################################################
# Shouldn't have to change anything below this line
#################################################################

def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

'''def stack(arrays, axis=0):
	arrays = [np.asanyarray(arr) for arr in arrays]
	if not arrays:
		raise ValueError('need at least one array to stack')
	
	shapes = set(arr.shape for arr in arrays)
	if len(shapes) != 1:
		raise ValueError('all input arrays must have the same shape')
	
	result_ndim = arrays[0].ndim + 1
	if not -result_ndim <= axis < result_ndim:
		msg = 'axis {0} out of bounds [-{1}, {1})'.format(axis, result_ndim)
		raise IndexError(msg)
	if axis < 0:
		axis += result_ndim
	
	sl = (slice(None),) * axis + (_nx.newaxis,)
	expanded_arrays = [arr[sl] for arr in arrays]
	return _nx.concatenate(expanded_arrays, axis=axis)'''

#################################################################
# Plot spectra
def plotCorrelations(filename):
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 3.0

	nCols = 2
	chosenCols = [0, 1, 3, 5, 7]

	# read in file
	data = loadtxt(filename, usecols=tuple(chosenCols))

	#print 'Loading data file...'
	ax.plot(data[:,0], data[:,1], linestyle='-', color='red', linewidth=lw, label='Lattice')
	ax.plot(data[:,0], data[:,2], linestyle='-', color='blue', linewidth=lw, label=r'$2\pi D_Q T = 0.5$')
	ax.plot(data[:,0], data[:,3], linestyle='-', color='green', linewidth=lw, label=r'$2\pi D_Q T = 1.0$')
	ax.plot(data[:,0], data[:,4], linestyle='-', color='purple', linewidth=lw, label=r'$2\pi D_Q T = 1.5$')
	ax.axhline(0.0, color='black', linewidth=1)
	
	ax.set_xlabel(r'$\Delta y$', fontsize = labelsize + 10)
	ax.set_ylabel(r'$\left< \delta \left( \frac{dN}{dy_1} \right) \delta \left( \frac{dN}{dy_2} \right) \right> \left< \frac{dN}{dy} \right>^{-1}$', fontsize = labelsize + 10)
	ax.legend(loc=0, ncol=1, prop={'size': plotfontsize+5})
	#text(0.9, 0.15, r'(b)', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes, size=30)
	#plt.title(pathname)
	
	#plt.show(block=False)
	#plt.show()
	outfilename = os.path.splitext(filename)[0] + '.pdf'
	plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	print 'Saved to', outfilename





def generate_all_plots():
	plotCorrelations('ecf_pion.out')
	plotCorrelations('ecf_proton.out')
	#pause()


if __name__ == "__main__":
	generate_all_plots()
	print 'Finished all.'

# End of file
