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

#################################################################
# Parameters to run script with
#################################################################
#file information
#filename = sys.argv[1]

#grid information
nDY = 241
nTrajectories = 1
nTauD = 3

qAxisColors = ['red', 'blue', 'green']
cmpStyles = ['-', '--']

panelLabels = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)']
panelCounter = 0

chosenParticle = 'pion'
chosenTauDs = ['0.0','0.5','1.0']

hbarC = 0.197327053

#################################################################
# Shouldn't have to change anything below this line
#################################################################

def pause():
    programPause = raw_input("Press the <ENTER> key to continue...")

def stack(arrays, axis=0):
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
	return _nx.concatenate(expanded_arrays, axis=axis)

#################################################################
# Plot spectra
def plotSpectra(filenamestem, tauDvals):
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 1.5

	filename1 = filenamestem + '_tauD_' + tauDvals[0] + '.out'
	filename2 = filenamestem + '_tauD_' + tauDvals[1] + '.out'
	filename3 = filenamestem + '_tauD_' + tauDvals[2] + '.out'
	nCols = 2
	chosenCols = [0, 1]
	dims = [nDY, nCols]

	# read in file
	data1 = loadtxt(filename1, usecols=tuple(chosenCols)).reshape(dims)
	data2 = loadtxt(filename2, usecols=tuple(chosenCols)).reshape(dims)
	data3 = loadtxt(filename3, usecols=tuple(chosenCols)).reshape(dims)
	data = stack([data1, data2, data3])

	#print 'Loading data file...'
	ax.plot(data[0,:,0], data[0,:,1], linestyle='-', color='red', linewidth=lw, marker='s', markevery=4, label=r'$\tau_D = %(td)s$ fm' % {'td': tauDvals[0]})
	ax.plot(data[1,:,0], data[1,:,1], linestyle='-', color='green', linewidth=lw,  marker='*', markevery=4, label=r'$\tau_D = %(td)s$ fm' % {'td': tauDvals[1]})
	ax.plot(data[2,:,0], data[2,:,1], linestyle='-', color='blue', linewidth=lw,  marker='o', markevery=4, label=r'$\tau_D = %(td)s$ fm' % {'td': tauDvals[2]})

	ax.axhline(0.0, color='black', linewidth=1)
	lowerLims1 = array([0.0, 2.0 * amin(data[0:nTauD,:,1])])
	lowerLims2 = array([0.0, 2.0 * amin(data[0:nTauD,:,1])])
	lowerLims3 = array([0.0, 2.0 * amin(data[0:nTauD,:,1])])
	upperLims1 = array([amax(data[0:nTauD,:,0]), 1.1 * amax(data[0:nTauD,:,1])])
	upperLims2 = array([amax(data[0:nTauD,:,0]), 1.1 * amax(data[0:nTauD,:,1])])
	upperLims3 = array([amax(data[0:nTauD,:,0]), 1.1 * amax(data[0:nTauD,:,1])])
	lowerLims = minimum( minimum(lowerLims1, lowerLims2), lowerLims3)
	upperLims = maximum( maximum(upperLims1, upperLims2), upperLims3)
	ax.axis(hstack( zip(lowerLims,upperLims) ))
	
	ax.set_xlabel(r'$\Delta y$', {'fontsize': plotfontsize + 5})
	ax.set_ylabel(r'$C(\Delta y)$', {'fontsize': plotfontsize + 5})
	ax.legend(loc=0, prop={'size': plotfontsize+5})
	#plt.title(filenamestem)
	
	#plt.show(block=False)
	#plt.show()
	outfilename = filenamestem + '_vs_tauD_wo_muB.pdf'
	plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	print 'Saved to', outfilename


#################################################################
# Plot HBT correlations
def plotHBT(filenamestem, tauDvals):
	# set-up
	plotfontsize = 12
	fig, ax = plt.subplots(1, 1)
	fig.subplots_adjust(wspace=0.0, hspace=0.0) 
	lw = 1.5

	filename1 = filenamestem + '_tauD_' + tauDvals[0] + '.out'
	filename2 = filenamestem + '_tauD_' + tauDvals[1] + '.out'
	filename3 = filenamestem + '_tauD_' + tauDvals[2] + '.out'
	nCols = 2
	chosenCols = [0, 6]
	dims = [nDY, nCols]

	# read in file
	data1 = loadtxt(filename1, usecols=tuple(chosenCols)).reshape(dims)
	data2 = loadtxt(filename2, usecols=tuple(chosenCols)).reshape(dims)
	data3 = loadtxt(filename3, usecols=tuple(chosenCols)).reshape(dims)
	data = stack([data1, data2, data3])

	#ax.plot(data[0,0,:,2], zeros(nDY)-1000.0, linestyle='-', color='black', linewidth=lw, label=r'$\tau_D = %(td)s$ fm' % {'td': tauDvals[0]})
	#ax.plot(data[0,1,:,2], zeros(nDY)-1000.0, linestyle='--', color='black', linewidth=lw, label=r'$\tau_D = %(td)s$ fm' % {'td': tauDvals[1]})
	#ax.plot(data[0,2,:,2], zeros(nDY)-1000.0, linestyle=':', color='black', linewidth=lw, label=r'$\tau_D = %(td)s$ fm' % {'td': tauDvals[2]})

	ax.plot(data[0,:,0], data[0,:,1], linestyle='-', color='red', linewidth=lw, marker='s', markevery=4, label=r'$\tau_D = %(td)s$ fm' % {'td': tauDvals[0]})
	ax.plot(data[1,:,0], data[1,:,1], linestyle='-', color='green', linewidth=lw, marker='*', markevery=4, label=r'$\tau_D = %(td)s$ fm' % {'td': tauDvals[1]})
	ax.plot(data[2,:,0], data[2,:,1], linestyle='-', color='blue', linewidth=lw, marker='o', markevery=4, label=r'$\tau_D = %(td)s$ fm' % {'td': tauDvals[2]})

	ax.axhline(0.0, color='black', linewidth=1)
	#lowerLims1 = array([0.0, 1.1 * amin(data[0,:,3])])
	#lowerLims2 = array([0.0, 1.1 * amin(data[1,:,3])])
	#lowerLims3 = array([0.0, 1.1 * amin(data[2,:,3])])
	#upperLims1 = array([amax(data[0,:,0]), 1.1 * amax(data[0,:,3])])
	#upperLims2 = array([amax(data[1,:,0]), 1.1 * amax(data[1,:,3])])
	#upperLims3 = array([amax(data[2,:,0]), 1.1 * amax(data[2,:,3])])
	#print lowerLims1, lowerLims2, lowerLims3
	#print upperLims1, upperLims2, upperLims3
	#lowerLims = minimum( minimum(lowerLims1, lowerLims2), lowerLims3)
	#upperLims = maximum( maximum(upperLims1, upperLims2), upperLims3)
	lowerLims = array([0.0, 2.0 * amin(data[0:nTauD,:,1])])
	upperLims = array([amax(data[0:nTauD,:,0]), 1.1 * amax(data[0:nTauD,:,1])])
	#print data.shape, lowerLims, upperLims
	ax.axis(hstack( zip(lowerLims,upperLims) ))
	
	ax.set_xlabel(r'$\Delta y$', {'fontsize': plotfontsize + 5})
	ax.set_ylabel(r'$C_{HBT}(\Delta y)$', {'fontsize': plotfontsize + 5})
	ax.legend(loc=0, ncol=1, prop={'size': plotfontsize+5})
	#plt.title(filename)
	
	#plt.show(block=False)
	#plt.show()
	outfilename = filenamestem + '_vs_tauD_wo_muB.pdf'
	plt.savefig(outfilename, format='pdf', bbox_inches='tight')
	print 'Saved to', outfilename


def generate_all_plots():
	plotSpectra('results_1sp/%(cp)s/%(cp)s_results_1sp_pointA' % {'cp': chosenParticle}, chosenTauDs)
	#plotSpectra('results_1sp/%(cp)s/%(cp)s_results_1sp_pointB' % {'cp': chosenParticle}, chosenTauDs)
	plotSpectra('results_1sp/%(cp)s/%(cp)s_results_1sp_midpoint' % {'cp': chosenParticle}, chosenTauDs)
	#plotSpectra('results_1sp/%(cp)s/%(cp)s_results_1sp_tauB_transportA' % {'cp': chosenParticle}, chosenTauDs)
	plotSpectra('results_1sp/%(cp)s/%(cp)s_results_1sp_tauA_transportB' % {'cp': chosenParticle}, chosenTauDs)
	plotHBT('results_HBT/%(cp)s/%(cp)s_results_HBT_pointA' % {'cp': chosenParticle}, chosenTauDs)
	#plotHBT('results_HBT/%(cp)s/%(cp)s_results_HBT_pointB' % {'cp': chosenParticle}, chosenTauDs)
	plotHBT('results_HBT/%(cp)s/%(cp)s_results_HBT_midpoint' % {'cp': chosenParticle}, chosenTauDs)
	#plotHBT('results_HBT/%(cp)s/%(cp)s_results_HBT_tauB_sigmaTA' % {'cp': chosenParticle}, chosenTauDs)
	plotHBT('results_HBT/%(cp)s/%(cp)s_results_HBT_tauA_transportB' % {'cp': chosenParticle}, chosenTauDs)
	pause()


if __name__ == "__main__":
	generate_all_plots()
	print 'Finished all.'

# End of file
