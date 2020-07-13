import math
import numpy as np
import matplotlib.pyplot as plt
import csv


'''
gG = np.genfromtxt(open("proxima/lcgoodg_bE.txt"), names=True, delimiter=" ")
gH = np.genfromtxt(open("proxima/lcgoodg_bF.txt"), names=True, delimiter=" ")
gl = np.genfromtxt(open("proxima/lcgoodg_bi.txt"), names=True, delimiter=" ")
go = np.genfromtxt(open("proxima/lcgoodg_bj.txt"), names=True, delimiter=" ")
gp = np.genfromtxt(open("proxima/lcgoodg_bm.txt"), names=True, delimiter=" ")
gn = np.genfromtxt(open("proxima/lcgoodg_bn.txt"), names=True, delimiter=" ")
'''

gE = np.genfromtxt(open("proxima/20f_lcprox_bE.txt"), names=True, delimiter=" ")
gF = np.genfromtxt(open("proxima/20f_lcprox_bF.txt"), names=True, delimiter=" ")
gi = np.genfromtxt(open("proxima/20f_lcprox_bi.txt"), names=True, delimiter=" ")
gj = np.genfromtxt(open("proxima/20f_lcprox_bj.txt"), names=True, delimiter=" ")
gm = np.genfromtxt(open("proxima/20f_lcprox_bm.txt"), names=True, delimiter=" ")
gn = np.genfromtxt(open("proxima/20f_lcprox_bn.txt"), names=True, delimiter=" ")



def bin_func(lcgVbe, sigma=2.0, corr=0.0, bin_width=1):
	"""
	Takes in a file (needs to be imported separately) where the 0th column is MJD, 
	5th column is magnitude, and 6th column is error. Bins the data by time (roughly)
	then takes the data points that lie within 2 sigma.
	Makes a scatter plot.
	Returns list containing lists of MJD, magnitudes, and errors.

	all these things have "Ve" in the names because my test case was lcgoodV_be.dat
	"""
	mjd_Ve = []
	mag_Ve = []
	mag_err_Ve = []

	for i in range(len(lcgVbe)):
		mjd_Ve.append(lcgVbe[i][0])
		mag_Ve.append(lcgVbe[i][5])
		mag_err_Ve.append(lcgVbe[i][6])

	vals_vE = [mjd_Ve, mag_Ve, mag_err_Ve]

	return(vals_vE)


Efile = bin_func(gE)
Ffile = bin_func(gF)
ifile = bin_func(gi)
jfile = bin_func(gj)
mfile = bin_func(gm)
nfile = bin_func(gn)


plt.scatter(Efile[0], Efile[1], color='r', marker='+')
plt.scatter(Ffile[0], Ffile[1], color='g', marker='+')
plt.scatter(ifile[0], ifile[1], color='b', marker='+')
plt.scatter(jfile[0], jfile[1], color='c', marker='+')
plt.scatter(mfile[0], mfile[1], color='m', marker='+')
plt.scatter(nfile[0], nfile[1], color='k', marker='+')
plt.show()







