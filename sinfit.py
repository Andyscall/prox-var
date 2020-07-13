import math
import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.optimize import curve_fit


sigmas = [1.0, 1.5, 2.0, 2.5, 3.0]
bins = 1
tail = [""]#["_old", "_new"]

for sigma in sigmas:
	file_1 = 'scaled_results_proxima_{}_{}_new.csv'.format(sigma, bins)
	#print(file_1)
	mjd_1 = []
	
	a_1 = []
	ae_1 = []
	
	title = []

	with open(file_1, 'rt') as csvfile:
		data = csv.reader(csvfile, delimiter=' ')
	
		for row in data:

			if len(row) < 16:
				title = row
				continue
			
			if float(row[19]) < 4:
				continue

			if float(row[0]) < 58440:
				continue
			
			else:
				mjd_1.append(float(row[0]))

				a_1.append(float(row[19]))
				ae_1.append(float(row[20]))

	x_0 = np.mean(mjd_1)
	
	def func(x, a, b, c, p):
		return(a * np.sin(b * (x - x_0 + p)) + c)

	popt, pcov = curve_fit(
		func, mjd_1, a_1, 
		p0=[0.01, 0.0728, 11.8, 20], 
		#sigma=ae_1, 
		bounds=([0.0, 0.0, 10.0, 0.0], 
			[1.0, 1.0, 20.0, 100.0]))#, 
		#method='trf')

	print(sigma, "params", popt)

	#popt = [0.04, 0.0628, 11.8, 25]

	curve = []
	for i in range(len(mjd_1)):
		w = func(mjd_1[i], *popt)
		curve.append(w)

	plt.plot(mjd_1, curve, 'r-', label='fit: a=%5.3f, b=%5.3f, c=%5.3f, p=%5.3f' % tuple(popt))
	plt.plot(mjd_1, a_1, 'b.', label='data')
	plt.xlabel("MJD")
	plt.ylabel("Mag")
	plt.legend()
	plt.title('fit: {} * sin({} * (x-x0 + {})) + {}'.format(popt[0], popt[1], popt[3], popt[2]), fontsize=8)
	#plt.show()
	plt.savefig("cf_plot_p19_{}_{}.pdf".format(sigma, bins))
	plt.close()