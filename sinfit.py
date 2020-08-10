import math
import numpy as np
import matplotlib.pyplot as plt
import csv
from scipy.optimize import curve_fit


#sigmas = [1.0, 1.5, 2.0, 2.5, 3.0]
sigmas = [2.0]
bins = 1
tail = [""]#["_old", "_new"]

for sigma in sigmas:
	file_1 = 'scaled_results_proxima_{}_{}.csv'.format(sigma, bins)
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
	

	#chi square loops and curve fitting
	def func(x, a, b, c, p):
		return(a * np.sin((2*math.pi/b) * (x - x_0 + p)) + c)

	params = [0.01, 86.3, 11.8, 20]

	#old b: 0.0728

	X = 99999

	for a in range(-5, 5):
		for b in range(-5, 5):
			for c in range(-5, 5):
				for p in range(0, 10):

					param = [0.01, 86.3, 11.8, 20]

					pa = param[0] + a*0.001
					pb = param[1] + b#*0.1
					pc = param[2] + c*0.1
					pp = param[3] + p

					popt, pcov = curve_fit(
						func, mjd_1, a_1, 
						p0=[pa, pb, pc, pp], 
						bounds=([0.0, 0.0, 5.0, 0.0], 
							[1.0, 200.0, 25.0, 100.0]))

					#print(sigma, "param", popt)

					cs = []
					for i in range(len(mjd_1)):
						w = func(mjd_1[i], *popt)
						cs.append(w)

					chis = []

					for i in range(len(mjd_1)):
						x = ((a_1[i] - cs[i])**2)/cs[i]
						chis.append(x)

					chisq = sum(chis)
					if chisq < X:
						#print("new")
						#print(a, b, c, p)
						X = chisq
						params[0] = pa
						params[1] = pb
						params[2] = pc
						params[3] = pp
					#print("p is {}".format(p))
				#print("c is {}".format(c))
			#print("b is {}".format(b))
		print("a is {}".format(a))
		#print(sigma, "param", popt)
		#print("chisq is {}".format(chisq))


	#curve fitting
	popt, pcov = curve_fit(
		func, mjd_1, a_1, 
		p0=params, 
		#sigma=ae_1, 
		bounds=([0.0, 0.0, 10.0, 0.0], 
			[1.0, 200.0, 20.0, 100.0]))#, 
		#method='trf')

	print(sigma, "params", popt)


	#popt = [0.04, 0.0628, 11.8, 25]

	curve = []
	ressq = []
	for i in range(len(mjd_1)):
		w = func(mjd_1[i], *popt)
		curve.append(w)
		r = (a_1[i] - w)**2
		ressq.append(r)


	#for sigma_Y we use the root mean square of the residuals
	sigma_Y = np.sqrt(np.mean(ressq))

	N = len(mjd_1)
	T = mjd_1[N-1] - mjd_1[0]

	sigma_a = (2/N)**(0.5) * sigma_Y

	sigma_f = (6/N)**(0.5) * sigma_Y/(math.pi * T * params[0])

	sigma_P = sigma_f * (params[1]**2)

	print("sigma_Y is {}, sigma_A is {}, sigma_P is {}".format(sigma_Y, sigma_a, sigma_P))
	


	plt.plot(mjd_1, curve, 'r-', label='fit: A=%5.3f, P=%5.3f, c=%5.3f, phi=%5.3f' % tuple(popt))
	plt.plot(mjd_1, a_1, 'b.', label='data')
	plt.xlabel("MJD")
	plt.ylabel("Mag")
	plt.legend()
	plt.title("sigma_Y = {}, sigma_a = {}, sigma_P = {}".format(sigma_Y, sigma_a, sigma_P), fontsize=8)
	plt.suptitle('{}*sin((2pi/{})*(x-x0+{}))+{}'.format(popt[0], popt[1], popt[3], popt[2]), fontsize=8)
	#plt.show()
	plt.savefig("cf_plot_p19_{}_{}.pdf".format(sigma, bins))
	plt.close()