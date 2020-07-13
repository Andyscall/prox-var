import math
import numpy as np
import matplotlib.pyplot as plt
import csv

sigmas = [1.0, 1.5, 2.0, 2.5, 3.0]
bins = 1
#star = ["HD82443", "HD20630", "gj358", "proxima", "proxima"]
#ylims = [(7.6, 8.2), (5.6, 6.0), (11.35, 11.55), (11.5, 12.0), (11.5, 12.0)]
#tail = ["", "", "", "old", "_new"]

star = ["gj358", "proxima", "proxima"]
ylims = [(11.35, 11.55), (11.6, 11.9), (11.6, 11.9)]
tail = ["", "_old", "_new"]

for i in range(len(star)):
	for sigma in sigmas:
		file_1 = 'scaled_results_{}_{}_{}{}.csv'.format(star[i], sigma, bins, tail[i])
		print(file_1)
		#print(file_1)
		mjd_1 = []
		c1_1 = []
		s1_1 = []
		e1_1 = []
		c2_1 = []
		s2_1 = []
		e2_1 = []
		c3_1 = []
		s3_1 = []
		e3_1 = []
		c4_1 = []
		s4_1 = []
		e4_1 = []
		c5_1 = []
		s5_1 = []
		e5_1 = []
		c6_1 = []
		s6_1 = []
		e6_1 = []
		a_1 = []
		ae_1 = []
		title = []

		with open(file_1, 'rt') as csvfile:
			data = csv.reader(csvfile, delimiter=' ')
			
			if star == 'proxima':
				for row in data:
					if len(row) < 16:
						title = row
						continue
					mjd_1.append(float(row[0]))

					c1_1.append(float(row[1]))
					e1_1.append(float(row[2]))
					s1_1.append(float(row[3]))

					c2_1.append(float(row[4]))
					e2_1.append(float(row[5]))
					s2_1.append(float(row[6]))

					c3_1.append(float(row[7]))
					e3_1.append(float(row[8]))
					s3_1.append(float(row[9]))

					c4_1.append(float(row[10]))
					e4_1.append(float(row[11]))
					s4_1.append(float(row[12]))

					c5_1.append(float(row[13]))
					e5_1.append(float(row[14]))
					s5_1.append(float(row[15]))

					c6_1.append(float(row[16]))
					e6_1.append(float(row[17]))
					s6_1.append(float(row[18]))

					a_1.append(float(row[19]))
					ae_1.append(float(row[20]))

				

			else:
				for row in data:
					if len(row) < 16:
						title = row
						continue
					mjd_1.append(float(row[0]))

					c1_1.append(float(row[1]))
					e1_1.append(float(row[2]))
					s1_1.append(float(row[3]))

					c2_1.append(float(row[4]))
					e2_1.append(float(row[5]))
					s2_1.append(float(row[6]))

					c3_1.append(float(row[7]))
					e3_1.append(float(row[8]))
					s3_1.append(float(row[9]))

					c4_1.append(float(row[10]))
					e4_1.append(float(row[11]))
					s4_1.append(float(row[12]))

					c5_1.append(float(row[13]))
					e5_1.append(float(row[14]))
					s5_1.append(float(row[15]))

					a_1.append(float(row[16]))
					ae_1.append(float(row[17]))
		
		if star != 'proxima':
			print(np.shape(e5_1))
			"""
			plt.scatter(mjd_1, c1_1, color='xkcd:lavender', marker='+', label='gF')
			plt.scatter(mjd_1, s1_1, color='xkcd:purpleish', marker='x', label='gF cut')
			plt.scatter(mjd_1, c2_1, color='xkcd:seafoam green', marker='+', label='gE')
			plt.scatter(mjd_1, s2_1, color='xkcd:blue green', marker='x', label='gE cut')
			plt.scatter(mjd_1, c3_1, color='xkcd:azure', marker='+', label='gi')
			plt.scatter(mjd_1, s3_1, color='xkcd:royal blue', marker='x', label='gi cut')
			plt.scatter(mjd_1, c4_1, color='xkcd:light orange', marker='+', label='gm')
			plt.scatter(mjd_1, s4_1, color='xkcd:dark orange', marker='x', label='gm cut')
			plt.scatter(mjd_1, c5_1, color='xkcd:cherry red', marker='+', label='gj')
			plt.scatter(mjd_1, s5_1, color='xkcd:cranberry', marker='x', label='gj cut')
			plt.scatter(mjd_1, c6_1, color='xkcd:wheat', marker='+', label='gn')
			plt.scatter(mjd_1, s6_1, color='xkcd:gold', marker='x', label='gn cut')
			xlim = (58200.0, 58900.0)
			plt.xlim(xlim[:,0])
			plt.ylim(ylims[i])
			plt.xlabel("MJD")
			plt.ylabel("Mag")
			plt.legend()
			plt.title(title, fontsize=8)
			#plt.show()
			plt.savefig("scaled_plot_{}_{}_{}{}.pdf".format(star[i], sigma, bins, tail[i]))
			plt.close()
			"""
			'''
			elif star != "gj358":
			#plt.scatter(mjd_1, a_1, color='k', marker='+', label='average')
			plt.scatter(mjd_1, c1_1, color='xkcd:lavender', marker='+', label='gA')
			plt.scatter(mjd_1, s1_1, color='xkcd:purpleish', marker='x', label='gA cut')
			plt.scatter(mjd_1, c2_1, color='xkcd:seafoam green', marker='+', label='gE')
			plt.scatter(mjd_1, s2_1, color='xkcd:blue green', marker='x', label='gE cut')
			plt.scatter(mjd_1, c3_1, color='xkcd:azure', marker='+', label='gi')
			plt.scatter(mjd_1, s3_1, color='xkcd:royal blue', marker='x', label='gi cut')
			plt.scatter(mjd_1, c4_1, color='xkcd:light orange', marker='+', label='gm')
			plt.scatter(mjd_1, s4_1, color='xkcd:dark orange', marker='x', label='gm cut')
			plt.scatter(mjd_1, c5_1, color='xkcd:cherry red', marker='+', label='gq')
			plt.scatter(mjd_1, s5_1, color='xkcd:cranberry', marker='x', label='gq cut')
			plt.ylim(ylims[i])
			plt.xlabel("MJD")
			plt.ylabel("Mag")
			plt.legend()
			plt.title(title, fontsize=8)
			#plt.show()
			plt.savefig("scaled_plot_{}_{}_{}{}.pdf".format(star[i], sigma, bins, tail[i]))
			plt.close()
			'''
		else:
			#plt.scatter(mjd_1, a_1, color='k', marker='+', label='average')
			plt.scatter(mjd_1, c1_1, color='xkcd:lavender', marker='+', label='gH')
			plt.scatter(mjd_1, s1_1, color='xkcd:purpleish', marker='x', label='gH cut')
			plt.scatter(mjd_1, c2_1, color='xkcd:seafoam green', marker='+', label='gG')
			plt.scatter(mjd_1, s2_1, color='xkcd:blue green', marker='x', label='gG cut')
			plt.scatter(mjd_1, c3_1, color='xkcd:azure', marker='+', label='gl')
			plt.scatter(mjd_1, s3_1, color='xkcd:royal blue', marker='x', label='gl cut')
			plt.scatter(mjd_1, c4_1, color='xkcd:light orange', marker='+', label='go')
			plt.scatter(mjd_1, s4_1, color='xkcd:dark orange', marker='x', label='go cut')
			plt.scatter(mjd_1, c5_1, color='xkcd:cherry red', marker='+', label='gp')
			plt.scatter(mjd_1, s5_1, color='xkcd:cranberry', marker='x', label='gp cut')
			plt.ylim(ylims[i])
			plt.xlabel("MJD")
			plt.ylabel("Mag")
			plt.legend()
			plt.title(title, fontsize=8)
			#plt.show()
			plt.savefig("scaled_plot_{}_{}_{}.pdf".format(star[i], sigma, bins))
			plt.close()

"""
sigmas = [1.0, 1.5, 2.0, 2.5, 3.0]
bins = 1

for sigma in sigmas:
	file_1 = 'scaled_results_HD20630_{}_{}.csv'.format(sigma, bins)
	#print(file_1)
	mjd_1 = []
	c1_1 = []
	s1_1 = []
	e1_1 = []
	c2_1 = []
	s2_1 = []
	e2_1 = []
	c3_1 = []
	s3_1 = []
	e3_1 = []
	c4_1 = []
	s4_1 = []
	e4_1 = []
	c5_1 = []
	s5_1 = []
	e5_1 = []
	c6_1 = []
	s6_1 = []
	e6_1 = []
	c7_1 = []
	s7_1 = []
	e7_1 = []
	a_1 = []
	ae_1 = []

	with open(file_1, 'rt') as csvfile:
		star = csv.reader(csvfile, delimiter=' ')
		
		for row in star:

			if len(row) < 16:
				title = row
				continue
			mjd_1.append(float(row[0]))
			c1_1.append(float(row[1]))
			e1_1.append(float(row[2]))
			s1_1.append(float(row[3]))
			c2_1.append(float(row[4]))
			e2_1.append(float(row[5]))
			s2_1.append(float(row[6]))
			c3_1.append(float(row[7]))
			e3_1.append(float(row[8]))
			s3_1.append(float(row[9]))
			c4_1.append(float(row[10]))
			e4_1.append(float(row[11]))
			s4_1.append(float(row[12]))
			c5_1.append(float(row[13]))
			e5_1.append(float(row[14]))
			s5_1.append(float(row[15]))
			c6_1.append(float(row[16]))
			e6_1.append(float(row[17]))
			s6_1.append(float(row[18]))
			c7_1.append(float(row[19]))
			e7_1.append(float(row[20]))
			s7_1.append(float(row[21]))
			a_1.append(float(row[22]))
			ae_1.append(float(row[23]))


	plt.scatter(mjd_1, c1_1, color='xkcd:lavender', marker='+', label='gA')
	plt.scatter(mjd_1, s1_1, color='xkcd:purpleish', marker='x', label='gA cut')
	plt.scatter(mjd_1, c2_1, color='xkcd:seafoam green', marker='+', label='gE')
	plt.scatter(mjd_1, s2_1, color='xkcd:blue green', marker='x', label='gE cut')
	plt.scatter(mjd_1, c3_1, color='xkcd:azure', marker='+', label='gi')
	plt.scatter(mjd_1, s3_1, color='xkcd:royal blue', marker='x', label='gi cut')
	plt.scatter(mjd_1, c4_1, color='xkcd:light orange', marker='+', label='gm')
	plt.scatter(mjd_1, s4_1, color='xkcd:dark orange', marker='x', label='gm cut')
	plt.scatter(mjd_1, c5_1, color='xkcd:cherry red', marker='+', label='gq')
	plt.scatter(mjd_1, s5_1, color='xkcd:cranberry', marker='x', label='gq cut')
	plt.scatter(mjd_1, c6_1, color='xkcd:baby pink', marker='+', label='Va')
	plt.scatter(mjd_1, s6_1, color='xkcd:dull pink', marker='x', label='Va cut')
	plt.scatter(mjd_1, c7_1, color='xkcd:maize', marker='+', label='Ve')
	plt.scatter(mjd_1, s7_1, color='xkcd:goldenrod', marker='x', label='Ve cut')
	#plt.scatter(mjd_1, a_1, color='k', marker='+', label='average')
	plt.ylim(5.6, 6.0)
	plt.xlabel("MJD")
	plt.ylabel("Mag")
	plt.legend()
	plt.title(title, fontsize=8)
	plt.savefig("scaled_plot_HD20630_{}_{}.pdf".format(sigma, bins))
	plt.close()
"""