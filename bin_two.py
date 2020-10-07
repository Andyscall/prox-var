import math
import numpy as np
import matplotlib.pyplot as plt
import csv

np.seterr(all='raise')



gE = np.genfromtxt(open("HD20630/lcgoodg_bE.dat"), names=True, delimiter=" ")
gi = np.genfromtxt(open("HD20630/lcgoodg_bi.dat"), names=True, delimiter=" ")
gA = np.genfromtxt(open("HD20630/lcgoodg_bA.dat"), names=True, delimiter=" ")
gm = np.genfromtxt(open("HD20630/lcgoodg_bm.dat"), names=True, delimiter=" ")
gq = np.genfromtxt(open("HD20630/lcgoodg_bq.dat"), names=True, delimiter=" ")


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

	b_mjd_Ve = []
	b_mag_Ve = []
	b_mag_err_Ve = []
	b_mag_err_Ve_w = []
	b_mjd_means = []

	for i in np.arange(58000, 59200, bin_width):
		
		b_mjds = []
		b_i = [] #temporary list to hold magnitudes
		b_e_i_p = [] #temporary list to hold errors in preparation for summing up, which is why the math is weird (error propagation)
		b_w = [] #temporary list to hold weights for magnitudes (1/(err)**2)

		if bin_width <= 1:
			bin_nums = 1/bin_width
			steps = bin_width
		else:
			bin_nums = 1
			steps = 1

		for b in np.arange(0, bin_nums, steps):

			for j in range(len(mjd_Ve)):

				if (int(mjd_Ve[j]) + b)==i:
					b_i.append(mag_Ve[j] + corr)
					b_e_i_p.append((mag_err_Ve[j]/(mag_Ve[j] + corr))**2)
					b_w.append(1/(mag_err_Ve[j])**2)
					b_mjds.append(mjd_Ve[j])
		
		if len(b_i)!=0:
			mjd_mean = np.average(b_mjds, weights=b_w)
			b_mjd_means.append(mjd_mean)
			b_mjd_Ve.append(i)
			w_mean = np.average(b_i, weights=b_w)
			b_mag_Ve.append(w_mean)
			b_mag_err_Ve.append(w_mean*math.sqrt((sum(b_e_i_p))))
			b_mag_err_Ve_w.append(1/((w_mean*(math.sqrt(sum(b_e_i_p))))**2))
		
		elif len(b_i) == 0:
			# adds a error-filled zero mag, so that the lists have the same number of data points
			b_mjd_Ve.append(i)
			b_mjd_means.append(0)
			b_mag_Ve.append(0)
			b_mag_err_Ve.append(9999)
			b_mag_err_Ve_w.append(0)

	mmag_Ve = np.average(b_mag_Ve, weights=b_mag_err_Ve_w)
	v_list = []

	for i in range(len(b_mag_Ve)):
		#does sums for the weighted standard deviation
		vV = b_mag_Ve[i] - mmag_Ve
		v_list.append(b_mag_err_Ve_w[i] * (vV**2))

	stdm_Ve = math.sqrt(sum(v_list)/sum(b_mag_err_Ve_w))
	mjds_Ve = []
	mags_Ve = []
	mags_err_Ve = []

	clipped_mag = []

	for i in range(len(b_mag_Ve)):
		#excludes values outside of however many sigma (default 2)
		if (mmag_Ve-(sigma*stdm_Ve)) < b_mag_Ve[i] and b_mag_Ve[i] < (mmag_Ve+(sigma*stdm_Ve)):
			mjds_Ve.append(b_mjd_Ve[i])
			mags_Ve.append(b_mag_Ve[i])
			mags_err_Ve.append(b_mag_err_Ve[i])
			clipped_mag.append(0)

		else:
			#adds an error-filled zero, so as to not have empty data points. Not sure this works?
			mjds_Ve.append(b_mjd_Ve[i])
			mags_Ve.append(0)
			mags_err_Ve.append(99999)
			clipped_mag.append(b_mag_Ve[i])

	vals_vE = [mjds_Ve, mags_Ve, mags_err_Ve, clipped_mag]

	return(vals_vE)


def mjd_func(lcgVbe, bin_width=1):
	mjd_Ve = []
	mag_err_Ve = []

	for i in range(len(lcgVbe)):
		mjd_Ve.append(lcgVbe[i][0])
		mag_err_Ve.append(lcgVbe[i][6])

	b_mjd_Ve = []
	b_mjd_means = []

	for i in np.arange(58000, 59200, bin_width):
		
		b_mjds = []
		b_w = [] #temporary list to hold weights for magnitudes (1/(err)**2)
		b_i = [] #temporary list to hold magnitudes

		if bin_width <= 1:
			bin_nums = 1/bin_width
			steps = bin_width
		else:
			bin_nums = 1
			steps = 1

		for b in np.arange(0, bin_nums, steps):

			for j in range(len(mjd_Ve)):

				if (int(mjd_Ve[j]) + b)==i:
					b_i.append(3)
					b_w.append(1/(mag_err_Ve[j]**2))
					#print("err", 1/(mag_err_Ve[j])**2)
					b_mjds.append(mjd_Ve[j])
		
		if len(b_i)!=0:
			mjd_mean = np.average(b_mjds, weights=b_w)
			b_mjd_means.append(mjd_mean)
			
		elif len(b_i) == 0:
			b_mjd_means.append(i)
	
	return(b_mjd_means)

def chi_func(mrange = range(-5, 5), qrange = range(-5, 5), Arange = range(-5, 5), Erange = range(-5, 5),
	Efac = 0.015, Afac = -0.124, qfac = -0.026, mfac = 0.093):
	bestE = 0.0
	besti = 0.0
	bestm = 0.0
	bestq = 0.0
	bestA = 0.0
	besta = 0.0
	beste = 0.0
	bestX = 9.99*10**10


	for A in Arange:
		for m in mrange:
			for q in qrange:
				for E in Erange:
					
					chisum = 0.0

					#technically these lists for the I data are redundant, but I put them in for consistency
					#The "scale_XX_w" list is for weights; it's equal to the inverse square err list.
					scalei = 0
					scale_gi = []
					scale_gi_err = []
					scale_gi_w = []

					scaleE = E*0.001 + Efac
					scale_gE = []
					scale_gE_err = []
					scale_gE_w = []

					scaleq = q*0.001 + qfac
					scale_gq = []
					scale_gq_err = []
					scale_gq_w = []

					scalem = m*0.001 + mfac
					scale_gm = []
					scale_gm_err = []
					scale_gm_w = []

					scaleA = A*0.001 + Afac
					scale_gA = []
					scale_gA_err = []
					scale_gA_w = []

					for j in range(len(binned_gA[1])):
						#This loop will multiply each of the data by their scale factor.
						#Not sure if the error needs any more manipulation done to it?
						if np.isnan(binned_gi[1][j]) or binned_gi[1][j] == 0.0:
							#this is necessary so that the program doesn't try to match everything to a bunch of zeroes
							continue
						else:
							sgi = binned_gi[1][j]
							sgie = binned_gi[2][j]
							sgiw = 1/(sgie**2)
						scale_gi.append(sgi)
						scale_gi_err.append(sgie)
						scale_gi_w.append(sgiw)

						if np.isnan(binned_gA[1][j]) or binned_gA[1][j] == 0.0:
							sgA = 0
							sgAe = 9999999
							sgAw = 0
						else:
							sgA = binned_gA[1][j]+scaleA
							sgAe = binned_gA[2][j]
							sgAw = 1/(sgAe**2)
						scale_gA.append(sgA)
						scale_gA_err.append(sgAe)
						scale_gA_w.append(sgAw)

						if np.isnan(binned_gE[1][j]) or binned_gE[1][j] == 0.0:
							sgE = 0
							sgEe = 9999999
							sgEw = 0
						else:
							sgE = binned_gE[1][j]+scaleE
							sgEe = binned_gE[2][j]
							sgEw = 1/(sgEe**2)
						scale_gE.append(sgE)
						scale_gE_err.append(sgEe)
						scale_gE_w.append(sgEw)
						
						if np.isnan(binned_gq[1][j]) or binned_gq[1][j] == 0.0:
							sgq = 0
							sgqe = 9999999
							sgqw = 0
						else:
							sgq = binned_gq[1][j]+scaleq
							sgqe = binned_gq[2][j]
							sgqw = 1/(sgqe**2)
						scale_gq.append(sgq)
						scale_gq_err.append(sgqe)
						scale_gq_w.append(sgqw)

						if np.isnan(binned_gm[1][j]) or binned_gm[1][j] == 0.0:
							sgm = 0
							sgme = 9999999
							sgmw = 0
						else:
							sgm = binned_gm[1][j]+scalem
							sgme = binned_gm[2][j]
							sgmw = 1/(sgme**2)
						scale_gm.append(sgm)
						scale_gm_err.append(sgme)
						scale_gm_w.append(sgmw)

					for j in range(len(binned_gi[1])):
						#Takes the weighted average of each bin; each scale_xX[j]*scale_xX_w[j] over the sum of the [j] weights
						try:
							ave = (scale_gE[j]*scale_gE_w[j] + scale_gi[j]*scale_gi_w[j] + scale_gA[j]*scale_gA_w[j] 
								+ scale_gq[j]*scale_gq_w[j] + scale_gm[j]*scale_gm_w[j]
								)/(scale_gE_w[j] + scale_gi_w[j] + scale_gq_w[j] + scale_gm_w[j]
								+ scale_gA_w[j])
						except:
							continue
						
						chi2 = 0
						chiE = ((scale_gE[j] - ave)**2)/(scale_gE_err[j]**2)
						chii = ((scale_gi[j] - ave)**2)/(scale_gi_err[j]**2)
						chiA = ((scale_gA[j] - ave)**2)/(scale_gA_err[j]**2)
						chiq = ((scale_gq[j] - ave)**2)/(scale_gq_err[j]**2)
						chim = ((scale_gm[j] - ave)**2)/(scale_gm_err[j]**2)
						
						chi2 = chiE + chii + chiq + chim + chiA
						if chi2 == 0:
							continue
						else:
							chisum = chisum + chi2
					
					if chisum < bestX:
						#the actual minimization bit; upates the values if necessary
						bestX = chisum
						besti = scalei
						bestA = scaleA
						bestE = scaleE
						bestq = scaleq
						bestm = scalem
						

				#print("q is", q)#, "bestX is", bestX)
			#print("m is", m)
		print("A is", A)

	print("Best chi:", bestX)
	print("scale E:", bestE)
	print("scale q:", bestq)
	print("scale m:", bestm)
	print("scale A:", bestA)
	
	return [bestX, bestE, bestq, bestm, bestA]


def fractioner(binfile):
	points = binfile[1]
	cuts = binfile[3]

	pcount = 0.0
	ccount = 0.0

	for i in range(len(points)):
		if points[i] != 0.0:
			pcount = pcount + 1.0
		if cuts[i] != 0.0:
			ccount = ccount+1.0

	f_data_loss = ccount / (ccount + pcount)

	return f_data_loss


#sigmas = [1.0, 1.5, 2.0, 2.5, 3.0]
sigmas = [2.0]
bins_width = 1

shifts_gA = []
shifts_gE = []
shifts_gi = []
shifts_gm = []
shifts_gq = []
chis = []

fdl_gA = []
fdl_gE = []
fdl_gi = []
fdl_gm = []
fdl_gq = []
mean_fdl = []


for sig in sigmas:

	binned_gE = bin_func(gE, sigma=sig, bin_width=bins_width)
	binned_gA = bin_func(gA, sigma=sig, bin_width=bins_width)
	binned_gi = bin_func(gi, sigma=sig, bin_width=bins_width)
	binned_gm = bin_func(gm, sigma=sig, bin_width=bins_width)
	binned_gq = bin_func(gq, sigma=sig, bin_width=bins_width)
	

	gA_mjd = mjd_func(gA, bin_width=bins_width)
	gE_mjd = mjd_func(gE, bin_width=bins_width)
	gi_mjd = mjd_func(gi, bin_width=bins_width)
	gm_mjd = mjd_func(gm, bin_width=bins_width)
	gq_mjd = mjd_func(gq, bin_width=bins_width)
	
	all_mjd = []


	b_gA = []
	b_gA_err = []
	b_gE = []
	b_gE_err = []
	b_gi = []
	b_gi_err = []
	b_gm = []
	b_gm_err = []
	b_gq = []
	b_gq_err = []
	

	for i in range(len(binned_gA[0])):
		if np.isnan(binned_gE[1][i]) or binned_gE[1][i] == 0.0:
			b_gE.append(0)
			b_gE_err.append(999999)
		else:
			sE = binned_gE[1][i]
			b_gE.append(sE)
			b_gE_err.append(binned_gE[2][i])

		if np.isnan(binned_gi[1][i]) or binned_gi[1][i] == 0.0:
			si = 0
			b_gi.append(0)
			b_gi_err.append(999999)
		else:
			si = binned_gi[1][i]
			b_gi.append(si)
			b_gi_err.append(binned_gi[2][i])

		if np.isnan(binned_gm[1][i]) or binned_gm[1][i] == 0.0:
			sm = 0
			b_gm.append(0)
			b_gm_err.append(999999)
		else:
			sm = binned_gm[1][i]
			b_gm.append(sm)
			b_gm_err.append(binned_gm[2][i])

		if np.isnan(binned_gq[1][i]) or binned_gq[1][i] == 0.0:
			sq = 0
			b_gq.append(0)
			b_gq_err.append(999999)
		else:
			sq = binned_gq[1][i]
			b_gq.append(sq)
			b_gq_err.append(binned_gq[2][i])

		if np.isnan(binned_gA[1][i]) or binned_gA[1][i] == 0.0:
			sA = 0
			b_gA.append(0)
			b_gA_err.append(999999)
		else:
			sA = binned_gA[1][i]
			b_gA.append(sA)
			b_gA_err.append(binned_gA[2][i])

	for i in range(len(gA_mjd)):
		w = []
		m = []
		m.append(gA_mjd[i])
		w.append(1/(b_gA_err[i]))
		m.append(gE_mjd[i])
		w.append(1/(b_gE_err[i]))
		m.append(gi_mjd[i])
		w.append(1/(b_gi_err[i]))
		m.append(gm_mjd[i])
		w.append(1/(b_gm_err[i]))
		m.append(gq_mjd[i])
		w.append(1/(b_gq_err[i]))
		
		try:
			ms = np.average(m, weights=w)
			all_mjd.append(ms)
		except:
			all_mjd.append(gA_mjd[i])


	fbgA = fractioner(binned_gA)
	fdl_gA.append(fbgA)
	fbgE = fractioner(binned_gE)
	fdl_gE.append(fbgE)
	fbgi = fractioner(binned_gi)
	fdl_gi.append(fbgi)
	fbgm = fractioner(binned_gm)
	fdl_gm.append(fbgm)
	fbgq = fractioner(binned_gq)
	fdl_gq.append(fbgq)
	avg_fdl = (fbgA+fbgE+fbgi+fbgm+fbgq)/5
	mean_fdl.append(avg_fdl)


	#writes to a file
	with open('oct20_bin_results_HD20630_{}_{}.csv'.format(sig, bins_width), 'wt') as f:
		csv_writer = csv.writer(f, delimiter=" ")

		#csv_writer.writerow(["MDJ", "mag_gA", "err_gA", "mag_gE", "err_gE", 
			#"mag_gi", "err_gi", "mag_gm", "err_gm", "mag_gq", "err_gq"])

		for i in range(len(binned_gA[0])):
			csv_writer.writerow([all_mjd[i], b_gA[i], b_gA_err[i], binned_gA[3][i],
				b_gE[i], b_gE_err[i], binned_gE[3][i],
				b_gi[i], b_gi_err[i], binned_gi[3][i],
				b_gm[i], b_gm_err[i], binned_gm[3][i],
				b_gq[i], b_gq_err[i], binned_gq[3][i],])


	#This next section does chi squared minimization, keeping the "i" data constant.


	rangem = range(-5, 5)
	rangeq = range(-5, 5)
	rangeA = range(-5, 5)
	rangeE = range(-6, 6)
	
	facE = 0.00
	facA = -0.120
	facq = -0.026
	facm = 0.074
	

	bests = chi_func(rangem, rangeq, rangeA, rangeE,
		facE, facA, facq, facm)


	bestX = bests[0]
	bestE = bests[1]
	bestq = bests[2]
	bestm = bests[3]
	bestA = bests[4]
	

	while bestA == (max(rangeA)*0.001 + facA):
		print("Amax")
		rangeA = [x+5 for x in rangeA]

		bests = chi_func(rangem, rangeq, rangeA, rangeE,
		facE, facA, facq, facm)

		bestX = bests[0]
		bestE = bests[1]
		besta = bests[2]
		bestm = bests[3]
		bestA = bests[4]
		
	while bestE == (max(rangeE)*0.001 + facE):
		print("Emax")
		rangeE = [x+5 for x in rangeE]

		
		bests = chi_func(rangem, rangeq, rangeA, rangeE,
		facE, facA, facq, facm)

		bestX = bests[0]
		bestE = bests[1]
		besta = bests[2]
		bestm = bests[3]
		bestA = bests[4]
		
	while bestm == (max(rangem)*0.001 + facm):
		print("mmax")
		rangem = [x+5 for x in rangem]
		
		bests = chi_func(rangem, rangeq, rangeA, rangeE,
		facE, facA, facq, facm)

		bestX = bests[0]
		bestE = bests[1]
		besta = bests[2]
		bestm = bests[3]
		bestA = bests[4]
		

	while bestq == (max(rangeq)*0.001 + facq):
		print("qmax")
		rangeq = [x+5 for x in rangeq]

		bests = chi_func(rangem, rangeq, rangeA, rangeE,
		facE, facA, facq, facm)

		bestX = bests[0]
		bestE = bests[1]
		besta = bests[2]
		bestm = bests[3]
		bestA = bests[4]
		


	while bestA == (min(rangeA)*0.001 + facA):
		print("Amin")
		rangeA = [x-5 for x in rangeA]

		
		bests = chi_func(rangem, rangeq, rangeA, rangeE,
		facE, facA, facq, facm)

		bestX = bests[0]
		bestE = bests[1]
		besta = bests[2]
		bestm = bests[3]
		bestA = bests[4]
		

	while bestE == (min(rangeE)*0.001 + facE):
		print("Emin")
		rangeE = [x-5 for x in rangeE]

		bests = chi_func(rangem, rangeq, rangeA, rangeE,
		facE, facA, facq, facm)

		bestX = bests[0]
		bestE = bests[1]
		besta = bests[2]
		bestm = bests[3]
		bestA = bests[4]
		

	while bestm == (min(rangem)*0.001 + facm):
		print("mmin")
		rangem = [x-5 for x in rangem]

		bests = chi_func(rangem, rangeq, rangeA, rangeE,
		facE, facA, facq, facm)

		bestX = bests[0]
		bestE = bests[1]
		besta = bests[2]
		bestm = bests[3]
		bestA = bests[4]
		

	while bestq == (min(rangeq)*0.001 + facq):
		print("qmin")
		rangeq = [x-5 for x in rangeq]

		bests = chi_func(rangem, rangeq, rangeA, rangeE,
		facE, facA, facq, facm)

		bestX = bests[0]
		bestE = bests[1]
		besta = bests[2]
		bestm = bests[3]
		bestA = bests[4]
		

	shifts_gA.append(format(bestA, '.3f'))
	shifts_gE.append(format(bestE, '.3f'))
	shifts_gi.append("const")
	shifts_gm.append(format(bestm, '.3f'))
	shifts_gq.append(format(bestq, '.3f'))
	chis.append(bestX)


	scaled_gE = []
	scaled_gA = []
	scaled_gm = []
	scaled_gq = []
	scaled_gi = []
	
	scaled_ave = []
	scaled_ave_err = []

	for i in range(len(binned_gA[0])):
		if np.isnan(binned_gE[1][i]) or binned_gE[1][i] == 0.0:
			sE = 0.0
			scaled_gE.append(0.0)
			E_err = 999999
		else:
			sE = binned_gE[1][i]+bestE
			scaled_gE.append(sE)
			E_err = binned_gE[2][i]

		if np.isnan(binned_gi[1][i]) or binned_gi[1][i] == 0.0:
			si = 0.0
			scaled_gi.append(0.0)
			i_err = 999999
		else:
			si = binned_gi[1][i]
			scaled_gi.append(si)
			i_err = binned_gi[2][i]

		if np.isnan(binned_gm[1][i]) or binned_gm[1][i] == 0.0:
			sm = 0.0
			scaled_gm.append(0.0)
			m_err = 999999
		else:
			sm = binned_gm[1][i]+bestm
			scaled_gm.append(sm)
			m_err = binned_gm[2][i]

		if np.isnan(binned_gq[1][i]) or binned_gq[1][i] == 0.0:
			sq = 0.0
			scaled_gq.append(0.0)
			q_err = 999999
		else:
			sq = binned_gq[1][i]+bestq
			scaled_gq.append(sq)
			q_err = binned_gq[2][i]

		if np.isnan(binned_gA[1][i]) or binned_gA[1][i] == 0.0:
			sA = 0.0
			scaled_gA.append(0.0)
			A_err = 999999
		else:
			sA = binned_gA[1][i]+bestA
			scaled_gA.append(sA)
			A_err = binned_gA[2][i]

		try:
			ave = (sA/(A_err**2) + sE/(E_err**2) + si/(i_err**2) + sm/(m_err**2) + sq/(q_err**2))/(
				1/(A_err**2) + 1/(E_err**2) + 1/(i_err**2) + 1/(m_err**2) + 1/(q_err**2))
		except:
			ave = 0
		
		scaled_ave.append(ave)

		try:
			sae = 1/(1/(A_err**2) + 1/(E_err**2) + 1/(i_err**2) + 1/(m_err**2) + 1/(q_err**2))
		except:
			sae = 9999999

		scaled_ave_err.append(sae)


	with open('oct20_scaled_results_HD20630_{}_{}.csv'.format(sig, bins_width), 'wt') as f:
		csv_writer = csv.writer(f, delimiter=" ")

		csv_writer.writerow(["scaling factors: gE: {} gA: {} gm: {} gq: {} Chi square: {}".format(bestE, bestA, bestm, bestq, bestX)])

		#csv_writer.writerow(["MDJ", "mag_gA", "err_gA", "mag_gE", "err_gE", 
			#"mag_gi", "err_gi", "mag_gm", "err_gm", "mag_gq", "err_gq", "mag_Va", "err_Va", "mag_Ve", "err_Ve", "scaled avg", "uncertainty"])

		for i in range(len(binned_gA[0])):
			csv_writer.writerow([all_mjd[i], scaled_gA[i], binned_gA[2][i], binned_gA[3][i]+bestA,
				scaled_gE[i], binned_gE[2][i], binned_gE[3][i]+bestE,
				scaled_gi[i], binned_gi[2][i], binned_gi[3][i],
				scaled_gm[i], binned_gm[2][i], binned_gm[3][i]+bestm,
				scaled_gq[i], binned_gq[2][i], binned_gq[3][i]+bestq,
				scaled_ave[i], scaled_ave_err[i]])

name2 = 'oct20_tabulated_results_HD20630.csv'
with open(name2, 'wt') as f:
	csv_writer = csv.writer(f, delimiter=" ")

	csv_writer.writerow(["SigmaClip", "gEshift", "gAshift", "gmshift", "gqshift", "Chi", "fractional data loss gA", "fdl gE", "fdl gi", "fdl gm", "fdl gq", "mean fdl"])

	for i in range(len(sigmas)):
		csv_writer.writerow([sigmas[i], shifts_gE[i], shifts_gA[i], shifts_gm[i], shifts_gq[i], chis[i], fdl_gA[i], fdl_gE[i], fdl_gi[i], fdl_gm[i], fdl_gq[i], mean_fdl[i]])

