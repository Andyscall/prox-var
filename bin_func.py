import math
import numpy as np
import matplotlib.pyplot as plt
import csv

np.seterr(all='raise')



Vbe = np.genfromtxt(open("HD20630/lcgoodV_be.dat"), names=True, delimiter=" ")

gE = np.genfromtxt(open("HD82443/lcgoodg_bE.txt"), delimiter=" ")
gi = np.genfromtxt(open("HD82443/lcgoodg_bi.txt"), delimiter=" ")
gA = np.genfromtxt(open("HD82443/lcgoodg_bA.txt"), delimiter=" ")
gm = np.genfromtxt(open("HD82443/lcgoodg_bm.txt"), delimiter=" ")
gq = np.genfromtxt(open("HD82443/lcgoodg_bq.txt"), delimiter=" ")


def bin_func(binfile, sigma=2.0, corr=0.0, bin_width=1):
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

	for i in range(len(binfile)):
		mjd_Ve.append(binfile[i][0])
		mag_Ve.append(binfile[i][5])
		mag_err_Ve.append(binfile[i][6])

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
			b_mag_err_Ve.append(99.999)
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
			mags_err_Ve.append(99.999)
			clipped_mag.append(b_mag_Ve[i])

	vals_vE = [mjds_Ve, mags_Ve, mags_err_Ve, clipped_mag]

	return(vals_vE)


def mjd_func(binfile, bin_width=1):
	mjd_Ve = []
	mag_err_Ve = []

	for i in range(len(binfile)):
		mjd_Ve.append(binfile[i][0])
		mag_err_Ve.append(binfile[i][6])

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



def chi_func(mrange = range(-5, 5), qrange = range(-5, 5), irange = range(-5, 5), Erange = range(-5, 5),
	Efac=0.092, ifac = 0.011, qfac = -0.009, mfac = 0.099):
	bestE = 0.0
	besti = 0.0
	bestm = 0.0
	bestq = 0.0
	bestX = 9.99*10**10

	for m in mrange:
		for q in qrange:
			for i in irange:
				for E in Erange:

					chisum = 0.0

					#technically these lists for the A data are redundant, but I put them in for consistency
					#The "scale_gX_w" list is for weights; it's equal to the inverse square err list.
					scale_gA = []
					scale_gA_err = []
					scale_gA_w = []

					scaleE = E*0.001 + Efac
					scale_gE = []
					scale_gE_err = []
					scale_gE_w = []

					scalei = i*0.001 + ifac
					scale_gi = []
					scale_gi_err = []
					scale_gi_w = []

					scaleq = q*0.001 + qfac
					scale_gq = []
					scale_gq_err = []
					scale_gq_w = []

					scalem = m*0.001 + mfac
					scale_gm = []
					scale_gm_err = []
					scale_gm_w = []


					for j in range(len(binned_gA[1])):
						#This loop will multiply each of the data by their scale factor.
						#Not sure if the error needs any more manipulation done to it?
						if np.isnan(binned_gA[1][j]) or binned_gA[1][j] == 0.0:
							#this is necessary so that the program doesn't try to match everything to a bunch of zeroes
							continue
							#sgA = 0
							#sgAe = 99.999
							#sgAw = 0
						else:
							sgA = binned_gA[1][j]
							sgAe = binned_gA[2][j]
							sgAw = 1/(sgAe**2)
						scale_gA.append(sgA)
						scale_gA_err.append(sgAe)
						scale_gA_w.append(sgAw)

						if np.isnan(binned_gE[1][j]) or binned_gE[1][j] == 0.0:
							sgE = 0
							sgEe = 99.999
							sgEw = 0
						else:
							sgE = binned_gE[1][j]+scaleE
							sgEe = binned_gE[2][j]
							sgEw = 1/(sgEe**2)
						scale_gE.append(sgE)
						scale_gE_err.append(sgEe)
						scale_gE_w.append(sgEw)

						if np.isnan(binned_gi[1][j]) or binned_gi[1][j] == 0.0:
							sgi = 0
							sgie =99.999
							sgiw = 0
						else:
							sgi = binned_gi[1][j]+scalei
							sgie = binned_gi[2][j]
							sgiw = 1/(sgie**2)
						scale_gi.append(sgi)
						scale_gi_err.append(sgie)
						scale_gi_w.append(sgiw)

						if np.isnan(binned_gq[1][j]) or binned_gq[1][j] == 0.0:
							sgq = 0
							sgqe = 99.999
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
							sgme =99.999
							sgmw = 0
						else:
							sgm = binned_gm[1][j]+scalem
							sgme = binned_gm[2][j]
							sgmw = 1/(sgme**2)
						scale_gm.append(sgm)
						scale_gm_err.append(sgme)
						scale_gm_w.append(sgmw)

						
					for j in range(len(binned_gA[1])):
						#Takes the weighted average of each bin; each scale_gX[j]*scale_gX_w[j] over the sum of the [j] weights
						try:
							ave = (scale_gE[j]*scale_gE_w[j] + scale_gi[j]*scale_gi_w[j] + scale_gA[j]*scale_gA_w[j] 
								+ scale_gq[j]*scale_gq_w[j] + scale_gm[j]*scale_gm_w[j]
								)/(scale_gA_w[j] + scale_gE_w[j] + scale_gi_w[j] + scale_gq_w[j] + scale_gm_w[j])
						except:
							continue
						
						chi2 = 0
						#print("e", scale_gE[j])
						#print("ave", ave)
						chiE = ((scale_gE[j] - ave)**2)/(scale_gE_err[j]**2)
						chii = ((scale_gi[j] - ave)**2)/(scale_gi_err[j]**2)
						chiA = ((scale_gA[j] - ave)**2)/(scale_gA_err[j]**2)
						
						chiq = ((scale_gq[j] - ave)**2)/(scale_gq_err[j]**2)
						chim = ((scale_gm[j] - ave)**2)/(scale_gm_err[j]**2)
						
						chi2 = chiE + chii + chiA + chiq + chim
						if chi2 == 0:
							continue
						else:
							chisum = chisum + chi2
							#print("chisum", chisum)
					
					if chisum < bestX:
						#the actual minimization bit; upates the values if necessary
						bestX = chisum
						besti = scalei
						bestE = scaleE
						bestq = scaleq
						bestm = scalem
				#print("i is", i, "bestX is", bestX)
			#print("q is", q)
		print("m is", m) #, "bestX is", bestX)

	print("Best chi:", bestX)
	print("scale E:", bestE)
	print("scale i:", besti)
	print("scale q:", bestq)
	print("scale m:", bestm)
	return [bestX, besti, bestE, bestq, bestm]


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

	binned_gA = bin_func(gA, sigma=sig, bin_width=bins_width)
	binned_gE = bin_func(gE, sigma=sig, bin_width=bins_width)
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
		if np.isnan(binned_gE[1][i]):
			b_gE.append(0)
			b_gE_err.append(99.999)
		else:
			sE = binned_gE[1][i]
			b_gE.append(sE)
			b_gE_err.append(binned_gE[2][i])

		if np.isnan(binned_gi[1][i]):
			si = 0
			b_gi.append(0)
			b_gi_err.append(99.999)
		else:
			si = binned_gi[1][i]
			b_gi.append(si)
			b_gi_err.append(binned_gi[2][i])

		if np.isnan(binned_gm[1][i]):
			sm = 0
			b_gm.append(0)
			b_gm_err.append(99.999)
		else:
			sm = binned_gm[1][i]
			b_gm.append(sm)
			b_gm_err.append(binned_gm[2][i])

		if np.isnan(binned_gq[1][i]):
			sq = 0
			b_gq.append(0)
			b_gq_err.append(99.999)
		else:
			sq = binned_gq[1][i]
			b_gq.append(sq)
			b_gq_err.append(binned_gq[2][i])

		if np.isnan(binned_gA[1][i]):
			sA = 0
			b_gA.append(0)
			b_gA_err.append(99.999)
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



	file_name = 'oct20_bin_results_HD82443_{}_{}.tsv'.format(sig, bins_width)

	with open(file_name, 'wt') as f:
		tsv_writer = csv.writer(f, delimiter="\t")

		tsv_writer.writerow(["MDJ", "mag_gA", "err_gA", "cut_gA", "mag_gE", "err_gE", "cut_gE",
			"mag_gi", "err_gi", "cut_gi", "mag_gm", "err_gm", "cut_gm", "mag_gq", "err_gq", "cut_gq"])

		for i in range(len(binned_gA[0])):
			tsv_writer.writerow([format(all_mjd[i], '.3f'), format(b_gA[i], '.3f'), format(b_gA_err[i], '.3f'), format(binned_gA[3][i], '.3f'),
				format(b_gE[i], '.3f'), format(b_gE_err[i], '.3f'), format(binned_gE[3][i], '.3f'),
				format(b_gi[i], '.3f'), format(b_gi_err[i], '.3f'), format(binned_gi[3][i], '.3f'),
				format(b_gm[i], '.3f'), format(b_gm_err[i], '.3f'), format(binned_gm[3][i], '.3f'),
				format(b_gq[i], '.3f'), format(b_gq_err[i], '.3f'), format(binned_gq[3][i], '.3f')])


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

	"""
	This next section does chi squared minimization, keeping the "A" data constant.
	"A" may not be the best choice, since it only covers some of the time range, but it looks prettiest
	in the range it does cover, so I chose it. Should be fairly easy to manually switch it out for one of the others.
	"""

	rangem = range(-5, 5)
	rangeq = range(-5, 5)
	rangei = range(-5, 5)
	rangeE = range(-5, 5)

	facE = 0.092
	faci = 0.011
	facq = -0.009
	facm = 0.099


	bests = chi_func(rangem, rangeq, rangei, rangeE,
		facE, faci, facq, facm)

	bestX = bests[0]
	besti = bests[1]
	bestE = bests[2]
	bestq = bests[3]
	bestm = bests[4]
	#print(besti, bestE, bestq, bestm)
	#print(min(rangem))
	#print(min(rangem) + facm)

	while besti == (max(rangei)*0.001 + faci):
		print("imax")
		rangei = [x+5 for x in rangei]

		bests = chi_func(rangem, rangeq, rangei, rangeE,
		facE, faci, facq, facm)

		bestX = bests[0]
		besti = bests[1]
		bestE = bests[2]
		bestq = bests[3]
		bestm = bests[4]

	while bestE == (max(rangeE)*0.001 + facE):
		print("Emax")
		rangeE = [x+5 for x in rangeE]

		bests = chi_func(rangem, rangeq, rangei, rangeE,
		facE, faci, facq, facm)

		bestX = bests[0]
		besti = bests[1]
		bestE = bests[2]
		bestq = bests[3]
		bestm = bests[4]

	while bestm == (max(rangem)*0.001 + facm):
		print("mmax")
		rangem = [x+5 for x in rangem]

		bests = chi_func(rangem, rangeq, rangei, rangeE,
		facE, faci, facq, facm)

		bestX = bests[0]
		besti = bests[1]
		bestE = bests[2]
		bestq = bests[3]
		bestm = bests[4]

	while bestq == (max(rangeq)*0.001 + facq):
		print("qmax")
		rangeq = [x+5 for x in rangeq]

		bests = chi_func(rangem, rangeq, rangei, rangeE,
		facE, faci, facq, facm)

		bestX = bests[0]
		besti = bests[1]
		bestE = bests[2]
		bestq = bests[3]
		bestm = bests[4]



	while besti == (min(rangei)*0.001 + faci):
		print("imin")
		rangei = [x-5 for x in rangei]

		bests = chi_func(rangem, rangeq, rangei, rangeE,
		facE, faci, facq, facm)

		bestX = bests[0]
		besti = bests[1]
		bestE = bests[2]
		bestq = bests[3]
		bestm = bests[4]

	while bestE == (min(rangeE)*0.001 + facE):
		print("Emin")
		rangeE = [x-5 for x in rangeE]

		bests = chi_func(rangem, rangeq, rangei, rangeE,
		facE, faci, facq, facm)

		bestX = bests[0]
		besti = bests[1]
		bestE = bests[2]
		bestq = bests[3]
		bestm = bests[4]

	while bestm == (min(rangem)*0.001 + facm):
		print("mmin")
		rangem = [x-5 for x in rangem]

		bests = chi_func(rangem, rangeq, rangei, rangeE,
		facE, faci, facq, facm)

		bestX = bests[0]
		besti = bests[1]
		bestE = bests[2]
		bestq = bests[3]
		bestm = bests[4]

	while bestq == (min(rangeq)*0.001 + facq):
		print("qmin")
		rangeq = [x-5 for x in rangeq]

		bests = chi_func(rangem, rangeq, rangei, rangeE,
		facE, faci, facq, facm)

		bestX = bests[0]
		besti = bests[1]
		bestE = bests[2]
		bestq = bests[3]
		bestm = bests[4]

	shifts_gA.append("const")
	shifts_gE.append(format(bestE, '.3f'))
	shifts_gi.append(format(besti, '.3f'))
	shifts_gm.append(format(bestm, '.3f'))
	shifts_gq.append(format(bestq, '.3f'))
	chis.append(bestX)


	scaled_gE = []
	scaled_gi = []
	scaled_gm = []
	scaled_gq = []
	scaled_gA = []

	scaled_ave = []
	scaled_ave_err = []

	for i in range(len(binned_gA[0])):
		if np.isnan(binned_gE[1][i]) or binned_gE[1][i] == 0.0:
			sE = 0.0
			scaled_gE.append(0.0)
			E_err = 99.999
		else:
			sE = binned_gE[1][i]+bestE
			scaled_gE.append(sE)
			E_err = binned_gE[2][i]

		if np.isnan(binned_gi[1][i]) or binned_gi[1][i] == 0.0:
			si = 0.0
			scaled_gi.append(0.0)
			i_err = 99.999
		else:
			si = binned_gi[1][i]+besti
			scaled_gi.append(si)
			i_err = binned_gi[2][i]

		if np.isnan(binned_gm[1][i]) or binned_gm[1][i] == 0.0:
			sm = 0.0
			scaled_gm.append(0.0)
			m_err = 99.999
		else:
			sm = binned_gm[1][i]+bestm
			scaled_gm.append(sm)
			m_err = binned_gm[2][i]

		if np.isnan(binned_gq[1][i]) or binned_gq[1][i] == 0.0:
			sq = 0.0
			scaled_gq.append(0.0)
			q_err = 99.999
		else:
			sq = binned_gq[1][i]+bestq
			scaled_gq.append(sq)
			q_err = binned_gq[2][i]

		if np.isnan(binned_gA[1][i]) or binned_gA[1][i] == 0.0:
			sA = 0.0
			scaled_gA.append(0.0)
			A_err = 99.999
		else:
			sA = binned_gA[1][i]
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
			sae = 99.999

		scaled_ave_err.append(sae)

	name = 'oct20_scaled_results_HD82443_{}_{}.tsv'.format(sig, bins_width)

	with open(name, 'wt') as f:
		tsv_writer = csv.writer(f, delimiter="\t")

		tsv_writer.writerow(["scaling factors: gA: const gE: {} gi: {} gm: {} gq: {} Chi square: {}".format(bestE, besti, bestm, bestq, bestX)])

		tsv_writer.writerow(["MDJ", "mag_gA", "err_gA", "cut_gA", "mag_gE", "err_gE", "cut_gE",
			"mag_gi", "err_gi", "cut_gi", "mag_gm", "err_gm", "cut_gm", "mag_gq", "err_gq", "cut_gq", 
			"scaled avg", "uncertainty"])
		
		for i in range(len(binned_gA[0])):
			tsv_writer.writerow([format(all_mjd[i], '.3f'), 
				format(scaled_gA[i], '.3f'), format(binned_gA[2][i], '.3f'), format(binned_gA[3][i], '.3f'),
				format(scaled_gE[i], '.3f'), format(binned_gE[2][i], '.3f'), format(binned_gE[3][i]+bestE, '.3f'),
				format(scaled_gi[i], '.3f'), format(binned_gi[2][i], '.3f'), format(binned_gi[3][i]+besti, '.3f'),
				format(scaled_gm[i], '.3f'), format(binned_gm[2][i], '.3f'), format(binned_gm[3][i]+bestm, '.3f'),
				format(scaled_gq[i], '.3f'), format(binned_gq[2][i], '.3f'), format(binned_gq[3][i]+bestq, '.3f'),
				format(scaled_ave[i], '.3f'), format(scaled_ave_err[i], '.3f')])


name2 = 'oct20_tabulated_results_HD82443.tsv'
with open(name2, 'wt') as f:
	tsv_writer = csv.writer(f, delimiter="\t")
	tsv_writer.writerow(["SigmaClip", "gEshift", "gishift", "gmshift", "gqshift", "Chi", "fdl gA", "fdl gE", "fdl gi", "fdl gm", "fdl gq", "mean fdl"])

	for i in range(len(sigmas)):
		tsv_writer.writerow([sigmas[i], shifts_gE[i], shifts_gi[i], 
			shifts_gm[i], shifts_gq[i], chis[i], fdl_gA[i], 
			fdl_gE[i], fdl_gi[i], fdl_gm[i], fdl_gq[i], mean_fdl[i]])
