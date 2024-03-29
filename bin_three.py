import math
import numpy as np
import matplotlib.pyplot as plt
import csv

np.seterr(all='raise')



gG = np.genfromtxt(open("gj358/lcgoodg_bG.txt"), names=True, delimiter=" ")
gH = np.genfromtxt(open("gj358/lcgoodg_bH.txt"), names=True, delimiter=" ")
gl = np.genfromtxt(open("gj358/lcgoodg_bl.txt"), names=True, delimiter=" ")
go = np.genfromtxt(open("gj358/lcgoodg_bo.txt"), names=True, delimiter=" ")
gp = np.genfromtxt(open("gj358/lcgoodg_bp.txt"), names=True, delimiter=" ")


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

	for i in np.arange(58000, 59000, bin_width):
		
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
			b_mjd_means.append(0.0)
			b_mag_Ve.append(0.0)
			b_mag_err_Ve.append(9999)
			b_mag_err_Ve_w.append(0.0)

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
			clipped_mag.append(0.0)

		else:
			#adds an error-filled zero, so as to not have empty data points. Not sure this works?
			mjds_Ve.append(b_mjd_Ve[i])
			mags_Ve.append(0.0)
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

	for i in np.arange(58000, 59000, bin_width):
		
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

def chi_func(orange = range(-5, 5), prange = range(-5, 5), Hrange = range(-5, 5), Grange = range(-5, 5),
	Gfac = 0.0, Hfac = 0.0, pfac = 0.0, ofac = 0.0):
	bestG = 0.0
	bestl = 0.0
	besto = 0.0
	bestp = 0.0
	bestH = 0.0
	bestX = 9.99*10**10


	for H in Hrange:
		for o in orange:
			for p in prange:
				for G in Grange:
					
					chisum = 0.0

					#technically these lists for the I data are redundant, but I put them in for consistency
					#The "scale_XX_w" list is for weights; it's equal to the inverse square err list.
					scalel = 0
					scale_gl = []
					scale_gl_err = []
					scale_gl_w = []

					scaleG = G*0.001 + Gfac
					scale_gG = []
					scale_gG_err = []
					scale_gG_w = []

					scalep = p*0.001 + pfac
					scale_gp = []
					scale_gp_err = []
					scale_gp_w = []

					scaleo = o*0.001 + ofac
					scale_go = []
					scale_go_err = []
					scale_go_w = []

					scaleH = H*0.001 + Hfac
					scale_gH = []
					scale_gH_err = []
					scale_gH_w = []

					for j in range(len(binned_gl[1])):
						#This loop will multiply each of the data by their scale factor.
						#Not sure if the error needs any more manipulation done to it?
						if np.isnan(binned_gl[1][j]) or binned_gl[1][j] == 0.0:
							#this is necessary so that the program doesn't try to match everything to a bunch of zeroes
							continue
						else:
							sgl = binned_gl[1][j]
							sgle = binned_gl[2][j]
							sglw = 1/(sgle**2)
						scale_gl.append(sgl)
						scale_gl_err.append(sgle)
						scale_gl_w.append(sglw)

						if np.isnan(binned_gH[1][j]) or binned_gH[1][j] == 0.0:
							sgH = 0.0
							sgHe = 9999999
							sgHw = 0.0
						else:
							sgH = binned_gH[1][j]+scaleH
							sgHe = binned_gH[2][j]
							sgHw = 1/(sgHe**2)
						scale_gH.append(sgH)
						scale_gH_err.append(sgHe)
						scale_gH_w.append(sgHw)

						if np.isnan(binned_gG[1][j]) or binned_gG[1][j] == 0.0:
							sgG = 0.0
							sgGe = 9999999
							sgGw = 0.0
						else:
							sgG = binned_gG[1][j]+scaleG
							sgGe = binned_gG[2][j]
							sgGw = 1/(sgGe**2)
						scale_gG.append(sgG)
						scale_gG_err.append(sgGe)
						scale_gG_w.append(sgGw)
						
						if np.isnan(binned_gp[1][j]) or binned_gp[1][j] == 0.0:
							sgp = 0.0
							sgpe = 9999999
							sgpw = 0.0
						else:
							sgp = binned_gp[1][j]+scalep
							sgpe = binned_gp[2][j]
							sgpw = 1/(sgpe**2)
						scale_gp.append(sgp)
						scale_gp_err.append(sgpe)
						scale_gp_w.append(sgpw)

						if np.isnan(binned_go[1][j]) or binned_go[1][j] == 0.0:
							sgo = 0.0
							sgoe = 9999999
							sgow = 0.0
						else:
							sgo = binned_go[1][j]+scaleo
							sgoe = binned_go[2][j]
							sgow = 1/(sgoe**2)
						scale_go.append(sgo)
						scale_go_err.append(sgoe)
						scale_go_w.append(sgow)

					for j in range(len(binned_gl[1])):
						#Takes the weighted average of each bin; each scale_xX[j]*scale_xX_w[j] over the sum of the [j] weights
						try:
							ave = (scale_gG[j]*scale_gG_w[j] + scale_gl[j]*scale_gl_w[j] + scale_gH[j]*scale_gH_w[j] 
								+ scale_gp[j]*scale_gp_w[j] + scale_go[j]*scale_go_w[j]
								)/(scale_gG_w[j] + scale_gl_w[j] + scale_gp_w[j] + scale_go_w[j]
								+ scale_gH_w[j])
						except:
							continue
						
						chi2 = 0
						chiG = ((scale_gG[j] - ave)**2)/(scale_gG_err[j]**2)
						chil = ((scale_gl[j] - ave)**2)/(scale_gl_err[j]**2)
						chiH = ((scale_gH[j] - ave)**2)/(scale_gH_err[j]**2)
						chip = ((scale_gp[j] - ave)**2)/(scale_gp_err[j]**2)
						chio = ((scale_go[j] - ave)**2)/(scale_go_err[j]**2)
						
						chi2 = chiG + chil + chip + chio + chiH
						chisum = chisum + chi2
					
					if chisum < bestX:
						#the actual minimization bit; upates the values if necessary
						bestX = chisum
						bestH = scaleH
						bestG = scaleG
						bestp = scalep
						besto = scaleo
						

				#print("q is", q)#, "bestX is", bestX)
			#print("m is", m)
		print("H is", H)

	print("Best chi:", bestX)
	print("scale G:", bestG)
	print("scale p:", bestp)
	print("scale o:", besto)
	print("scale H:", bestH)
	
	return [bestX, bestG, bestp, besto, bestH]


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


sigmas = [1.0, 1.5, 2.0, 2.5, 3.0]
#sigmas = [2.5]
bins_width = 1

shifts_gH = []
shifts_gG = []
shifts_gl = []
shifts_go = []
shifts_gp = []
chis = []

fdl_gH = []
fdl_gG = []
fdl_gl = []
fdl_go = []
fdl_gp = []
mean_fdl = []


for sig in sigmas:

	binned_gG = bin_func(gG, sigma=sig, bin_width=bins_width)
	binned_gH = bin_func(gH, sigma=sig, bin_width=bins_width)
	binned_gl = bin_func(gl, sigma=sig, bin_width=bins_width)
	binned_go = bin_func(go, sigma=sig, bin_width=bins_width)
	binned_gp = bin_func(gp, sigma=sig, bin_width=bins_width)
	

	gH_mjd = mjd_func(gH, bin_width=bins_width)
	gG_mjd = mjd_func(gG, bin_width=bins_width)
	gl_mjd = mjd_func(gl, bin_width=bins_width)
	go_mjd = mjd_func(go, bin_width=bins_width)
	gp_mjd = mjd_func(gp, bin_width=bins_width)
	
	all_mjd = []


	b_gH = []
	b_gH_err = []
	b_gG = []
	b_gG_err = []
	b_gl = []
	b_gl_err = []
	b_go = []
	b_go_err = []
	b_gp = []
	b_gp_err = []


	for i in range(len(binned_gH[0])):
		if np.isnan(binned_gG[1][i]) or binned_gG[1][i] == 0.0:
			b_gG.append(0.0)
			b_gG_err.append(999999)
		else:
			sG = binned_gG[1][i]
			b_gG.append(sG)
			b_gG_err.append(binned_gG[2][i])

		if np.isnan(binned_gl[1][i]) or binned_gl[1][i] == 0.0:
			sl = 0.0
			b_gl.append(0.0)
			b_gl_err.append(999999)
		else:
			sl = binned_gl[1][i]
			b_gl.append(sl)
			b_gl_err.append(binned_gl[2][i])

		if np.isnan(binned_go[1][i]) or binned_go[1][i] == 0.0:
			so = 0.0
			b_go.append(0.0)
			b_go_err.append(999999)
		else:
			so = binned_go[1][i]
			b_go.append(so)
			b_go_err.append(binned_go[2][i])

		if np.isnan(binned_gp[1][i]) or binned_gp[1][i] == 0.0:
			sp = 0.0
			b_gp.append(0.0)
			b_gp_err.append(999999)
		else:
			sp = binned_gp[1][i]
			b_gp.append(sp)
			b_gp_err.append(binned_gp[2][i])

		if np.isnan(binned_gH[1][i]) or binned_gH[1][i] == 0.0:
			sH = 0.0
			b_gH.append(0.0)
			b_gH_err.append(999999)
		else:
			sH = binned_gH[1][i]
			b_gH.append(sH)
			b_gH_err.append(binned_gH[2][i])


	for i in range(len(gH_mjd)):
		w = []
		m = []
		m.append(gH_mjd[i])
		w.append(1/(b_gH_err[i]))
		m.append(gG_mjd[i])
		w.append(1/(b_gG_err[i]))
		m.append(gl_mjd[i])
		w.append(1/(b_gl_err[i]))
		m.append(go_mjd[i])
		w.append(1/(b_go_err[i]))
		m.append(gp_mjd[i])
		w.append(1/(b_gp_err[i]))
		
		try:
			ms = np.average(m, weights=w)
			all_mjd.append(ms)
		except:
			all_mjd.append(gH_mjd[i])


	fbgH = fractioner(binned_gH)
	fdl_gH.append(fbgH)
	fbgG = fractioner(binned_gG)
	fdl_gG.append(fbgG)
	fbgl = fractioner(binned_gl)
	fdl_gl.append(fbgl)
	fbgo = fractioner(binned_go)
	fdl_go.append(fbgo)
	fbgp = fractioner(binned_gp)
	fdl_gp.append(fbgp)
	avg_fdl = (fbgH+fbgG+fbgl+fbgo+fbgp)/5
	mean_fdl.append(avg_fdl)


	#writes to a file
	with open('bin_results_gj358_{}_{}.csv'.format(sig, bins_width), 'wt') as f:
		csv_writer = csv.writer(f, delimiter=" ")

		#csv_writer.writerow(["MDJ", "mag_gA", "err_gA", "mag_gE", "err_gE", 
			#"mag_gi", "err_gi", "mag_gm", "err_gm", "mag_gq", "err_gq"])

		for i in range(len(binned_gH[0])):
			csv_writer.writerow([all_mjd[i], b_gH[i], b_gH_err[i], binned_gH[3][i],
				b_gG[i], b_gG_err[i], binned_gG[3][i],
				b_gl[i], b_gl_err[i], binned_gl[3][i],
				b_go[i], b_go_err[i], binned_go[3][i],
				b_gp[i], b_gp_err[i], binned_gp[3][i],])


	#This next section does chi squared minimization, keeping the "i" data constant.


	rangeo = range(-5, 5)
	rangep = range(-5, 5)
	rangeH = range(-5, 5)
	rangeG = range(-5, 5)
	
	facG = -0.029
	facH = 0.00
	facp = 0.00
	faco = 0.009
	

	bests = chi_func(rangeo, rangep, rangeH, rangeG,
		facG, facH, facp, faco)


	bestX = bests[0]
	bestG = bests[1]
	bestp = bests[2]
	besto = bests[3]
	bestH = bests[4]
	

	while bestH == (max(rangeH)*0.001 + facH):
		print("Hmax")
		rangeH = [x+4 for x in rangeH]
		
		bests = chi_func(rangeo, rangep, rangeH, rangeG, facG, facH, facp, faco)

		bestX = bests[0]
		bestG = bests[1]
		bestp = bests[2]
		besto = bests[3]
		bestH = bests[4]
		
	while bestG == (max(rangeG)*0.001 + facG):
		print("Gmax")
		rangeG = [x+4 for x in rangeG]

		bests = chi_func(rangeo, rangep, rangeH, rangeG, facG, facH, facp, faco)

		bestX = bests[0]
		bestG = bests[1]
		bestp = bests[2]
		besto = bests[3]
		bestH = bests[4]
		
	while besto == (max(rangeo)*0.001 + faco):
		print("omax")
		rangeo = [x+4 for x in rangeo]
		
		bests = chi_func(rangeo, rangep, rangeH, rangeG, facG, facH, facp, faco)

		bestX = bests[0]
		bestG = bests[1]
		bestp = bests[2]
		besto = bests[3]
		bestH = bests[4]
		
	while bestp == (max(rangep)*0.001 + facp):
		print("pmax")
		rangep = [x+4 for x in rangep]

		bests = chi_func(rangeo, rangep, rangeH, rangeG, facG, facH, facp, faco)

		bestX = bests[0]
		bestG = bests[1]
		bestp = bests[2]
		besto = bests[3]
		bestH = bests[4]
		


	while bestH == (min(rangeH)*0.001 + facH):
		print("Hmin")
		rangeH = [x-4 for x in rangeH]

		bests = chi_func(rangeo, rangep, rangeH, rangeG, facG, facH, facp, faco)

		bestX = bests[0]
		bestG = bests[1]
		bestp = bests[2]
		besto = bests[3]
		bestH = bests[4]
		
	while bestG == (min(rangeG)*0.001 + facG):
		print("Gmin")
		rangeG = [x-4 for x in rangeG]

		bests = chi_func(rangeo, rangep, rangeH, rangeG, facG, facH, facp, faco)

		bestX = bests[0]
		bestG = bests[1]
		bestp = bests[2]
		besto = bests[3]
		bestH = bests[4]
		
	while besto == (min(rangeo)*0.001 + faco):
		print("omin")
		rangeo = [x-4 for x in rangeo]

		bests = chi_func(rangeo, rangep, rangeH, rangeG, facG, facH, facp, faco)

		bestX = bests[0]
		bestG = bests[1]
		bestp = bests[2]
		besto = bests[3]
		bestH = bests[4]
		
	while bestp == (min(rangep)*0.001 + facp):
		print("pmin")
		print("before", bestp)
		rangep = [x-4 for x in rangep]

		bests = chi_func(rangeo, rangep, rangeH, rangeG, facG, facH, facp, faco)
		#bests = chi_func(rangeo, rangep, rangeH, rangeG, bestG, bestH, bestp, besto)
		#bests = chi_func(rangeo, rangep, rangeH, rangeG, facG, facH, facp, faco)


		bestX = bests[0]
		bestG = bests[1]
		bestp = bests[2]
		besto = bests[3]
		bestH = bests[4]

		print("after", bestp)
		


	shifts_gH.append(format(bestH, '.3f'))
	shifts_gG.append(format(bestG, '.3f'))
	shifts_gl.append("const")
	shifts_go.append(format(besto, '.3f'))
	shifts_gp.append(format(bestp, '.3f'))
	chis.append(bestX)


	scaled_gG = []
	scaled_gH = []
	scaled_go = []
	scaled_gp = []
	scaled_gl = []
	
	scaled_ave = []
	scaled_ave_err = []

	for i in range(len(binned_gH[0])):
		if np.isnan(binned_gG[1][i]) or binned_gG[1][i] == 0.0:
			sG = 0.0
			scaled_gG.append(0.0)
			G_err = 999999
		else:
			sG = binned_gG[1][i]+bestG
			scaled_gG.append(sG)
			G_err = binned_gG[2][i]

		if np.isnan(binned_gl[1][i]) or binned_gl[1][i] == 0.0:
			sl = 0.0
			scaled_gl.append(0.0)
			l_err = 999999
		else:
			sl = binned_gl[1][i]
			scaled_gl.append(sl)
			l_err = binned_gl[2][i]

		if np.isnan(binned_go[1][i]) or binned_go[1][i] == 0.0:
			so = 0.0
			scaled_go.append(0.0)
			o_err = 999999
		else:
			so = binned_go[1][i]+besto
			scaled_go.append(so)
			o_err = binned_go[2][i]

		if np.isnan(binned_gp[1][i]) or binned_gp[1][i] == 0.0:
			sp = 0.0
			scaled_gp.append(0.0)
			p_err = 999999
		else:
			sp = binned_gp[1][i]+bestp
			scaled_gp.append(sp)
			p_err = binned_gp[2][i]

		if np.isnan(binned_gH[1][i]) or binned_gH[1][i] == 0.0:
			sH = 0.0
			scaled_gH.append(0.0)
			H_err = 999999
		else:
			sA = binned_gH[1][i]+bestH
			scaled_gH.append(sH)
			H_err = binned_gH[2][i]

		try:
			ave = (sH/(H_err**2) + sG/(G_err**2) + sl/(l_err**2) + so/(o_err**2) + sp/(p_err**2))/(
				1/(H_err**2) + 1/(G_err**2) + 1/(l_err**2) + 1/(o_err**2) + 1/(p_err**2))
		except:
			ave = 0
		
		scaled_ave.append(ave)

		try:
			sae = 1/(1/(H_err**2) + 1/(G_err**2) + 1/(l_err**2) + 1/(o_err**2) + 1/(p_err**2))
		except:
			sae = 9999999

		scaled_ave_err.append(sae)


	with open('scaled_results_gj358_{}_{}_oo.csv'.format(sig, bins_width), 'wt') as f:
		csv_writer = csv.writer(f, delimiter=" ")

		csv_writer.writerow(["scaling factors: gG: {} gH: {} go: {} gp: {} Chi square: {}".format(bestG, bestH, besto, bestp, bestX)])

		#csv_writer.writerow(["MDJ", "mag_gA", "err_gA", "mag_gE", "err_gE", 
			#"mag_gi", "err_gi", "mag_gm", "err_gm", "mag_gq", "err_gq", "mag_Va", "err_Va", "mag_Ve", "err_Ve", "scaled avg", "uncertainty"])

		for i in range(len(binned_gH[0])):
			csv_writer.writerow([all_mjd[i], scaled_gH[i], binned_gH[2][i], binned_gH[3][i]+bestH,
				scaled_gG[i], binned_gG[2][i], binned_gG[3][i]+bestG,
				scaled_gl[i], binned_gl[2][i], binned_gl[3][i],
				scaled_go[i], binned_go[2][i], binned_go[3][i]+besto,
				scaled_gp[i], binned_gp[2][i], binned_gp[3][i]+bestp,
				scaled_ave[i], scaled_ave_err[i]])

name2 = 'tabulated_results_gj358_oo.csv'
with open(name2, 'wt') as f:
	csv_writer = csv.writer(f, delimiter=" ")

	csv_writer.writerow(["SigmaClip", "gGshift", "gHshift", "goshift", "gpshift", "Chi", "fractional data loss gH", "fdl gG", "fdl gl", "fdl go", "fdl gp", "mean fdl"])

	for i in range(len(sigmas)):
		csv_writer.writerow([sigmas[i], shifts_gG[i], shifts_gH[i], shifts_go[i], shifts_gp[i], chis[i], fdl_gH[i], fdl_gG[i], fdl_gl[i], fdl_go[i], fdl_gp[i], mean_fdl[i]])

