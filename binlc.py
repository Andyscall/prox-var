import math
import numpy as np
import matplotlib.pyplot as plt


lcgVbe = np.genfromtxt(open("lcgoodV_be.dat"), names=True, delimiter=" ")
lcgVba = np.genfromtxt(open("lcgoodV_ba.dat"), names=True, delimiter=" ")
lcggbq = np.genfromtxt(open("lcgoodg_bq.dat"), names=True, delimiter=" ")
lcggbm = np.genfromtxt(open("lcgoodg_bm.dat"), names=True, delimiter=" ")
lcggbi = np.genfromtxt(open("lcgoodg_bi.dat"), names=True, delimiter=" ")
lcggbE = np.genfromtxt(open("lcgoodg_bE.dat"), names=True, delimiter=" ")
lcggbA = np.genfromtxt(open("lcgoodg_bA.dat"), names=True, delimiter=" ")


mjd_Ve = []
mag_Ve = []
mag_err_Ve = [] 

mjd_Va = []
mag_Va = []
mag_err_Va = []

mjd_gq = []
mag_gq = []
mag_err_gq = []

mjd_gm = []
mag_gm = []
mag_err_gm = []

mjd_gi = []
mag_gi = []
mag_err_gi = []

mjd_gE = []
mag_gE = []
mag_err_gE = []

mjd_gA = []
mag_gA = []
mag_err_gA = []

for i in range(len(lcgVbe)):
	mjd_Ve.append(lcgVbe[i][0])
	mag_Ve.append(lcgVbe[i][5])
	mag_err_Ve.append(lcgVbe[i][6])


for i in range(len(lcgVba)):
	mjd_Va.append(lcgVba[i][0])
	mag_Va.append(lcgVba[i][5])
	mag_err_Va.append(lcgVba[i][6])

for i in range(len(lcggbq)):
	mjd_gq.append(lcggbq[i][0])
	mag_gq.append(lcggbq[i][5])
	mag_err_gq.append(lcggbq[i][6])

for i in range(len(lcggbm)):
	mjd_gm.append(lcggbm[i][0])
	mag_gm.append(lcggbm[i][5])
	mag_err_gm.append(lcggbm[i][6])

for i in range(len(lcggbi)):
	mjd_gi.append(lcggbi[i][0])
	mag_gi.append(lcggbi[i][5])
	mag_err_gi.append(lcggbi[i][6])

for i in range(len(lcggbE)):
	mjd_gE.append(lcggbE[i][0])
	mag_gE.append(lcggbE[i][5])
	mag_err_gE.append(lcggbE[i][6])

for i in range(len(lcggbA)):
	mjd_gA.append(lcggbA[i][0])
	mag_gA.append(lcggbA[i][5])
	mag_err_gA.append(lcggbA[i][6])


b_mjd_Ve = []
b_mag_Ve = []
b_mag_err_Ve = []

b_mjd_Va = []
b_mag_Va = []
b_mag_err_Va = []

b_mjd_gq = []
b_mag_gq = []
b_mag_err_gq = []

b_mjd_gm = []
b_mag_gm = []
b_mag_err_gm = []

b_mjd_gi = []
b_mag_gi = []
b_mag_err_gi = []

b_mjd_gE = []
b_mag_gE = []
b_mag_err_gE = []

b_mjd_gA = []
b_mag_gA = []
b_mag_err_gA = []


for i in range(math.floor(mjd_Ve[0]*100), math.floor(mjd_Ve[len(mjd_Ve)-1]*100)):
	b_i = []
	b_e_i_p = []
	for j in range(len(mjd_Ve)):
		if math.floor(mjd_Ve[j]*100)==i:
			b_i.append(mag_Ve[j])
			b_e_i_p.append((mag_err_Ve[j]/mag_Ve[j])**2)
	if len(b_i)!=0:
		b_mjd_Ve.append(i*0.01)
		b_mag_Ve.append(np.mean(b_i))
		b_mag_err_Ve.append(np.mean(b_i)*((sum(b_e_i_p))**(1/2)))

for i in range(math.floor(mjd_Va[0]*100), math.floor(mjd_Va[len(mjd_Va)-1]*100)):
	b_i = []
	b_e_i_p = []
	for j in range(len(mjd_Va)):
		if math.floor(mjd_Va[j]*100)==i:
			b_i.append(mag_Va[j])
			b_e_i_p.append((mag_err_Va[j]/mag_Va[j])**2)
	if len(b_i)!=0:
		b_mjd_Va.append(i*0.01)
		b_mag_Va.append(np.mean(b_i))
		b_mag_err_Va.append(np.mean(b_i)*((sum(b_e_i_p))**(1/2)))


for i in range(math.floor(mjd_gq[0]*100), math.floor(mjd_gq[len(mjd_gq)-1]*100)):
	b_i = []
	b_e_i_p = []
	for j in range(len(mjd_gq)):
		if math.floor(mjd_gq[j]*100)==i:
			b_i.append(mag_gq[j])
			b_e_i_p.append((mag_err_gq[j]/mag_gq[j])**2)
	if len(b_i)!=0:
		b_mjd_gq.append(i*0.01)
		b_mag_gq.append(np.mean(b_i))
		b_mag_err_gq.append(np.mean(b_i)*((sum(b_e_i_p))**(1/2)))

for i in range(math.floor(mjd_gm[0]*100), math.floor(mjd_gm[len(mjd_gm)-1]*100)):
	b_i = []
	b_e_i_p = []
	for j in range(len(mjd_gm)):
		if math.floor(mjd_gm[j]*100)==i:
			b_i.append(mag_gm[j])
			b_e_i_p.append((mag_err_gm[j]/mag_gm[j])**2)
	if len(b_i)!=0:
		b_mjd_gm.append(i*0.01)
		b_mag_gm.append(np.mean(b_i))
		b_mag_err_gm.append(np.mean(b_i)*((sum(b_e_i_p))**(1/2)))

for i in range(math.floor(mjd_gi[0]*100), math.floor(mjd_gi[len(mjd_gi)-1]*100)):
	b_i = []
	b_e_i_p = []
	for j in range(len(mjd_gi)):
		if math.floor(mjd_gi[j]*100)==i:
			b_i.append(mag_gi[j])
			b_e_i_p.append((mag_err_gi[j]/mag_gi[j])**2)
	if len(b_i)!=0:
		b_mjd_gi.append(i*0.01)
		b_mag_gi.append(np.mean(b_i))
		b_mag_err_gi.append(np.mean(b_i)*((sum(b_e_i_p))**(1/2)))

for i in range(math.floor(mjd_gE[0]*100), math.floor(mjd_gE[len(mjd_gE)-1]*100)):
	b_i = []
	b_e_i_p = []
	for j in range(len(mjd_gE)):
		if math.floor(mjd_gE[j]*100)==i:
			b_i.append(mag_gE[j])
			b_e_i_p.append((mag_err_gE[j]/mag_gE[j])**2)
	if len(b_i)!=0:
		b_mjd_gE.append(i*0.01)
		b_mag_gE.append(np.mean(b_i))
		b_mag_err_gE.append(np.mean(b_i)*((sum(b_e_i_p))**(1/2)))

for i in range(math.floor(mjd_gA[0]*100), math.floor(mjd_gA[len(mjd_gA)-1]*100)):
	b_i = []
	b_e_i_p = []
	for j in range(len(mjd_gA)):
		if math.floor(mjd_gA[j]*100)==i:
			b_i.append(mag_gA[j])
			b_e_i_p.append((mag_err_gA[j]/mag_gA[j])**2)
	if len(b_i)!=0:
		b_mjd_gA.append(i*0.01)
		b_mag_gA.append(np.mean(b_i))
		b_mag_err_gA.append(np.mean(b_i)*((sum(b_e_i_p))**(1/2)))



mmag_Ve = np.mean(b_mag_Ve)
varm_Ve = np.var(b_mag_Ve)
stdm_Ve = np.std(b_mag_Ve)
mjds_Ve = []
mags_Ve = []

mmag_Va = np.mean(b_mag_Va)
varm_Va = np.var(b_mag_Va)
stdm_Va = np.std(b_mag_Va)
mjds_Va = []
mags_Va = []

mmag_gq = np.mean(b_mag_gq)
varm_gq = np.var(b_mag_gq)
stdm_gq = np.std(b_mag_gq)
mjds_gq = []
mags_gq = []

mmag_gm = np.mean(b_mag_gm)
varm_gm = np.var(b_mag_gm)
stdm_gm = np.std(b_mag_gm)
mjds_gm = []
mags_gm = []

mmag_gi = np.mean(b_mag_gi)
varm_gi = np.var(b_mag_gi)
stdm_gi = np.std(b_mag_gi)
mjds_gi = []
mags_gi = []

mmag_gE = np.mean(b_mag_gE)
varm_gE = np.var(b_mag_gE)
stdm_gE = np.std(b_mag_gE)
mjds_gE = []
mags_gE = []

mmag_gA = np.mean(b_mag_gA)
varm_gA = np.var(b_mag_gA)
stdm_gA = np.std(b_mag_gA)
mjds_gA = []
mags_gA = []


for i in range(len(b_mag_Ve)):
	if (mmag_Ve-2*stdm_Ve) <= b_mag_Ve[i] <= (mmag_Ve+2*stdm_Ve):
		mjds_Ve.append(b_mjd_Ve[i])
		mags_Ve.append(b_mag_Ve[i])

plt.scatter(mjds_Ve, mags_Ve, marker="+")

for i in range(len(b_mag_Va)):
	if (mmag_Va-2*stdm_Va) <= b_mag_Va[i] <= (mmag_Va+2*stdm_Va):
		mjds_Va.append(b_mjd_Va[i])
		mags_Va.append(b_mag_Va[i])

plt.scatter(mjds_Va, mags_Va, marker="+")

for i in range(len(b_mag_gq)):
	if (mmag_gq-2*stdm_gq) <= b_mag_gq[i] <= (mmag_gq+2*stdm_gq):
		mjds_gq.append(b_mjd_gq[i])
		mags_gq.append(b_mag_gq[i])

plt.scatter(mjds_gq, mags_gq, marker="+")

for i in range(len(b_mag_gm)):
	if (mmag_gm-2*stdm_gm) <= b_mag_gm[i] <= (mmag_gm+2*stdm_gm):
		mjds_gm.append(b_mjd_gm[i])
		mags_gm.append(b_mag_gm[i])

plt.scatter(mjds_gm, mags_gm, marker="+")

for i in range(len(b_mag_gi)):
	if (mmag_gi-2*stdm_gi) <= b_mag_gi[i] <= (mmag_gi+2*stdm_gi):
		mjds_gi.append(b_mjd_gi[i])
		mags_gi.append(b_mag_gi[i])

plt.scatter(mjds_gi, mags_gi, marker="+")

for i in range(len(b_mag_gE)):
	if (mmag_gE-2*stdm_gE) <= b_mag_gE[i] <= (mmag_gE+2*stdm_gE):
		mjds_gE.append(b_mjd_gE[i])
		mags_gE.append(b_mag_gE[i])

plt.scatter(mjds_gE, mags_gE, marker="+")

for i in range(len(b_mag_gA)):
	if (mmag_gA-2*stdm_gA) <= b_mag_gA[i] <= (mmag_gA+2*stdm_gA):
		mjds_gA.append(b_mjd_gA[i])
		mags_gA.append(b_mag_gA[i])

plt.scatter(mjds_gA, mags_gA, marker="+")
