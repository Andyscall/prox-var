import math
import numpy as np
import matplotlib.pyplot as plt
import csv

#goodg = np.genfromtxt(open("proxima/lcgoodcombinedg.txt"), names=True, delimiter=" ")

gi = []
gj = []
gE = []
gF = []
gn = []
gm = []


with open("proxima/lcgoodcombinedg.txt", 'rb') as csvfile:
	data = csv.reader(csvfile, delimiter=' ')

	for row in data:
		if row[2] == "bi":
			gi.append(row)
		elif row[2] == "bj":
			gj.append(row)
		elif row[2] == "bE":
			gE.append(row)
		elif row[2] == "bF":
			gF.append(row)
		elif row[2] == "bn":
			gn.append(row)
		elif row[2] == "bm":
			gm.append(row)
		else:
			print("problem row", row)



with open('lcgoodg_bi.txt', 'wt') as f:
	csv_writer = csv.writer(f, delimiter=" ")

	for i in gi:
		csv_writer.writerow(i)


with open('lcgoodg_bj.txt', 'wt') as f:
	csv_writer = csv.writer(f, delimiter=" ")

	for i in gj:
		csv_writer.writerow(i)


with open('lcgoodg_bE.txt', 'wt') as f:
	csv_writer = csv.writer(f, delimiter=" ")

	for i in gE:
		csv_writer.writerow(i)


with open('lcgoodg_bF.txt', 'wt') as f:
	csv_writer = csv.writer(f, delimiter=" ")

	for i in gF:
		csv_writer.writerow(i)


with open('lcgoodg_bn.txt', 'wt') as f:
	csv_writer = csv.writer(f, delimiter=" ")

	for i in gn:
		csv_writer.writerow(i)


with open('lcgoodg_bm.txt', 'wt') as f:
	csv_writer = csv.writer(f, delimiter=" ")

	for i in gm:
		csv_writer.writerow(i)



