import numpy as np
import math
import csv


#E, F, j, i, m, n

data_bE = []
data_bF = []
data_bj = []
data_bi = []
data_bm = []
data_bn = []

file_1 = "lc_p_2019.csv"

with open(file_1, 'rt') as csvfile:
	data = csv.reader(csvfile, delimiter=',')
	next(data, None)
	for row in data:
		row[0] = format(float(row[0]) - 2400000.5, '.5f')
		row[6] = float(row[6])
		row[7] = float(row[7])
		if row[6] > 20:
			continue
		if row[2] == 'bE':
			data_bE.append(row)
		if row[2] == 'bF':
			data_bF.append(row)
		if row[2] == 'bj':
			data_bj.append(row)
		if row[2] == 'bi':
			data_bi.append(row)
		if row[2] == 'bm':
			data_bm.append(row)
		if row[2] == 'bn':
			data_bn.append(row)

with open("19_lcprox_bE.txt", 'wt') as csvfile:
	csv_writer = csv.writer(csvfile, delimiter=' ')
	for row in data_bE:
		csv_writer.writerow(row)

with open("19_lcprox_bF.txt", 'wt') as csvfile:
	csv_writer = csv.writer(csvfile, delimiter=' ')
	for row in data_bF:
		csv_writer.writerow(row)

with open("19_lcprox_bj.txt", 'wt') as csvfile:
	csv_writer = csv.writer(csvfile, delimiter=' ')
	for row in data_bj:
		csv_writer.writerow(row)

with open("19_lcprox_bi.txt", 'wt') as csvfile:
	csv_writer = csv.writer(csvfile, delimiter=' ')
	for row in data_bi:
		csv_writer.writerow(row)

with open("19_lcprox_bm.txt", 'wt') as csvfile:
	csv_writer = csv.writer(csvfile, delimiter=' ')
	for row in data_bm:
		csv_writer.writerow(row)

with open("19_lcprox_bn.txt", 'wt') as csvfile:
	csv_writer = csv.writer(csvfile, delimiter=' ')
	for row in data_bn:
		csv_writer.writerow(row)



