import numpy as np
import math
import csv


#E, F, j, i, m, n

data_bm = []
data_bA = []
data_bq = []

data_bi = []
data_bm = []
data_bn = []

data_Vbe = []
data_Vbf = []

#file_1 = "lc_p_2019.csv"
#file_1 = "Proxima/lc_p_9_2020.csv"
file_1 = "HD82443/hd82443_2020_data.csv"

with open(file_1, 'rt') as csvfile:
	data = csv.reader(csvfile, delimiter=',')
	next(data, None)
	for row in data:
		if len(row) < 5:
			continue
		row[0] = format(float(row[0]) - 2400000.5, '.5f')
		row[5] = float(row[5])
		row[6] = float(row[6])
		if row[6] > 20:
			continue
		if row[2] == 'bm':
			data_bm.append(row)
		if row[2] == 'bA':
			data_bA.append(row)
		if row[2] == 'bq':
			data_bq.append(row)
		'''
		if row[2] == 'bi':
			data_bi.append(row)
		if row[2] == 'bm':
			data_bm.append(row)
		if row[2] == 'bn':
			data_bn.append(row)
		if row[2] == 'be':
			data_Vbe.append([row[0], row[5], row[6]])
		if row[2] == 'bf':
			data_Vbf.append([row[0], row[5], row[6]])
		'''

with open("20_oct_HD82443_bm.txt", 'wt') as csvfile:
	csv_writer = csv.writer(csvfile, delimiter=' ')
	for row in data_bm:
		csv_writer.writerow(row)

with open("20_oct_HD82443_bA.txt", 'wt') as csvfile:
	csv_writer = csv.writer(csvfile, delimiter=' ')
	for row in data_bA:
		csv_writer.writerow(row)

with open("20_oct_HD82443_bq.txt", 'wt') as csvfile:
	csv_writer = csv.writer(csvfile, delimiter=' ')
	for row in data_bq:
		csv_writer.writerow(row)
"""

with open("20_sep_lcprox_bF.txt", 'wt') as csvfile:
	csv_writer = csv.writer(csvfile, delimiter=' ')
	for row in data_bF:
		csv_writer.writerow(row)

with open("20_sep_lcprox_bj.txt", 'wt') as csvfile:
	csv_writer = csv.writer(csvfile, delimiter=' ')
	for row in data_bj:
		csv_writer.writerow(row)

with open("20_sep_lcprox_bi.txt", 'wt') as csvfile:
	csv_writer = csv.writer(csvfile, delimiter=' ')
	for row in data_bi:
		csv_writer.writerow(row)

with open("20_sep_lcprox_bm.txt", 'wt') as csvfile:
	csv_writer = csv.writer(csvfile, delimiter=' ')
	for row in data_bm:
		csv_writer.writerow(row)

with open("20_sep_lcprox_bn.txt", 'wt') as csvfile:
	csv_writer = csv.writer(csvfile, delimiter=' ')
	for row in data_bn:
		csv_writer.writerow(row)
"""
"""

with open("lcprox_Vbe.txt", 'wt') as csvfile:
	csv_writer = csv.writer(csvfile, delimiter=' ')
	for row in data_Vbe:
		csv_writer.writerow(row)

with open("lcprox_Vbf.txt", 'wt') as csvfile:
	csv_writer = csv.writer(csvfile, delimiter=' ')
	for row in data_Vbf:
		csv_writer.writerow(row)
"""
