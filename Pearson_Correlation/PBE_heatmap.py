import numpy as np
import csv
import copy
import random
import matplotlib.pyplot as plt
import pandas
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Read Data

ifile = open('test.csv', "rt")
reader = csv.reader(ifile)
csvdata = []
for row in reader:
    csvdata.append(row)
ifile.close()
numrow = len(csvdata)
numcol = len(csvdata[0])
csvdata = np.array(csvdata).reshape(numrow, numcol)
data = csvdata[:, :]

# m = 50
# m = 14
m = 36

Corr = [[0.0 for a in range(m)] for b in range(4)]

Labels = pandas.read_excel('Corr.xlsx', 'Label')
# print(Labels)
Prop = ['Lattice Constant', 'Decomposition Energy', 'Band Gap','SLME @ 5μm']

x = np.arange(m)
xx = [0.0] * m
for i in range(0, m):
    xx[i] = x[i] + 0.5
y = [0.5, 1.5, 2.5, 3.5]
f = 16
r = 75

for i in range(0, m):
    for j in range(0, 4):
        # Corr[j][i] = np.float16(data[j, i])
       Corr[j][i]  = np.float16(data[j,i+14])

scale = ['linear']
plotposition = [131, 132, 133]

# fig=plt.figure(figsize=(14,4),dpi=450)
# plt.rcParams.update({'font.size': 16})
# plt.rc('font', family='Arial narrow')
# plt.subplots_adjust(left=0.18, right=1.06, top=0.85, bottom=0.23, wspace=0.2, hspace=0.2)

# fig=plt.figure(figsize=(8,4),dpi=450)
# plt.rcParams.update({'font.size': 16})
# plt.rc('font', family='Arial narrow')
# plt.subplots_adjust(left=0.31, right=0.97, top=0.85, bottom=0.24, wspace=0.2, hspace=0.2)

fig = plt.figure(figsize=(12, 4), dpi=450)
plt.rcParams.update({'font.size': 16})
plt.rc('font', family='Arial narrow')
plt.subplots_adjust(left=0.21, right=1.04, top=0.85, bottom=0.24, wspace=0.2, hspace=0.2)

ax = plt.plot(plotposition[0])
plt.plot(plotposition[0])
plt.xscale(scale[0])
plt.yscale(scale[0])
plt.xlim([0, m])
plt.ylim([0, 4])

# plt.xticks(xx[0:m], Labels.Label[0:m], rotation=75, fontsize=16)
# plt.xticks(xx[0:m], Labels.Label[0:m], rotation=90, fontsize=16)
plt.xticks(xx[0:m], Labels.Label[14:m+14], rotation=90, fontsize=16)

plt.yticks(y[:], Prop[:], rotation=30, fontsize=20)
# plt.title('Correlation: Descriptors vs PBE Properties', fontname='Arial narrow', size=24, horizontalalignment='center',
#           pad=15)
plt.pcolor(Corr, cmap='seismic',vmax=1.0, vmin=-1.0)
cbar = plt.colorbar(orientation="vertical", pad=0.02)
cbar.set_label(label='Correlation Coefficient', size=22)

plt.savefig('PBE_prop_remake.png', dpi=450)
