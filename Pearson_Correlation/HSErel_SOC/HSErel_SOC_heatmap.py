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
# number of descriptors plotted
#1-14 are composition descriptors
#15-50 are element properties descriptors
# m = 50
m = 14
# m = 36

Corr = [[0.0 for a in range(m)] for b in range(3)]

Labels = pandas.read_excel('Corr.xlsx', 'Label', engine='openpyxl')
Prop = [ 'Decomposition Energy', 'Band Gap', 'SLME @ 5μm']

# generating plot grids
# xx for grids in x axis
# y for grids in y axis, representing 4 correlated properties

x = np.arange(m)
xx = [0.0] * m
for i in range(0, m):
    xx[i] = x[i] + 0.5
y = [0.5, 1.5, 2.5]
f = 16
r = 75

for i in range(0, m):
    for j in range(0, 3):
        #use this when m=14 or 50
        Corr[j][i]  = np.float16(data[j,i])
        # use this when m=36
        # Corr[j][i] = np.float16(data[j, i + 14])

scale = ['linear']
plotposition = [131, 132, 133]


fig = plt.figure(figsize=(12, 4), dpi=450)
plt.rcParams.update({'font.size': 16})
plt.rc('font', family='Arial narrow')
plt.subplots_adjust(left=0.19, right=1.04, top=0.85, bottom=0.24, wspace=0.2, hspace=0.2)

ax = plt.plot(plotposition[0])
plt.plot(plotposition[0])
plt.xscale(scale[0])
plt.yscale(scale[0])
plt.xlim([0, m])
plt.ylim([0, 3])

# 3 different  plot settings for x ticks rotation
plt.xticks(xx[0:m], Labels.Label[0:m], rotation=90, fontsize=14)
# plt.xticks(xx[0:m], Labels.Label[0:m], rotation=75, fontsize=18)

# this plot setting is only used for elemental properties descriptor when m=36
# plt.xticks(xx[0:m], Labels.Label[14:m + 14], rotation=90, fontsize=14)

plt.yticks(y[:], Prop[:], rotation=30, fontsize=20)
plt.pcolor(Corr, cmap='seismic')
cbar = plt.colorbar(orientation="vertical", pad=0.02)
cbar.set_label(label='Correlation Coefficient', size=22)

plt.savefig('HSE_SOC_elem_frac_remake.png', dpi=450, bbox_inches='tight')
