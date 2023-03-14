import numpy as np
import pandas
import copy
import scipy
from scipy import stats
import csv

    # Read Data
ifile  = open('HSE_PBE_SOC_data.csv', "rt")
reader = csv.reader(ifile)
csvdata=[]
for row in reader:
        csvdata.append(row)   
ifile.close()
numrow=len(csvdata)
numcol=len(csvdata[0]) 
csvdata = np.array(csvdata).reshape(numrow,numcol)
Formula = csvdata[:,0]
Mixing = csvdata[:,1]
HSE_decomp_en = csvdata[:,2]
HSE_gap  = csvdata[:,3]
HSE_slme = csvdata[:,4]
Comp_desc = csvdata[:,5:19]
Elem_desc = csvdata[:,19:]

X = csvdata[:,5:]

n = Formula.size
m = int(X.size/n)

Corr = [[0.0 for a in range(m)] for b in range(3)]

xx = [0.0]*n
prop = [0.0]*n

for i in range(0,m):
    for j in range(0,n):
        xx[j] = float(X[j,i])
        prop[j] = float(HSE_decomp_en[j])
    x = stats.pearsonr(xx[:],prop[:])
    Corr[0][i] = x[0]

for i in range(0,m):
    for j in range(0,n):
        xx[j] = float(X[j,i])
        prop[j] = float(HSE_gap[j])
    x = stats.pearsonr(xx[:],prop[:])
    Corr[1][i] = x[0]

for i in range(0,m):
    for j in range(0,n):
        xx[j] = float(X[j,i])
        prop[j] = float(HSE_slme[j])
    x = stats.pearsonr(xx[:],prop[:])
    Corr[2][i] = x[0]


f_out=open("test.csv",'w')
writer=csv.writer(f_out)
writer.writerows(Corr)
f_out.close()











