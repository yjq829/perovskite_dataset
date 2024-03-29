import numpy as np
import pandas
import copy
import scipy
from scipy import stats
import csv

    # Read Data
ifile  = open('PBE_data.csv', "rt")
reader = csv.reader(ifile)
csvdata=[]
for row in reader:
        csvdata.append(row)   
ifile.close()
numrow=len(csvdata)
numcol=len(csvdata[0]) 
csvdata = np.array(csvdata).reshape(numrow,numcol)

Formula = csvdata[0:514,0]
Mixing = csvdata[0:514,1]
PBE_latt = csvdata[0:514,2]
PBE_gap  = csvdata[0:514,3]
PBE_form_en = csvdata[0:514,4]
PBE_decomp_en = csvdata[0:514,5]
PBE_slme5  = csvdata[0:514,6]
PBE_fom  = csvdata[0:514,7]
Comp_desc = csvdata[0:514,8:22]
Elem_desc = csvdata[0:514,22:]


X = csvdata[:,8:]

n = Formula.size
m = int(X.size/n)
print(m)
Corr = [[0.0 for a in range(m)] for b in range(4)]

xx = [0.0]*n
prop = [0.0]*n

for i in range(0,m):
    for j in range(0,n):
        xx[j] = float(X[j,i])
        prop[j] = float(PBE_latt[j])
    x = stats.pearsonr(xx[:],prop[:])
    Corr[0][i] = x[0]


for i in range(0,m):
    for j in range(0,n):
        xx[j] = float(X[j,i])
        prop[j] = float(PBE_decomp_en[j])
    x = stats.pearsonr(xx[:],prop[:])
    Corr[1][i] = x[0]

for i in range(0,m):
    for j in range(0,n):
        xx[j] = float(X[j,i])
        prop[j] = float(PBE_gap[j])
    x = stats.pearsonr(xx[:],prop[:])
    Corr[2][i] = x[0]



for i in range(0,m):
    for j in range(0,n):
        xx[j] = float(X[j,i])
        prop[j] = float(PBE_slme5[j])
    x = stats.pearsonr(xx[:],prop[:])
    Corr[3][i] = x[0]


print(len(Corr[3]))


f_out=open("test.csv",'w')
writer=csv.writer(f_out)
writer.writerows(Corr)
f_out.close()







