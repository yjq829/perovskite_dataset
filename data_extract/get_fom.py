from __future__ import print_function
import numpy as np    
import csv
import copy
import random
import matplotlib.pyplot as plt
import pandas
import os
from numpy import linspace
from scipy.interpolate import UnivariateSpline



ifile  = open('absorption.dat', "rt")
reader = csv.reader(ifile)
csvdata=[]
for row in reader:
        csvdata.append(row)   
ifile.close()
numrow=len(csvdata)
numcol=len(csvdata[0]) 
csvdata = np.array(csvdata).reshape(numrow,numcol)
abs_data = csvdata[1:,0]

n = len(abs_data)
energy = [0.0]*n
lambdas = [0.0]*n
alphas  = [0.0]*n

for i in range(0,n):
    xx = abs_data[i].split(" ")
    energy[i] = float(xx[0])
    if energy[i] != 0.0:
        lambdas[i] = 1239.84193/energy[i]
    alphas[i] = (float(xx[1]) + float(xx[2]) + float(xx[3]))/3



AM_data = pandas.read_excel('Data_AM.xlsx')

m = len(AM_data.Wavelength[:])
xs = [0.0]*m
for i in range(0,m):
    xs[i] = 1239.84193/AM_data.Wavelength[i]

s1 = UnivariateSpline(energy[:], alphas[:], s=0)
y1 = s1(xs)


num = 0
den = 0
for i in range(0,m-1):
    num = num + y1[i]*AM_data.Etr[i]*(xs[i+1]-xs[i])
    den = den + AM_data.Etr[i]*(xs[i+1]-xs[i])

fom = num/den

text_file = open("fom.txt", "w")
text_file.write("fom: %s" % fom)
text_file.close()



