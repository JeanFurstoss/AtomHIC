import numpy as np
import pandas as pd
import glob
import os
import os.path
import sys
from sklearn import cluster
from sklearn.mixture import GaussianMixture as GM
from sklearn.model_selection import GridSearchCV
import re
from sklearn.mixture._gaussian_mixture import _compute_precision_cholesky
import matplotlib.pyplot as plt

path = "../../tests/FitGMM/ForTest/Lab"

data = []
for i in range(5):
    data.append(np.transpose(np.genfromtxt(path+str(i+1)+"/data1")))

path_gmm = "../../tests/FitGMM/Ref_output.dat"
fp = open(path_gmm)
lines = fp.readlines()
mu = []
cov = []
dim = 2
for i in range(5):
    mu.append([])
    cov.append([])
    current_line = 0
    Found = False
    nbClust = 0
    for l in range(len(lines)):
        if "Lab"+str(i+1) in lines[l]:
            current_line = l
            Found = True
        if Found and l == current_line+2:
            nbClust = int(lines[l].split(' ')[1].strip())
    if( Found and nbClust != 0):
        for n in range(nbClust):
            mu[i].append(np.zeros(dim))
            cov[i].append(np.zeros((dim,dim)))
            for d in range(dim):
                mu[i][n][d] = float(lines[current_line+5+n*(4+dim)].split(' ')[1+d].strip())
                for d2 in range(dim):
                    cov[i][n][d][d2] = float(lines[current_line+7+n*(4+dim)+d].split(' ')[d2].strip())

colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99']

theta = np.linspace(0,2*np.pi,1000)
facval = 5.

ps = 2
ls = 2

fig = plt.figure()

for i in range(len(data)):
    # plot points
    plt.scatter(data[i][0],data[i][1],c=colors[i],s=ps)
    plt.scatter(-10,-10,c=colors[i],s=10,label="Label "+str(i+1))
    # plot ellispis
    for n in range(len(mu[i])):
        eigval, eigvec = np.linalg.eig(cov[i][n])
        ellipsis = ( np.sqrt(facval * eigval[None,:]) * eigvec ) @ [ np.sin(theta), np.cos(theta) ]
        ellipsis[0] += mu[i][n][0]
        ellipsis[1] += mu[i][n][1]
        plt.plot(ellipsis[0,:],ellipsis[1,:],c=colors[i],linewidth=ls)
plt.xlim([-0.25, 2.25])
plt.ylim([-0.25, 2.25])
plt.legend()
plt.show()
