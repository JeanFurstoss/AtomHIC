#!/usr/bin/python3

from joblib import Parallel, delayed
import numpy as np
import glob
import os
import os.path
import sys
from sklearn.mixture import GaussianMixture as GM
from sklearn.model_selection import GridSearchCV
import re

#python3 ComputeProbability.py PathToDirectoryContainingBinaryFilesOfGMMModel
dirModel = os.path.join(sys.argv[1], '')

#Find out whether previous training was successful
trained = False
if os.path.exists(dirModel+"training.out"):
    text_file = open(dirModel+"training.out", "r")
    for line in text_file:
        if re.search('All clusters have been properly recognized', line):
            trained = True
            print('Reading previously trained model')
            print()
    text_file.close()

if( trained ):
    #Get number of clusters
    text_file = open(dirModel+"training.out", "r")
    for line in text_file:
        if re.search('Clusters Mg ---  ', line):
            print(line)
            numClustersMg = int(line.split('---  ')[1])
        elif re.search('Clusters Si ---  ', line):
            print(line)
            numClustersSi = int(line.split('---  ')[1])
        elif re.search('Clusters Ox ---  ', line):
            print(line)
            numClustersOx = int(line.split('---  ')[1])
        elif re.search('Categories: ', line):
            categories = line.split('[', 1)[1].split(']')[0].replace("'", "").split(', ')
    text_file.close()

    #Read the values
    meansMg = np.load(dirModel+'meansMg.npy')
    meansSi = np.load(dirModel+'meansSi.npy')
    meansOx = np.load(dirModel+'meansOx.npy')
    covarMg = np.load(dirModel+'covariancesMg.npy')
    covarSi = np.load(dirModel+'covariancesSi.npy')
    covarOx = np.load(dirModel+'covariancesOx.npy')
    weightsMg = np.load(dirModel+'weightsMg.npy')
    weightsSi = np.load(dirModel+'weightsSi.npy')
    weightsOx = np.load(dirModel+'weightsOx.npy')
    
    #Read clusters id's from training output
    clusterIDsMg = ['']*numClustersMg
    clusterIDsSi = ['']*numClustersSi
    clusterIDsOx = ['']*numClustersOx
    text_file = open(dirModel+"training.out", "r")
    for line in text_file:
        if re.search('Identified clusters Mg: ', line):
            clusterIDsMg = line.split('[', 1)[1].split(']')[0].replace("'", "").split(', ')
        elif re.search('Identified clusters Si: ', line):
            clusterIDsSi = line.split('[', 1)[1].split(']')[0].replace("'", "").split(', ')
        elif re.search('Identified clusters Ox: ', line):
            clusterIDsOx = line.split('[', 1)[1].split(']')[0].replace("'", "").split(', ')
    text_file.close()
    
    structures = np.unique(clusterIDsMg)
    dim = len(meansMg[0])
    
    for i in range(len(structures)):
        f = open(structures[i]+'.dat','w')
        f.write("NUMBER_OF_DIMENSION "+str(dim)+"\n")
        f.write("ATOM_TYPE Mg\n")
        nbclust = 0
        for j in range(len(clusterIDsMg)):
            if( clusterIDsMg[j] == structures[i] ):
                nbclust += 1
        f.write("NUMBER_OF_CLUSTER "+str(nbclust)+"\n")
        for j in range(len(clusterIDsMg)):
            if( clusterIDsMg[j] == structures[i] ):
                f.write("WEIGHT "+str(weightsMg[j])+"\n")
                f.write("DETERMINANT "+str(np.linalg.det(covarMg[j]))+"\n")
                f.write("ESPERANCE")
                for n in range(dim):
                    f.write(" "+str(meansMg[j][n]))
                f.write("\n")
                inv_cov = np.linalg.inv(covarMg[j])
                f.write("INVERSE OF COVARIANT MATRIX\n")
                for n1 in range(dim):
                    for n2 in range(dim):
                        f.write(str(inv_cov[n1][n2])+" ")
                    f.write("\n")
        f.write("ATOM_TYPE Si\n")
        nbclust = 0
        for j in range(len(clusterIDsSi)):
            if( clusterIDsSi[j] == structures[i] ):
                nbclust += 1
        f.write("NUMBER_OF_CLUSTER "+str(nbclust)+"\n")
        for j in range(len(clusterIDsSi)):
            if( clusterIDsSi[j] == structures[i] ):
                f.write("WEIGHT "+str(weightsSi[j])+"\n")
                f.write("DETERMINANT "+str(np.linalg.det(covarSi[j]))+"\n")
                f.write("ESPERANCE")
                for n in range(dim):
                    f.write(" "+str(meansSi[j][n]))
                f.write("\n")
                inv_cov = np.linalg.inv(covarSi[j])
                f.write("INVERSE OF COVARIANT MATRIX\n")
                for n1 in range(dim):
                    for n2 in range(dim):
                        f.write(str(inv_cov[n1][n2])+" ")
                    f.write("\n")
        f.write("ATOM_TYPE O\n")
        nbclust = 0
        for j in range(len(clusterIDsOx)):
            if( clusterIDsOx[j] == structures[i] ):
                nbclust += 1
        f.write("NUMBER_OF_CLUSTER "+str(nbclust)+"\n")
        for j in range(len(clusterIDsOx)):
            if( clusterIDsOx[j] == structures[i] ):
                f.write("WEIGHT "+str(weightsOx[j])+"\n")
                f.write("DETERMINANT "+str(np.linalg.det(covarOx[j]))+"\n")
                f.write("ESPERANCE")
                for n in range(dim):
                    f.write(" "+str(meansOx[j][n]))
                f.write("\n")
                inv_cov = np.linalg.inv(covarOx[j])
                f.write("INVERSE OF COVARIANT MATRIX\n")
                for n1 in range(dim):
                    for n2 in range(dim):
                        f.write(str(inv_cov[n1][n2])+" ")
                    f.write("\n")


