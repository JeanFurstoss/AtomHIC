from sklearn import datasets
from sklearn import decomposition
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import sys
import argparse
import os.path
from os import path

pca = PCA(n_components=2)

parser = argparse.ArgumentParser()
parser.add_argument("CrystalName", help="Name of the crystal")
parser.add_argument("ave_cutoff or not", help="0 if the PCA is performed on Steinhardt params or 1 if the PCA is performed on averaged Steinhardt parameters")

crystalName=sys.argv[1]
AC=int(sys.argv[2])

fullpath="../../data/Steinhardt/"+crystalName+"/PCA_Analysis/"

if(path.exists(fullpath)):
    Defect_names = []
    St_Params = []
    RefNames = []
    # read the defects in the database
    for root, subdirs, files in os.walk(fullpath):
        for filename in files:
            if( filename[0] != "." ):
                if( AC == 1 and filename[-14:] == "ave_cutoff.dat" ):
                    Defect_names.append(filename)
                if( AC == 0 and filename[-14:] != "ave_cutoff.dat" ):
                    Defect_names.append(filename)
    for i in range(len(Defect_names)):
        file_path = fullpath+Defect_names[i]
        St_Params.append([])
        RefNames.append([])
        f = open(file_path, 'r')
        Lines = f.readlines()
        count = 0
        for line in Lines:
            if( count == 1 ):
                l_sph=int(line.split()[1])
            elif( count == 3 ):
                n_ref=int(line.split()[1])
                for r in range(n_ref):
                    St_Params[i].append([])
                    RefNames[i].append("")
            elif( count > 3 ):
                curRef = int(line.split()[1])-1
                RefNames[i][curRef] = line.split()[0]
                for r in range(l_sph):
                    St_Params[i][curRef].append(float(line.split()[r+2]))
            count += 1
        f.close()

    datas = []
    for r in range(n_ref):
        datas.append([])
        for i in range(len(Defect_names)):
            datas[r].append(np.array(St_Params[i][r]))
        #print(datas[r])
        pca_dat = pca.fit_transform(datas[r])
        for i in range(len(pca_dat)):
            plt.scatter(pca_dat[i][0],pca_dat[i][1], s=200) 
            plt.text(pca_dat[i][0],pca_dat[i][1], Defect_names[i][0:-15])
        plt.title(RefNames[0][r])
        plt.savefig(RefNames[0][r]+"_PCA.png")
        plt.clf()

else:
    print("The directory \""+fullpath+"\" does not exist")

