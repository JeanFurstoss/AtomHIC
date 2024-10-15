#!/usr/bin/python3
from joblib import Parallel, delayed
import numpy as np
from numpy import linalg as la
import glob
import os
import os.path
import sys
from sklearn.mixture import GaussianMixture as GM
from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt
import seaborn
import re
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.collections import PolyCollection
from matplotlib.ticker import MaxNLocator, FuncFormatter, MultipleLocator
import matplotlib.patches as mpatches
import random
import math
import argparse

dirModel = '../../data/MachineLearningModels/GaussianMixtureModel/'
database_present = glob.glob(dirModel+"/*/")
help_print = ''
for i in range(len(database_present)):
    help_print += database_present[i].replace(dirModel,'')+'\n'

parser = argparse.ArgumentParser()
parser.add_argument("AtomHICDatabaseName", help="Name of the database (already stored in AtomHIC /data/), available basis :\n"+help_print)
parser.add_argument("DumpsDatabasePath", help="Path of the directory containing the dumps used for the fitting of the GMM")
parser.add_argument("FilteringType", help="Filtering type (it should be one of the field in the .cfg dumps file)")
parser.add_argument("DescriptorName", help="Name of the descriptor (it should be one of the field in the .cfg dumps file)")
parser.add_argument("PathWhereSaveFigs", help="Path of the directory where to save figures")
parser.parse_args()

ath_database_name = sys.argv[1]
dump_database_path = sys.argv[2]
FilterType = sys.argv[3]
DescriptorName = sys.argv[4]
path_to_plot = sys.argv[5]


NbDimEachPlot = 3 # Number of dimensions of each plot
nsamp = 25000 # number of sampled points in the GMM

def readQ(input_file, FilterValue, nbDim, FilterType=FilterType, DescriptorName=DescriptorName):
    startCol = 0
    PosFilter = 0
    filter_val = []
    with open(input_file,'r') as file:
        found = False
        for line in file:
            if 'ITEM: ATOMS' in line:
                colvals = line.split(' ')
                DesFound = False
                for i in range(len(colvals)):
                    s = ''.join(x for x in colvals[i] if x.isdigit())
                    for_des = colvals[i].replace("["+s+"]",'')
                    if for_des == DescriptorName and not DesFound:
                        startCol = i-2
                        DesFound = True
                    if colvals[i] == FilterType :
                        PosFilter = i-2
                found = True
            elif found:
                filter_val.append((line.split(' ')[PosFilter]).strip())

    pos_tmp=np.genfromtxt(input_file,skip_header=9)
    Descriptors = []
    for f in range(len(FilterValue)):
        Descriptors.append([])
        Descriptors[-1] = [pos_tmp[i,startCol:startCol+nbDim] for i in range(len(pos_tmp)) if filter_val[i] == FilterValue[f] ]
    return Descriptors 

print("Reading "+ath_database_name+" AtomHIC database")

os.makedirs(dirModel, exist_ok=True)
ModelName = os.path.join(ath_database_name, '')

if not os.path.exists(dirModel+ModelName):
    print("The database does not exists in /data/MachineLearningModels/GaussianMixtureModel/")
    print("Possible databases :")
    print(help_print)

Labels = glob.glob(dirModel+ModelName+"/*.ath")
LabelsName = []
nbLabels = len(Labels)
for i in range(nbLabels):
    tmp = Labels[i].replace(dirModel+ModelName,'')
    LabelsName.append(tmp.replace('.ath',''))

# Open the first one to count how many filter there are
FilterValues = []
NbDim = 0
with open(Labels[0],'r') as file: 
    for line in file:
        if 'FILTER_VALUE' in line:
            FilterValues.append((line.split(' ')[1]).strip())
        if 'NUMBER_OF_DIMENSION' in line:
            NbDim = int((line.split(' ')[1]).strip())

nbFilter = len(FilterValues)
# Read properties of GMM
numClusters = np.zeros(nbFilter)
means = []
invcovars = []
weights = []
ClusterLabels = []
for i in range(nbFilter):
    means.append([])
    invcovars.append([])
    weights.append([])
    ClusterLabels.append([])

clust_count_fil = []
for i in range(nbFilter):
    clust_count_fil.append(0)


for i in range(len(Labels)):
    InputFile = open(Labels[i],'r')
    Lines = InputFile.readlines()
    count = 0
    countNbClust = 1000
    currentFilter = 1000
    clust_count = 1000
    for line in Lines:
        if 'FILTER_VALUE' in line:
            for l in range(nbFilter):
                if FilterValues[l] == (line.split(' ')[1]).strip():
                    currentFilter = l
                    break
            countNbClust = 1000
            clust_count = 1000
        if 'NUMBER_OF_CLUSTER' in line:
            currentNbClust = int((line.split(' ')[1]).strip())
            countNbClust = count
            for c in range(currentNbClust):
                means[currentFilter].append([])
                invcovars[currentFilter].append([])
            clust_count = 0
        if count == (countNbClust+1)+clust_count*(4+NbDim):
            weights[currentFilter].append(float((line.split(' ')[1]).strip()))
            ClusterLabels[currentFilter].append(LabelsName[i])
        if count == (countNbClust+3)+clust_count*(4+NbDim):
            for dim in range(NbDim):
                means[currentFilter][clust_count_fil[currentFilter]].append(float((line.split(' ')[1+dim]).strip()))
        if count >= (countNbClust+5)+clust_count*(4+NbDim) and count <= (countNbClust+5)+clust_count*(4+NbDim)+NbDim-1 :
            invcovars[currentFilter][clust_count_fil[currentFilter]].append([])
            for dim in range(NbDim):
                invcovars[currentFilter][clust_count_fil[currentFilter]][-1].append(float((line.split(' ')[dim]).strip()))
            if( count == (countNbClust+5)+clust_count*(4+NbDim)+NbDim-1 ):
                clust_count_fil[currentFilter] += 1
                clust_count += 1
        count += 1

print("Done")

print("Labelled GMM properties:")
for f in range(nbFilter):
    print("\tFilter : "+str(FilterValues[f]))
    for l in range(nbLabels):
        nbClust = 0
        for k in range(len(weights[f])):
            if ClusterLabels[f][k] == LabelsName[l]:
                nbClust += 1
        print("\t\t"+LabelsName[l]+": "+str(int(nbClust))+" clusters")

def nearestPD(A):
    """Find the nearest positive-definite matrix to input

    A Python/Numpy port of John D'Errico's `nearestSPD` MATLAB code [1], which
    credits [2].

    [1] https://www.mathworks.com/matlabcentral/fileexchange/42885-nearestspd

    [2] N.J. Higham, "Computing a nearest symmetric positive semidefinite
    matrix" (1988): https://doi.org/10.1016/0024-3795(88)90223-6
    """

    B = (A + A.T) / 2
    _, s, V = la.svd(B)

    H = np.dot(V.T, np.dot(np.diag(s), V))

    A2 = (B + H) / 2

    A3 = (A2 + A2.T) / 2

    if isPD(A3):
        return A3

    spacing = np.spacing(la.norm(A))
    # The above is different from [1]. It appears that MATLAB's `chol` Cholesky
    # decomposition will accept matrixes with exactly 0-eigenvalue, whereas
    # Numpy's will not. So where [1] uses `eps(mineig)` (where `eps` is Matlab
    # for `np.spacing`), we use the above definition. CAVEAT: our `spacing`
    # will be much larger than [1]'s `eps(mineig)`, since `mineig` is usually on
    # the order of 1e-16, and `eps(1e-16)` is on the order of 1e-34, whereas
    # `spacing` will, for Gaussian random matrixes of small dimension, be on
    # othe order of 1e-16. In practice, both ways converge, as the unit test
    # below suggests.
    I = np.eye(A.shape[0])
    k = 1
    while not isPD(A3):
        mineig = np.min(np.real(la.eigvals(A3)))
        A3 += I * (-mineig * k**2 + spacing)
        k += 1

    return A3


def isPD(B):
    """Returns true when input is positive-definite, via Cholesky"""
    try:
        _ = la.cholesky(B)
        return True
    except la.LinAlgError:
        return False


#Creating and load the GM model
GMMModels = []
for i in range(nbFilter):
    for k in range(len(invcovars[i])):
        temp = nearestPD(np.array(invcovars[i][k]))
        invcovars[i][k] = temp
    GMMModels.append(GM(len(weights[i])))
    GMMModels[-1].weights_ = np.array(weights[i])
    GMMModels[-1].means_ = np.array(means[i])
    GMMModels[-1].precisions_cholesky_ = np.linalg.cholesky(np.array(invcovars[i]))
    GMMModels[-1].covariances_ = np.linalg.inv(np.array(invcovars[i]))


templist = set(ClusterLabels[0])
categories = (list(templist))
categories = sorted(categories)

#Get corresponding array index for probabilities
dictLab = []
for i in range(nbFilter):
    dictLab.append([None]*len(categories))
    for j in range(len(categories)):
        dictLab[i][j] = []
        for b in range(len(ClusterLabels[i])):
            if ClusterLabels[i][b] == categories[j]:
                dictLab[i][j] += [b]

shifting_plot = 1./(nbLabels+1)

q_art = []
lab = []
q_GMM = []
for l in range(nbFilter):
    q_art_temp, lab_temp = GMMModels[l].sample(nsamp)
    q_art.append(q_art_temp)
    lab.append(lab_temp)
    q_GMM.append([])
    for i in range(len(categories)):
        q_GMM[-1].append([])

for i in range(nsamp):
    found = []
    cat = []
    for l in range(nbFilter):
        found.append(False)
        cat.append(-1)
    for b in range(len(categories)):
        for l in range(nbFilter):
            if not found[l]:
                for s in range(len(dictLab[l][b])):
                    if( lab[l][i] == dictLab[l][b][s] ):
                        cat[l] = b
                        found[l] = True
                        break
        all_found = True
        for l in range(nbFilter):
            all_found *= found[l]
        if all_found:
            break
    for l in range(nbFilter):
        if found[l]:
            q_GMM[l][cat[l]].append(q_art[l][i])

NbVp = int(NbDim/NbDimEachPlot)
NbVpComp = math.ceil((NbDim/NbDimEachPlot)-int(NbDim/NbDimEachPlot))
NbLastL = NbDim-(NbVp*NbDimEachPlot)
LastLoopBeg = NbVp*NbDimEachPlot
LastLoopEnd = NbVp*NbDimEachPlot+NbLastL

dat_vp = []
for i in range(NbVp+NbVpComp):
    dat_vp.append([])
    for l in range(nbFilter):
        dat_vp[i].append([])

print("Reading dumps in "+dump_database_path+" directory")

for base, dirs, files in os.walk(dump_database_path):
    for directories in dirs:
        index_cat = -1
        for i in range(len(categories)):
            if( directories == categories[i] ):
                index_cat = i
                break
        list_file = sorted(glob.glob(base+directories+"/*.cfg"))
        # read first file to estimate number of point and reduce if necessary
        MeanNbDes = np.zeros(len(FilterValues))
        for input_file in list_file:
            Descriptors = readQ(input_file, FilterValues, NbDim)
            for i in range(len(FilterValues)):
                MeanNbDes[i] += len(Descriptors[i])
        
        MinNbDes = np.min(MeanNbDes)
        MinNbDes /= len(list_file)
        if( len(list_file)*MinNbDes > nsamp ):
            step_to_rem = round(nsamp/MinNbDes)
            list_file_tmp = []
            for f in range(step_to_rem):
                list_file_tmp.append(list_file[round((len(list_file)-1)/step_to_rem)])
            list_file = []
            for f in range(step_to_rem):
                list_file.append(list_file_tmp[f])

        for input_file in list_file:
            Descriptors = readQ(input_file, FilterValues, NbDim)
            for l in range(nbFilter):
                for q_l in range(len(Descriptors[l])):
                    for v in range(NbVp):
                        for i in range(v*NbDimEachPlot,(v+1)*NbDimEachPlot):
                            randind = random.randint(0,len(q_GMM[l][index_cat])-1)
                            dat_vp[v][l].append(["database",i+1+(index_cat+1)*shifting_plot,Descriptors[l][q_l][i]])
                            dat_vp[v][l].append(["GMM",i+1+(index_cat+1)*shifting_plot,q_GMM[l][index_cat][randind][i]])
                    for i in range(LastLoopBeg,LastLoopEnd):
                        randind = random.randint(0,len(q_GMM[l][index_cat])-1)
                        dat_vp[-1][l].append(["database",i+1+(index_cat+1)*shifting_plot,Descriptors[l][q_l][i]])
                        dat_vp[-1][l].append(["GMM",i+1+(index_cat+1)*shifting_plot,q_GMM[l][index_cat][randind][i]])


color_panel = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'] # from colorbrewer2.org
colors_struct = []
for l in range(nbLabels):
    index = round(len(color_panel)*l/nbLabels)
    if( index == len(color_panel) ):
        index -= 1
    colors_struct.append(color_panel[index])


legend_patch = []
for i in range(len(categories)):
    legend_patch.append(mpatches.Patch(color=colors_struct[i], label=categories[i]))

colors = seaborn.color_palette('Set2')

fs = 15
print("Plotting figures in "+path_to_plot+" :")
for i in range(nbFilter):
    print("\tFor filter "+FilterValues[i])
    for d in range(NbVp+NbVpComp):
        if( NbVpComp != 0 and d == NbVp+NbVpComp-1 ):
            print("\t\tFor dimensions in between "+str(LastLoopBeg+1)+" and "+str(LastLoopEnd))
        else:
            print("\t\tFor dimensions in between "+str(int(d*NbDimEachPlot)+1)+" and "+str(int((d+1)*NbDimEachPlot)))
        df = pd.DataFrame(dat_vp[d][i],columns=["artif","l","q"])
        ylims = [min(df["q"]),max(df["q"])]
        fig = plt.figure(figsize=(15,7))
        ax = fig.add_subplot(111)
        Ax = seaborn.violinplot(x="l",y="q",hue="artif",data=df, palette=['.4', '.7'], density_norm='count', split=True)
        count = 0
        for ind, violin in enumerate(Ax.findobj(PolyCollection)):
            for_color = int((ind-(int(ind/(2*nbLabels))*2*nbLabels))/2)
            rgb = to_rgb(colors_struct[for_color])
            violin.set_facecolor(rgb)
            if ind % 2 != 0:
                violin.set_alpha(.5)
        Ax.xaxis.set_major_locator(MultipleLocator(nbLabels))
        nbtab = int(40.*nbLabels/NbDimEachPlot)
        if( NbVpComp != 0 and d == NbVp+NbVpComp-1 ):
            nbtab = int(40.*nbLabels/(NbDim-(NbDimEachPlot*NbVp)))
        shiftval = d*NbDimEachPlot
        def custom_formatter(x, pos, nbtab=nbtab, shiftval=shiftval):
            tabs = ' '*nbtab
            return tabs+str(pos+shiftval)
        Ax.xaxis.set_major_formatter(FuncFormatter(custom_formatter))
        Ax.set_title("Filter : "+str(FilterValues[i]),fontweight="bold",fontsize=fs)
        Ax.set_ylim(ylims)
        plt.yticks(fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.ylabel(r'$\bar{q}_l$',fontsize=fs)
        plt.xlabel(r'l',fontsize=fs)
        plt.grid(axis='x')
        Ax.xaxis.label.set_visible(True)
        l0 = plt.legend(fontsize=fs,loc='upper right')
        l1 = plt.legend(handles=legend_patch,loc='upper left',ncol=3,fontsize=fs)
        plt.gca().add_artist(l1)
        plt.gca().add_artist(l0)
        if( NbVpComp != 0 ):
            if( d == NbVp+NbVpComp-1 ):
                plt.savefig(path_to_plot+'/'+str(FilterValues[i])+'Filter_'+str(LastLoopBeg+1)+'_'+str(LastLoopEnd)+'_Violin.png',dpi=500)
            else:
                plt.savefig(path_to_plot+'/'+str(FilterValues[i])+'Filter_'+str(int(d*NbDimEachPlot)+1)+'_'+str(int((d+1)*NbDimEachPlot))+'_Violin.png',dpi=500)
        else:
            plt.savefig(path_to_plot+'/'+str(FilterValues[i])+'Filter_'+str(int(d*NbDimEachPlot)+1)+'_'+str(int((d+1)*NbDimEachPlot))+'_Violin.png',dpi=500)
        print("\t\tDone")

