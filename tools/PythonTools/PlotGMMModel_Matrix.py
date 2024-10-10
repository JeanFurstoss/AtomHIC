#!/usr/bin/python3
import argparse
import numpy as np
import glob
import os
import os.path
import sys
from sklearn.mixture import GaussianMixture as GM
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import re

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

maxNbDat = 500

# Read the fitted GMM
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


#Creating and load the GM model
GMMModels = []
for i in range(nbFilter):
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
                dictLab[i][j].append(b)

def plot_ellipse(semimaj=1,semimin=1,phi=0,x_cent=0,y_cent=0,theta_num=1e3,ax=None,plot_kwargs=None,\
                    fill=False,fill_kwargs=None,data_out=False,cov=None,mass_level=0.68):
    # Get Ellipse Properties from cov matrix
    if cov is not None:
        eig_vec,eig_val,u = np.linalg.svd(cov)
        # Make sure 0th eigenvector has positive x-coordinate
        if eig_vec[0][0] < 0:
            eig_vec[0] *= -1
        semimaj = np.sqrt(eig_val[0])
        semimin = np.sqrt(eig_val[1])
        if mass_level is None:
            multiplier = np.sqrt(2.279)
        else:
            distances = np.linspace(0,20,20001)
            chi2_cdf = chi2.cdf(distances,df=2)
            multiplier = np.sqrt(distances[np.where(np.abs(chi2_cdf-mass_level)==np.abs(chi2_cdf-mass_level).min())[0][0]])
        semimaj *= multiplier
        semimin *= multiplier
        phi = np.arccos(np.dot(eig_vec[0],np.array([1,0])))
        if eig_vec[0][1] < 0 and phi > 0:
            phi *= -1

    # Generate data for ellipse structure
    theta = np.linspace(0,2*np.pi,int(theta_num))
    r = 1 / np.sqrt((np.cos(theta))**2 + (np.sin(theta))**2)
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    data = np.array([x,y])
    S = np.array([[semimaj,0],[0,semimin]])
    R = np.array([[np.cos(phi),-np.sin(phi)],[np.sin(phi),np.cos(phi)]])
    T = np.dot(R,S)
    data = np.dot(T,data)
    data[0] += x_cent
    data[1] += y_cent

    # Output data?
    if data_out == True:
        return data

    # Plot!
    return_fig = False
    if ax is None:
        return_fig = True
        fig,ax = plt.subplots()

    if plot_kwargs is None:
        ax.plot(data[0],data[1],color='b',linestyle='-')
    else:
        ax.plot(data[0],data[1],**plot_kwargs)

    if fill == True:
        ax.fill(data[0],data[1],**fill_kwargs)

    if return_fig == True:
        return fig

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


dataPoints = []
for f in range(nbFilter):
    dataPoints.append([])

print("Reading data")
for base, dirs, files in os.walk(dump_database_path):
    if( len(dirs) != 0 and len(dirs) != len(categories) ):
        print(dirs)
        print(categories)
        print("Inconsistant label size between GMM model and dump database")
    for directories in range(len(dirs)):
        for f in range(nbFilter):
            dataPoints[f].append([])
            dataPoints[f][directories] = np.empty((0, 20))
    for directories in range(len(dirs)):
        index_in_cat = -1
        for ind in range(len(categories)):
            if( categories[ind] == dirs[directories] ):
                index_in_cat = ind
                break
        list_file = sorted(glob.glob(base+dirs[directories]+"/*.cfg"))
        if len(list_file) == 0:
            continue
        
        # read first file to estimate number of point and reduce if necessary
        MeanNbDes = np.zeros(len(FilterValues))
        for input_file in list_file:
            Descriptors = readQ(input_file, FilterValues, NbDim)
            for i in range(len(FilterValues)):
                MeanNbDes[i] += len(Descriptors[i])
        
        MinNbDes = np.min(MeanNbDes)
        MinNbDes /= len(list_file)
        if( len(list_file)*MinNbDes > maxNbDat ):
            step_to_rem = round(maxNbDat/MinNbDes)
            list_file_tmp = []
            for f in range(step_to_rem):
                list_file_tmp.append(list_file[round((len(list_file)-1)/step_to_rem)])
            list_file = []
            for f in range(step_to_rem):
                list_file.append(list_file_tmp[f])
        for input_file in list_file:
            Descriptors = readQ(input_file, FilterValues, NbDim)
            for f in range(nbFilter):
                dataPoints[f][index_in_cat] = np.append(dataPoints[f][index_in_cat], Descriptors[f], axis=0)

color_panel = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'] # from colorbrewer2.org
colors_struct = []
for l in range(nbLabels):
    index = round(len(color_panel)*l/nbLabels)
    if( index == len(color_panel) ):
        index -= 1
    colors_struct.append(color_panel[index])

print("Plotting figures in "+path_to_plot+" :")
covs = np.empty((2,2))
means = np.empty((2))
for s in range(len(dataPoints)):
    print("\tFor filter "+FilterValues[s])
    fig, axs = plt.subplots(NbDim, NbDim, figsize=(90, 90))
    for i in range(NbDim):
        for b in range(NbDim):
            if b < i:
                continue
            for key in range(len(dataPoints[s])):
                plot_kwargs = {'color':colors_struct[key]}
                fill_kwargs = {'color':colors_struct[key],'alpha':0.3}
                axs[i, b].scatter(dataPoints[s][key][:,i], dataPoints[s][key][:,b], label=categories[key], c=colors_struct[key], linewidth=0.5, alpha=0.5)
                for dic in range(len(dictLab[s][key])):
                    for x_p in range(2):
                        means[x_p] = GMMModels[s].means_[dictLab[s][key][dic]][-(x_p-1)*i+x_p*b]
                        for y_p in range(2):
                            covs[x_p][y_p] = GMMModels[s].covariances_[dictLab[s][key][dic]][-(x_p-1)*i+x_p*b][-(y_p-1)*i+y_p*b]
                    U, s_ell, Vt = np.linalg.svd(covs)
                    angle = np.degrees(np.arctan2(U[1, 0], U[0, 0]))
                    width, height = 2 * np.sqrt(s_ell)
                    plot_ellipse(semimaj=width/2.,semimin=height/2.,phi=angle*np.pi/180,x_cent=means[0],y_cent=means[1],ax=axs[i, b],plot_kwargs=plot_kwargs,fill=True,fill_kwargs=fill_kwargs)
            if( i == 0 and b == 0 ):
                axs[i, b].legend()
            if( i == 0 and b == 1 ):
                axs[i, b].set_title("Filter : "+str(FilterValues[s]), fontsize=75)
            if i+1 < 21:
                axs[i, b].set_xlabel('q'+str(i+1))
            else:
                print("case")
                axs[i, b].set_xlabel('q'+str(i+1-11)+' mono')
            if b+1 < 21:
                axs[i, b].set_ylabel('q'+str(b+1))
            else:
                print("case")
                axs[i, b].set_ylabel('q'+str(b+2-11)+' mono')
    fig.tight_layout()
    fig.savefig(path_to_plot+'/'+str(FilterValues[s])+'Filter_Matrix.png')
    print("\tDone")
