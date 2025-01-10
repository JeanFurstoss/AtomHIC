#!/usr/bin/python3
import argparse
import numpy as np
from numpy import linalg as la
import glob
import os
import os.path
import sys
from sklearn.mixture import GaussianMixture as GM
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import re
# TODO when filter is none
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

maxNbDat = 1000000
minWeightForGMMPlot = 0.001
ptsize=1

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
                dictLab[i][j].append(b)

#def plot_MLC(means

def plot_ellipse_bis(mean, cov, ax, plot_kwargs, fill_kwargs):
    eigval, eigvec = np.linalg.eig(cov)
    faceigval = 1.
    #print("NEW")
    #print(cov)
    #print(eigval)
    #print(eigvec)
    theta = np.linspace(0,2*np.pi,1000)
    if eigval[0] > 0 and eigval[1] > 0:
        eigval *= faceigval
        ellipsis = ( np.sqrt(eigval[None,:]) * eigvec ) @ [ np.sin(theta), np.cos(theta) ]
        ellipsis[0] += mean[0]
        ellipsis[1] += mean[1]
        ax.plot(ellipsis[0,:],ellipsis[1,:],**plot_kwargs)
        ax.fill(ellipsis[0,:],ellipsis[1,:],**fill_kwargs)

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
    skip_head = 0
    count_line = 0
    with open(input_file,'r') as file:
        found = False
        for line in file:
            count_line += 1
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
                skip_head = count_line
            elif found:
                filter_val.append((line.split(' ')[PosFilter]).strip())

    pos_tmp=np.genfromtxt(input_file,skip_header=skip_head)
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
            dataPoints[f][directories] = np.empty((0, NbDim))
    print("test")
    print(dataPoints[0][0])
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

### TEST
#q_art = []
#lab = []
#q_GMM = []
#nsamp = 1000
#for l in range(nbFilter):
#    q_art_temp, lab_temp = GMMModels[l].sample(nsamp)
#    q_art.append(q_art_temp)
#    lab.append(lab_temp)
#    q_GMM.append([])
#    for i in range(len(categories)):
#        q_GMM[-1].append([])
#
#for i in range(nsamp):
#    found = []
#    cat = []
#    for l in range(nbFilter):
#        found.append(False)
#        cat.append(-1)
#    for b in range(len(categories)):
#        for l in range(nbFilter):
#            if not found[l]:
#                for s in range(len(dictLab[l][b])):
#                    if( lab[l][i] == dictLab[l][b][s] ):
#                        cat[l] = b
#                        found[l] = True
#                        break
#        all_found = True
#        for l in range(nbFilter):
#            all_found *= found[l]
#        if all_found:
#            break
#    for l in range(nbFilter):
#        if found[l]:
#            q_GMM[l][cat[l]].append(q_art[l][i])
### END TEST


print("Only clusters with weight higher than "+str(minWeightForGMMPlot)+" will be plotted")
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
                axs[i, b].scatter(dataPoints[s][key][:,i], dataPoints[s][key][:,b], label=categories[key], c=colors_struct[key], linewidth=ptsize, alpha=0.5)
                #to_plot_i = []
                #to_plot_b = []
                #for samp in range(len(q_GMM[s][key])):
                #    to_plot_i.append(q_GMM[s][key][samp][i])
                #    to_plot_b.append(q_GMM[s][key][samp][b])
                #axs[i, b].scatter(to_plot_i, to_plot_b, label=categories[key], c=colors_struct[key], linewidth=ptsize, alpha=0.5)
                for dic in range(len(dictLab[s][key])):
                    if( GMMModels[s].weights_[dictLab[s][key][dic]] > minWeightForGMMPlot ):
                        for x_p in range(2):
                            means[x_p] = GMMModels[s].means_[dictLab[s][key][dic]][-(x_p-1)*i+x_p*b]
                            for y_p in range(2):
                                covs[x_p][y_p] = GMMModels[s].covariances_[dictLab[s][key][dic]][-(x_p-1)*i+x_p*b][-(y_p-1)*i+y_p*b]
                        plot_ellipse_bis(means, covs, axs[i,b], plot_kwargs, fill_kwargs)
                        #U, s_ell, Vt = np.linalg.svd(covs)
                        #angle = np.degrees(np.arctan2(U[1, 0], U[0, 0]))
                        #width, height = 2 * np.sqrt(s_ell)
                        #plot_ellipse(semimaj=width/2.,semimin=height/2.,phi=angle*np.pi/180,x_cent=means[0],y_cent=means[1],ax=axs[i, b],plot_kwargs=plot_kwargs,fill=True,fill_kwargs=fill_kwargs)
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
    fig.savefig(path_to_plot+'/'+str(FilterValues[s])+'Filter_Matrix.pdf')
    print("\tDone")
