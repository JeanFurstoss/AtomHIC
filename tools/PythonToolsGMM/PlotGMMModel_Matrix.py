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
import matplotlib.patches as mpatches
import re

def main():
    # TODO when filter is none
    PathOfThisScript = os.path.dirname(__file__)+"/"
    dirModel = '../../data/MachineLearningModels/GaussianMixtureModel/'
    database_present = glob.glob(PathOfThisScript+dirModel+"/*/")
    help_print = ''
    for i in range(len(database_present)):
        help_print += database_present[i].replace(dirModel,'')+'\n'
    
    parser = argparse.ArgumentParser()
    parser.add_argument("AtomHICDatabaseName", help="Name of the database (already stored in AtomHIC /data/ or could be \"none\" if the ML model should not be plotted), available basis :\n"+help_print)
    parser.add_argument("MinClusterWeight", help="Minimum weight for a given cluster to be plotted")
    parser.add_argument("DumpsDatabasePath", help="Path of the directory containing the dumps used for the fitting of the GMM (could be \"none\" if the datapoints should not be plotted)")
    parser.add_argument("MaxNumberOfPlottedData", help="Maximum number of datapoints to be plotted")
    parser.add_argument("FilteringType", help="Filtering type (it should be one of the field in the .cfg dumps file)")
    parser.add_argument("DescriptorName", help="Name of the descriptor (it should be one of the field in the .cfg dumps file)")
    parser.add_argument("PathWhereSaveFigs", help="Path of the directory where to save figures")
    parser.parse_args()
    
    ath_database_name = sys.argv[1]
    minWeightForGMMPlot = float(sys.argv[2])
    dump_database_path = sys.argv[3]
    maxNbDat = int(sys.argv[4])
    FilterType = sys.argv[5]
    DescriptorName = sys.argv[6]
    path_to_plot = sys.argv[7]
    
    ptsize=1
    
    # Read the fitted GMM
    if( ath_database_name != "none" ):
        print("Reading "+ath_database_name+" AtomHIC database")
        
        os.makedirs(PathOfThisScript+dirModel, exist_ok=True)
        ModelName = os.path.join(ath_database_name, '')
        
        if not os.path.exists(PathOfThisScript+dirModel+ModelName):
            print("The database does not exists in /data/MachineLearningModels/GaussianMixtureModel/")
            print("Possible databases :")
            print(help_print)
        
        Labels = glob.glob(PathOfThisScript+dirModel+ModelName+"/*.ath")
        LabelsName = []
        nbLabels = len(Labels)
        for i in range(nbLabels):
            tmp = Labels[i].replace(PathOfThisScript+dirModel+ModelName,'')
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
    
   
    
    if( ath_database_name != "none" ):
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
    
   
    if( dump_database_path != "none" and ath_database_name == "none" ):
        categories = []
        nbFilter = 0
        NbDim = -1
        FilterValues = []
        for base, dirs, files in os.walk(dump_database_path):
            for directories in range(len(dirs)):
                categories.append(dirs[directories])
                list_file = sorted(glob.glob(base+dirs[directories]+"/*.cfg"))
                for input_file in list_file:
                    NewNbFilter, FilterValue_temp, NewNbDim = GetFilterAndDim(input_file, FilterType, DescriptorName)
                    for i in range(NewNbFilter):
                        FilterValues.append(FilterValue_temp[i])
                    if( NbDim == -1 ):
                        NbDim = NewNbDim
                    elif( NbDim != NewNbDim ):
                        print("Warning, the number of dimension of the descriptor is not consistent between all file of database")
                    if( NewNbFilter > nbFilter ):
                        nbFilter = NewNbFilter
        categories = (list(categories))
        categories = sorted(categories)
        nbLabels = len(categories)
        print("The different labels are:")
        print(categories)
        print("Filter values :")
        FilterValues = np.unique(np.array(FilterValues))
        print(FilterValues)
        print("Descriptor dimension: "+str(NbDim))
    
    if( dump_database_path != "none" ):
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
                    Descriptors = readQ(input_file, FilterValues, NbDim, FilterType, DescriptorName)
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
                    Descriptors = readQ(input_file, FilterValues, NbDim, FilterType, DescriptorName)
                    for f in range(nbFilter):
                        dataPoints[f][index_in_cat] = np.append(dataPoints[f][index_in_cat], Descriptors[f], axis=0)
    
    color_panel = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928'] # from colorbrewer2.org
    colors_struct = []
    for l in range(nbLabels):
        index = round(len(color_panel)*l/nbLabels)
        if( index == len(color_panel) ):
            index -= 1
        colors_struct.append(color_panel[index])
    
    if( ath_database_name != "none" ):
        print("Only clusters with weight higher than "+str(minWeightForGMMPlot)+" will be plotted")
    
    print("Plotting figures in "+path_to_plot+" :")
    covs = np.empty((2,2))
    means = np.empty((2))
    for s in range(nbFilter):
        print("\tFor filter "+FilterValues[s])
        fig, axs = plt.subplots(NbDim, NbDim, figsize=(90, 90))
        for i in range(NbDim):
            for b in range(NbDim):
                if b < i:
                    continue
                for key in range(nbLabels):
                    plot_kwargs = {'color':colors_struct[key]}
                    fill_kwargs = {'color':colors_struct[key],'alpha':0.3}
                    if( dump_database_path != "none" ):
                        axs[i, b].scatter(dataPoints[s][key][:,i], dataPoints[s][key][:,b], label=categories[key], c=colors_struct[key], linewidth=ptsize, alpha=0.5)
    
                    if( ath_database_name != "none" ):
                        for dic in range(len(dictLab[s][key])):
                            if( GMMModels[s].weights_[dictLab[s][key][dic]] > minWeightForGMMPlot ):
                                for x_p in range(2):
                                    means[x_p] = GMMModels[s].means_[dictLab[s][key][dic]][-(x_p-1)*i+x_p*b]
                                    for y_p in range(2):
                                        covs[x_p][y_p] = GMMModels[s].covariances_[dictLab[s][key][dic]][-(x_p-1)*i+x_p*b][-(y_p-1)*i+y_p*b]
                                plot_ellipse_bis(means, covs, axs[i,b], plot_kwargs, fill_kwargs)
                if( i == 0 and b == 0 ):
                    if( dump_database_path == "none" ):
                        patches = []
                        for key in range(nbLabels):
                            patches.append(mpatches.Patch(color=colors_struct[key], label=categories[key]))
                        axs[i,b].legend(handles=patches)
                    else:
                        axs[i, b].legend()
                if( i == 0 and b == 1 ):
                    axs[i, b].set_title("Filter : "+str(FilterValues[s]), fontsize=75)
                axs[i, b].set_xlabel('q'+str(i+1))
                axs[i, b].set_ylabel('q'+str(b+1))
        fig.tight_layout()
        #fig.savefig(path_to_plot+'/'+str(FilterValues[s])+'Filter_Matrix.pdf')
        fig.savefig(path_to_plot+'/'+str(FilterValues[s])+'Filter_Matrix.png',dpi=100)
        print("\tDone")

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
 
def plot_ellipse_bis(mean, cov, ax, plot_kwargs, fill_kwargs):
    eigval, eigvec = np.linalg.eig(cov)
    faceigval = 1.
    theta = np.linspace(0,2*np.pi,1000)
    if eigval[0] > 0 and eigval[1] > 0:
        eigval *= faceigval
        ellipsis = ( np.sqrt(eigval[None,:]) * eigvec ) @ [ np.sin(theta), np.cos(theta) ]
        ellipsis[0] += mean[0]
        ellipsis[1] += mean[1]
        ax.plot(ellipsis[0,:],ellipsis[1,:],**plot_kwargs)
        ax.fill(ellipsis[0,:],ellipsis[1,:],**fill_kwargs)

def readQ(input_file, FilterValue, nbDim, FilterType, DescriptorName):
    startCol = 0
    PosFilter = 0
    filter_val = []
    skip_head = 0
    count_line = 0
    clean_vals = []
    with open(input_file,'r') as file:
        found = False
        for line in file:
            count_line += 1
            line = line.replace('\t',' ')
            line = line.strip()
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

def GetFilterAndDim(input_file, FilterType, DescriptorName):
    startCol = 0
    PosFilter = 0
    filter_val = []
    skip_head = 0
    count_line = 0
    clean_vals = []
    for_nb_dim = []
    with open(input_file,'r') as file:
        found = False
        for line in file:
            count_line += 1
            line = line.replace('\t',' ')
            line = line.strip()
            if 'ITEM: ATOMS' in line:
                colvals = line.split(' ')
                for i in range(len(colvals)):
                    s = ''.join(x for x in colvals[i] if x.isdigit())
                    for_des = colvals[i].replace("["+s+"]",'')
                    if for_des == DescriptorName:
                        for_nb_dim.append(int(s))
                    if colvals[i] == FilterType :
                        PosFilter = i-2
                found = True
            elif found:
                filter_val.append((line.split(' ')[PosFilter]).strip())
    filter_val = np.unique(np.array(filter_val))
    for_nb_dim = np.array(for_nb_dim)
    return len(filter_val), filter_val, np.max(for_nb_dim)
 
if __name__ == '__main__':
    main()
