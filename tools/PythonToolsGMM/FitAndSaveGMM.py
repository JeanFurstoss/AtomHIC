#!/usr/bin/python3

from joblib import Parallel, delayed
import argparse
import numpy as np
import glob
import os
import os.path
import sys
from sklearn.mixture import GaussianMixture as GM
from sklearn.model_selection import GridSearchCV
import re

dirModelToPrint = '../../data/MachineLearningModels/GaussianMixtureModel/'

parser = argparse.ArgumentParser()
parser.add_argument("DumpsDatabasePath", help="Path of the directory containing the dumps used for the fitting of the GMM, this directory should contain subdirectories having the name of the labels and a DescriptorProperties.ath file (an example file could be found in /data/ExampleFiles/) describing the properties of the descriptors used and filtering type")
parser.add_argument("DescriptorName", help="Name of the descriptor (it should be one of the field in the .cfg dumps file)")
parser.add_argument("nclust_min", help="low bound of the range for searching the optimal GMM number of cluster")
parser.add_argument("nclust_max", help="upper bound of the range for searching the optimal GMM number of cluster")
parser.add_argument("AtomHICDatabaseName", help="Name of the database (will be stored in AtomHIC /data/MachineLearningModels/GaussianMixtureModel/)")
parser.parse_args()

dump_database_path = sys.argv[1]
DescriptorName = sys.argv[2]
n_clust_min = int(sys.argv[3])
n_clust_max = int(sys.argv[4])
ath_database_name = sys.argv[5]

# Create the atomhic database directory
PathOfThisScript = os.path.dirname(__file__)
database_present = glob.glob(PathOfThisScript+dirModelToPrint+"/*/")
for i in range(len(database_present)):
    database_present[i] = database_present[i].replace(PathOfThisScript+dirModelToPrint,'')
    database_present[i] = database_present[i].replace('/','')

# change of name if already exist
BaseNameOk = False
ext=0
while not BaseNameOk:
    BaseNameOk = True
    for n in range(len(database_present)):
        if( ath_database_name == database_present[n] ):
            BaseNameOk = False
            ext += 1
            if ext != 1:
                if ext < 11:
                    ath_database_name = ath_database_name[:-2]
                if ext > 10 and ext < 101:
                    ath_database_name = ath_database_name[:-3]
                elif ext > 100 and ext < 1001:
                    ath_database_name = ath_database_name[:-4]
                elif ext > 1000 and ext < 10001:
                    ath_database_name = ath_database_name[:-5]
            ath_database_name += '_'+str(int(ext))
            break

os.mkdir(PathOfThisScript+dirModelToPrint+ath_database_name)

print("Reading "+dump_database_path+" database")

DesProp = glob.glob(dump_database_path+"/*.ath")
IsDesProp = False

FilterType = "none"
NbDim = 0
DescriptorProperties = []

for i in range(len(DesProp)):
    if( DesProp[i] == dump_database_path+"DescriptorProperties.ath" ):
        print("Reading DesscriptorProperties.ath file")
        with open(DesProp[i],'r') as file: 
            for line in file:
                DescriptorProperties.append(line)
                if 'FILTER_TYPE' in line:
                    FilterType = (line.split(' ')[1]).strip()
                if 'NUMBER_OF_DIMENSION' in line:
                    NbDim = int((line.split(' ')[1]).strip())
        IsDesProp = True
        break

def getNbDimFromDump(input_file, DescriptorName=DescriptorName):
    line_header = 0
    dim = 0
    found = False
    DesFound = False
    with open(input_file,'r') as file:
        for line in file:
            if 'ITEM: ATOMS' in line:
                colvals = line.split(' ')
                for i in range(len(colvals)):
                    s = ''.join(x for x in colvals[i] if x.isdigit())
                    for_des = (colvals[i].replace("["+s+"]",'')).strip()
                    if for_des == DescriptorName:
                        dim += 1
                        DesFound = True
                found = True
                break
    if not found:
        print("The provided dump file does not contain \"ITEM: ATOMS\" tag of cfg file needed for identifying the properties of atoms")
    if not DesFound:
        print("The provided descriptor name does not correspond to an auxiliary property of the dump file")
    return dim 

def getFilterValueFromDump(input_file, FilterType=FilterType):
    PosFilter = 0
    filter_val = []
    FilterValue = []
    found = False
    FilFound = False
    with open(input_file,'r') as file:
        for line in file:
            if 'ITEM: ATOMS' in line:
                colvals = line.split(' ')
                for i in range(len(colvals)):
                    if colvals[i] == FilterType :
                        PosFilter = i-2
                        FilFound = True
                found = True
            elif found:
                filter_val.append((line.split(' ')[PosFilter]).strip())
    if not found:
        print("The provided dump file does not contain \"ITEM: ATOMS\" tag of cfg file needed for identifying the properties of atoms")
    if not FilFound:
        print("The provided filtering type does not correspond to an auxiliary property of the dump file")
    # make unique filter_val for having FilterValue array
    templist = set(filter_val)
    FilterValue = (list(templist))
    return FilterValue 

FilterValue = []
if not IsDesProp:
    print("DescriptorProperties.ath file not found in "+dump_database_path+", we consider that descriptors are not filtered and will print no properties of the descriptors in the GMM database")
    FilterValue.append("none")

if not IsDesProp or NbDim == 0:
    print("Searching number of dimension from dump database")
    # search dimension by reading the files
    first = True
    oldNbDim = 0
    for base, dirs, files in os.walk(dump_database_path):
        for directories in dirs:
            list_file = sorted(glob.glob(base+directories+"/*.cfg"))
            if len(list_file) == 0:
                continue
            for input_file in list_file:
                if first:
                    oldNbDim = getNbDimFromDump(input_file, DescriptorName=DescriptorName)
                    first = False
                else:
                    newNbDim = getNbDimFromDump(input_file, DescriptorName=DescriptorName)
                    if oldNbDim != newNbDim:
                        print("Warning, the descriptor dimension are inconsistant between the different dumps in the database (particularly in "+input_file+")")
    NbDim = oldNbDim
    print("From dump database, we infer a descriptor dimension of "+str(NbDim))
    if IsDesProp:
        for line in range(len(DescriptorProperties)):
            if 'NUMBER_OF_DIMENSION' in DescriptorProperties[line]:
                DescriptorProperties[line] = 'NUMBER_OF_DIMENSION '+str(NbDim)+"\n"
                break

if IsDesProp:
    # get the FilterValue array
    print("Collecting filter values from dump database")
    for base, dirs, files in os.walk(dump_database_path):
        for directories in dirs:
            list_file = sorted(glob.glob(base+directories+"/*.cfg"))
            if len(list_file) == 0:
                continue
            for input_file in list_file:
                filter_value_temp = getFilterValueFromDump(input_file, FilterType=FilterType)
                for i in range(len(filter_value_temp)):
                    already = False
                    for j in range(len(FilterValue)):
                        if( filter_value_temp[i] == FilterValue[j] ):
                            already = True
                            break
                    if not already:
                        FilterValue.append(filter_value_temp[i])
    print("From dump database, the filter can take the following values :")
    print(FilterValue)


nbFilter=len(FilterValue)

def readQ(input_file, FilterValue, nbDim, FilterType=FilterType, DescriptorName=DescriptorName):
    startCol = 0
    PosFilter = 0
    filter_val = []
    line_header = 0
    Descriptors = []
    count_line = 0
    found = False
    DesFound = False
    FilFound = False
    with open(input_file,'r') as file:
        for line in file:
            count_line += 1
            if 'ITEM: ATOMS' in line:
                colvals = line.split(' ')
                for i in range(len(colvals)):
                    s = ''.join(x for x in colvals[i] if x.isdigit())
                    for_des = colvals[i].replace("["+s+"]",'')
                    if for_des == DescriptorName and not DesFound:
                        startCol = i-2
                        DesFound = True
                    if colvals[i] == FilterType :
                        PosFilter = i-2
                        FilFound = True
                found = True
                line_header = count_line
            elif found:
                filter_val.append((line.split(' ')[PosFilter]).strip())
    if not found:
        print("The provided dump file does not contain \"ITEM: ATOMS\" tag of cfg file needed for identifying the properties of atoms")
    if not DesFound:
        print("The provided descriptor name does not correspond to an auxiliary property of the dump file")
    if not FilFound and FilterType != "none":
        print("The provided filtering type does not correspond to an auxiliary property of the dump file")
 
    pos_tmp=np.genfromtxt(input_file,skip_header=line_header)
    if FilterType != "none":
        for f in range(len(FilterValue)):
            Descriptors.append([])
            Descriptors[-1] = [pos_tmp[i,startCol:startCol+nbDim] for i in range(len(pos_tmp)) if filter_val[i] == FilterValue[f] ]
    else:
        Descriptors.append([])
        Descriptors[-1] = [pos_tmp[i,startCol:startCol+nbDim] for i in range(len(pos_tmp))]
    return Descriptors 


print('Number of atoms in each structure')
#Read structures from database
categories = []
Structs = []
X = []
for i in range(nbFilter):
    Structs.append([])
    X.append(np.empty((0,NbDim)))

for base, dirs, files in os.walk(dump_database_path):
    for directories in dirs:
        list_file = sorted(glob.glob(base+directories+"/*.cfg"))
        if len(list_file) == 0:
            continue
        for f in range(nbFilter):
            Structs[f].append(np.empty((0,NbDim)))
        for input_file in list_file:
            Descriptors = readQ(input_file,FilterValue,NbDim,FilterType,DescriptorName)
            for f in range(nbFilter):
                X[f] = np.append(X[f], Descriptors[f], axis=0)
                Structs[f][-1] = np.append(Structs[f][-1], Descriptors[f], axis=0)
        print(directories)
        for f in range(nbFilter):
            print(FilterValue[f]+': ', Structs[f][-1].shape[0])
        print('------------------------------------')
        categories.append(directories)


#Save training information to file
class Tee(object):
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush()
    def flush(self):
        for f in self.files:
            f.flush()


def gmm_bic_score(estimator, X):
    return -estimator.bic(X)

#Find out number of clusters
numClusters = []
param_grid = {"n_components": range(n_clust_min, n_clust_max),}
for f in range(nbFilter):
    grid_seach = GridSearchCV(GM(), param_grid=param_grid, scoring=gmm_bic_score)
    grid_seach.fit(X[f])
    print(f)
    print(grid_seach.best_estimator_)
    numClusters.append(grid_seach.best_params_['n_components'])

print('Number of clusters used in each model: ')
for f in range(nbFilter):
    print('Clusters '+FilterValue[f]+' --- ', numClusters[f])
print('--------------------------------')
for f in range(nbFilter):
    print('Datapoints '+FilterValue[f]+' --- ', X[f].shape)
print()
print('Categories: ', categories)

#Creating and training the GM model
model = []
clusterIDs = []
prevClust = []
fout = open(PathOfThisScript+dirModelToPrint+ath_database_name+'/labelling.out','w')
for f in range(nbFilter):
    model.append(GM(numClusters[f], n_init=100))
    model[f].fit(X[f])
    clusterIDs.append(['']*numClusters[f])
    prevClust.append(np.zeros(numClusters[f]))
    for i in range(len(Structs[f])):
        testProbs = model[f].predict_proba(Structs[f][i])
        aveProb = np.mean(testProbs, axis=0)
        print('-----------------------')
        print('Test structure: '+categories[i])
        print('Probabilities '+FilterValue[f]+': ', aveProb)
        fout.write('-----------------------'+'\n')
        fout.write('Test structure: '+categories[i]+'\n')
        fout.write('Probabilities '+FilterValue[f]+': '+str(aveProb)+'\n')
        for b in range(len(clusterIDs[f])):
            if aveProb[b] > prevClust[f][b]:
                clusterIDs[f][b] = categories[i]
                prevClust[f][b] = aveProb[b]

    print('Identified clusters '+FilterValue[f]+': ', clusterIDs[f])
    fout.write('Identified clusters '+FilterValue[f]+': '+str(clusterIDs[f])+'\n')
    for i in range(len(clusterIDs[f])):
        if clusterIDs[f][i] == '':
            print('Not all clusters have been recognized!')
            exit()
    print('All clusters have been properly recognized')
    fout.write('All clusters have been properly recognized\n')
fout.close()

print("Writting database in "+ath_database_name)
for i in range(len(categories)):
    fo = open(PathOfThisScript+dirModelToPrint+ath_database_name+'/'+categories[i]+'.ath','w')
    if not IsDesProp: 
        fo.write("NUMBER_OF_DIMENSION "+str(NbDim)+"\n")
        fo.write("FILTER_TYPE none\n")
    else:
        for l in range(len(DescriptorProperties)):
            fo.write(DescriptorProperties[l])
    for f in range(nbFilter):
        fo.write("FILTER_VALUE "+FilterValue[f]+"\n")
        nbclust = 0
        for j in range(len(clusterIDs[f])):
            if( clusterIDs[f][j] == categories[i] ):
                nbclust += 1
        fo.write("NUMBER_OF_CLUSTER "+str(nbclust)+"\n")
        for j in range(len(clusterIDs[f])):
            if( clusterIDs[f][j] == categories[i] ):
                fo.write("WEIGHT "+str(model[f].weights_[j])+"\n")
                fo.write("DETERMINANT_OF_COV_MATRIX "+str(np.linalg.det(model[f].covariances_[j]))+"\n")
                fo.write("ESPERANCE")
                for n in range(NbDim):
                    fo.write(" "+str(model[f].means_[j][n]))
                fo.write("\n")
                inv_cov = np.linalg.inv(model[f].covariances_[j])
                fo.write("INVERSE_OF_COVARIANCE_MATRIX\n")
                for n1 in range(NbDim):
                    for n2 in range(NbDim):
                        fo.write(str(inv_cov[n1][n2])+" ")
                    fo.write("\n")
    fo.close()

