import pandas as pd
import numpy as np
import scipy.stats as st
import random
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

colours = ['purple', 'green', 'orange', 'brown']
n_clusters = 2

def finding_threshold(score_data):
    overall_range=[-20,20]
    accuracy = 1000
    density = st.gaussian_kde(score_data)
    minima_range= [-2,2]
    data = []

    for i in range(overall_range[0]*accuracy,(overall_range[1]*accuracy)+1):
        data.append(density(i/accuracy)[0])       
    lb_index = int(abs((overall_range[0]*accuracy)-(minima_range[0]*accuracy)))
    len_minima_range = int((minima_range[1]-minima_range[0])*accuracy)
        
    minima_index = argrelextrema(np.array(data[lb_index:lb_index+len_minima_range]), np.less)

    threshold_value = ((minima_index[0]/accuracy)+minima_range[0])[0]
    return threshold_value

def emtscoring (gene_data):
    print(type(gene_data))
    zeb1 = np.array(gene_data['ZEB1'])
    miR200= np.array(gene_data['miR-200'])
    emtscore= st.zscore(np.subtract(zeb1, miR200))
    return emtscore

def stemnessscoring (gene_data):
    lin28 = np.array(gene_data['LIN28'])
    let7= np.array(gene_data['let-7'])
    stemnessscore= st.zscore(np.subtract(lin28, let7))
    return stemnessscore

def resistancescoring (gene_data):
    arv7 = np.array(gene_data['AR-v7'])
    ar= np.array(gene_data['AR'])
    resistancescore= st.zscore(np.add(ar, arv7))
    return resistancescore

def clustering (n_clusters, gene_data):
    kmeans = KMeans(n_clusters=n_clusters)
    clusterer=kmeans.fit_predict(gene_data)
    labels = kmeans.labels_
    print(labels)
    indices = np.array([np.where(labels==i) for i in range (0, n_clusters)], dtype = object)
    return indices

df = pd.read_table('TS_runA_10_4.dat', names = ["Model No.","Stable States", "ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"], usecols=["ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"])
dfz = pd.DataFrame(st.zscore(df), columns=["ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"])

emt = emtscoring(dfz)
resistance = resistancescoring(dfz)
stemness = stemnessscoring(dfz)

emtthres = 0.3340000000000001
resthres = -0.4870000000000001
stemthres = -0.548
indices = clustering(n_clusters, dfz)
quad=[0]*4
for i in range(0,dfz.shape[0]):
    if (emt[i]>emtthres and resistance[i]>resthres):
        quad[0]=quad[0]+1
    if (emt[i]<emtthres and resistance[i]>resthres):
        quad[1]=quad[1]+1
    if (emt[i]<emtthres and resistance[i]<resthres):
        quad[2]=quad[2]+1
    if (emt[i]>emtthres and resistance[i]<resthres):
        quad[3]=quad[3]+1
    
print("QUADRANTS: ", quad)
'''
plt.title('Stemness vs. EMT Scores')
plt.xlabel('EMT Scores')
plt.axvline(x = emtthres, color = 'red')
plt.ylabel('Stemness Scores')
plt.axhline(y = stemthres, color = 'black')
for i in range(0, n_clusters):
    index = list(indices[i])
    plt.scatter(emt[index],stemness[index], s = 0.2, color = colours[i] )
plt.savefig('Clust3MappedStemnessvsEMT.png')
plt.close()'''

plt.title('Resistance vs. EM Scores',  font = 'Arial', fontsize = 19)
plt.xlabel('EM Score', font = 'Arial', fontsize = 16)
plt.axvline(x = emtthres, color = 'red')
plt.ylabel('Resistance Score',  font = 'Arial', fontsize = 16)
plt.xticks(font='Arial', fontsize=12)
plt.yticks(font='Arial', fontsize=12)
plt.axhline(y = resthres, color = 'black')
for i in range(0, n_clusters):
    index = list(indices[i])
    plt.scatter(emt[index],resistance[index], s = 0.2, color = colours[i] )
plt.savefig('Clust2MappedResistancevsEM.png')
plt.close()

'''
plt.title('Stemness vs. Resistance Scores')
plt.xlabel('Resistance Scores')
plt.axvline(x = resthres, color = 'red')
plt.ylabel('Stemness Scores')
plt.axhline(y = stemthres, color = 'black')
for i in range(0, n_clusters):
    index = list(indices[i])
    plt.scatter(resistance[index],stemness[index], s = 0.2, color = colours[i] )
plt.savefig('Clust3MappedStemnessvsResistance.png')
plt.close()'''
