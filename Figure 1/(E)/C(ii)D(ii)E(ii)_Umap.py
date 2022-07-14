import os
import pandas as pd
import matplotlib.pyplot as plt
import umap
from sklearn.preprocessing import StandardScaler
import sklearn
import umap
import scipy
from scipy import stats
import numpy as np

df = pd.read_table ('TS_RunA_10_4.dat', names = ["Model No.","Stable States", "ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"], usecols=["ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"])
state_dataframe = pd.DataFrame(stats.zscore(df), columns = ["ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"])

 

sub_dataframe = state_dataframe.sample(n=30000)
emt_df_1 = sub_dataframe['ZEB1']-sub_dataframe['miR-200']
resistance_df_1 =sub_dataframe['AR']+sub_dataframe['AR-v7']
emt_df = pd.DataFrame(stats.zscore(emt_df_1))
resistance_df =pd.DataFrame(stats.zscore(resistance_df_1))
emt_mat = emt_df.to_numpy()
resistance_mat = resistance_df.to_numpy()

 

print(len(resistance_mat))
print(len(emt_mat))
print(emt_mat[0])

 

quad_mat = np.zeros(len(emt_mat))

 


for i in range (len(resistance_mat)):
    if ((emt_mat[i] < 0.3340000000000001) and ( resistance_mat[i]< -0.4870000000000001)):
      quad_mat[i] = 3
    if ((emt_mat[i] > 0.3340000000000001) and ( resistance_mat[i]> -0.4870000000000001)):
       quad_mat[i] = 1
    if ((emt_mat[i] < 0.3340000000000001) and ( resistance_mat[i]> -0.4870000000000001)):
      quad_mat[i] = 2
    if ((emt_mat[i] > 0.3340000000000001) and ( resistance_mat[i]< -0.4870000000000001)):
      quad_mat[i] = 4
print(quad_mat)

 

quad_df = pd.DataFrame(quad_mat)

 

def Umap_analysis(data,n_neighbors=100, min_dist=0.8, n_components=2, metric='euclidean'):
    fit = umap.UMAP(random_state = 42)
    umap_data = fit.fit_transform(data)

    return umap_data

#This function plots the output from the Umap analysis. 
def Umap_scatter(state_dataframe):
    sub_dataframe = state_dataframe.sample(n=30000)
    embedding = Umap_analysis(sub_dataframe,n_neighbors=100)
    print(embedding)
    Umap = pd.DataFrame()
    emt_df_1 = sub_dataframe['ZEB1']-sub_dataframe['miR-200']
    resistance_df_1 =sub_dataframe['AR']+sub_dataframe['AR-v7']
    emt_df = pd.DataFrame(stats.zscore(emt_df_1))
    resistance_df =pd.DataFrame(stats.zscore(resistance_df_1))
    emt_mat = emt_df.to_numpy()
    resistance_mat = resistance_df.to_numpy()


    quad_mat = np.zeros(len(emt_mat))
    quad_1 = []
    quad_2 = []
    quad_3 = []
    quad_4 = []
    for i in range (len(resistance_mat)):
      if ((emt_mat[i] < 0.3340000000000001) and ( resistance_mat[i]< -0.4870000000000001)):
        quad_mat[i] = 3
        quad_3.append([i])
      if ((emt_mat[i] > 0.3340000000000001) and ( resistance_mat[i]> -0.4870000000000001)):
        quad_mat[i] = 1
        quad_1.append([i])
      if ((emt_mat[i] < 0.3340000000000001) and ( resistance_mat[i]> -0.4870000000000001)):
        quad_mat[i] = 2
        quad_2.append([i])
      if ((emt_mat[i] > 0.3340000000000001) and ( resistance_mat[i]< -0.4870000000000001)):
        quad_mat[i] = 4
        quad_4.append([i])
    #quad_1= np.where(quad_mat==1)
    #quad_2= np.where(quad_mat==2)
    #quad_3= np.where(quad_mat==3)
    #quad_4= np.where(quad_mat==4)
    print(quad_mat)
    print('***********QUADRANTS***********')
    print(quad_1, quad_2, quad_3, quad_4)

    plt.scatter(embedding[quad_1,0], embedding[quad_1,1], s= 1, color = 'blue', label = "Quad1")
    plt.scatter(embedding[quad_2,0], embedding[quad_2,1], s= 1, color = 'green', label = 'Quad 2')
    plt.scatter(embedding[quad_3,0], embedding[quad_3,1], s= 1, color = 'red', label = 'Quad 3')
    plt.scatter(embedding[quad_4,0], embedding[quad_4,1], s= 1, color = 'brown', label = 'Quad 4')
    plt.title('UMAP based on EM and Resistance Scores')
    plt.xlabel('UMAP 1')
    plt.ylabel('UMAP 2')
    plt.legend()
    plt.savefig('ScatterUMAPWORed')
    plt.show()

 

    quad_df = pd.DataFrame(quad_mat)

 
    Umap['UMAP_1'] = embedding[:,0]
    #coord_1_df = embedding[:,0]
    #coord_2_df = embedding[:,1]
    Umap['UMAP_2'] = embedding[:,1]
    Umap.plot.scatter(x = 'UMAP_1',y='UMAP_2',s=1,c=sub_dataframe['ZEB1']-sub_dataframe['miR-200'],cmap=plt.cm.nipy_spectral)
    plt.title('EMT score UMAP',fontsize=14)
    print('hiiii')
    plt.xlabel('UMAP_1',fontsize=14)
    plt.ylabel('UMAP_2',fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig('EM Score UMAP')
    plt.show()
    plt.close()

 


    Umap.plot.scatter(x='UMAP_1',y='UMAP_2',s=1,c=sub_dataframe['AR']+sub_dataframe['AR-v7'],cmap=plt.cm.nipy_spectral)
    plt.title('Resistance score UMAP',fontsize=14)
    plt.xlabel('UMAP_1',fontsize=14)
    plt.ylabel('UMAP_2',fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig('Resistance score UMAP')
    plt.close()

 
'''
    Umap.plot.scatter(x='UMAP_1',y='UMAP_2',s=1,c=quad_df,cmap=plt.cm.nipy_spectral)
    plt.title('UMAP',fontsize=14)
    plt.xlabel('UMAP_1',fontsize=14)
    plt.ylabel('UMAP_2',fontsize=14)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig('UMAP_1')
    plt.close()'''
    
     

Umap_scatter(state_dataframe)

  

print('hi')