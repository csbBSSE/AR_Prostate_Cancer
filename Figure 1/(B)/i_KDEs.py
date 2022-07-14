''' Plot the kde plot with histograms for each node:
Take the data of each node
Perform z-score analysis on the data
Plot the histogram and kde
Save the image
'''


import pandas as pd 
import seaborn as sns; sns.set_theme(color_codes=True)
import scipy
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt



fig, ax = plt.subplots(nrows=3, ncols=3, sharex=True, sharey=True, figsize=(9,9))

fig.text(0.5, 0.02, 'Z-normalised score', ha='center', fontsize = 22)
fig.text(0.03, 0.5, 'Kernel Density Estimate', va='center', rotation='vertical', fontsize = 22)

genes = ["ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"]

df = pd.read_table ('TS_RunA_10_4.dat', names = ["Model No.","Stable States", "ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"], usecols=genes)
dfz= pd.DataFrame(stats.zscore(df))

i=0
arr =np.transpose(np.array([[0.0,0.751,0.989,0.989,0.0,0.0,0.0,0.0,0.931],[0.581,0.435,0.430,0.441,0.642,0.630,0.669,0.578,0.425]]))
print(arr)
bimod = pd.DataFrame(data=arr, columns= ['HDT','BC'], index=genes)
print(bimod.loc['ZEB1','HDT'])

dfz = np.array(dfz)

for i in range(0, 3):
    for j in range(0,3):
        g = genes[i*3+j]
        ax[i,j].hist(x= dfz[:,(i*3+j)], density = True, bins = 24, color = 'lightsteelblue', label = g)
        sns.kdeplot(data = dfz[:,(i*3+j)], ax = ax[i][j])
        xt = ax[i,j].get_xticks()
        ax[i,j].set_xticklabels(xt, fontsize = 16)
        yt = ax[i,j].get_yticks()
        ax[i,j].set_yticklabels(xt, fontsize = 16)
        ax[i,j].set_ylabel('')
        ax[i,j].set_title(g, fontsize = 16)
        textstr="BC: "+str(bimod.loc[g,'BC'])
        ax[i,j].text(-5 ,0.55, textstr, fontsize = 14)
        j=j+1
    i=i+1

plt.savefig('Common KDE Larger')
plt.show()
