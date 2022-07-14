import seaborn as sns
import pandas as pd
from scipy import stats
import random
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_table('TS_runA_10_4.dat', names = ["Model No.","Stable States", "ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"], usecols=["ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"])
dfz= stats.zscore(df)
row = dfz.shape[0]
column = dfz.shape[1]

list = range (0, row)
pick = 25000
indexesrand= random.sample(list, pick)
print(indexesrand)
dflist= dfz.tolist()
arr = [[0]*column]*pick
i=0
for x in indexesrand:
    arr[i]= dflist[x]

    i=i+1
heatdat=pd.DataFrame(arr, columns = ["ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"])
heatmap=sns.clustermap(heatdat, metric= "euclidean", cmap="seismic",robust=True,figsize=(6, 7), dendrogram_ratio=(.1, .2), yticklabels=False, cbar_pos = (0.02, 0.85, 0.05, 0.18), font='Arial', fontsize = 12)
heatmap.savefig("Heatmap of Run A")