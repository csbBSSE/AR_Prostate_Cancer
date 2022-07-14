from string import ascii_letters
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import statistics
import fontstyle


sns.set_theme(style="white")
''' df_1 is the dataframe, which has the columns in the order of "Model No.","Stable States", "PD-L1", "miR-200","SNAI1","let-7","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7","AR"], usecols=["PD-L1", "miR-200","SNAI1","let-7","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7","AR"'''
df_1 = pd.read_table ('TS_CircuitB_arranged.dat', names = ["Model No.","Stable States", "PD-L1", "miR-200","SNAI1","let-7","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7","AR"], usecols=["PD-L1", "miR-200","SNAI1","let-7","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7","AR"])

#print(df_1)
corr1 = df_1.corr()
#print(corr1)
corrmat1 = corr1.to_numpy()


mean_mat = np.zeros(shape=(10, 10))
stdev_mat = np.zeros(shape=(10, 10))
mean_mat= corrmat1

#print(mean_mat)

annot_mat = mean_mat
annot_matrix = np.round_(annot_mat, decimals = 2)
annot_lis = annot_matrix.tolist()



def corr_sig(df=None):
    p_matrix = np.zeros(shape=(df.shape[1],df.shape[1]))
    for col in df.columns:
        for col2 in df.drop(col,axis=1).columns:
            _ , p = stats.pearsonr(df[col],df[col2])
            p_matrix[df.columns.to_list().index(col),df.columns.to_list().index(col2)] = p
    return p_matrix



# Generate a mask for the upper triangle
mask = np.triu(np.ones_like(mean_mat, dtype=bool))

# Set up the matplotlib figure
f, ax = plt.subplots(figsize=(8,8))

# Generate a custom diverging colormap
cmap = sns.diverging_palette(230, 20, as_cmap=True)

# Draw the heatmap with the mask and correct aspect ratio
xticklabels = ["PD-L1", "miR-200","SNAI1","let-7","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7","AR"]
yticklabels = ["PD-L1", "miR-200","SNAI1","let-7","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7","AR"]
plt.yticks(fontsize=14)
sns.heatmap(mean_mat, mask=mask, cmap=cmap, vmax=.3, center=0,annot=True,yticklabels=yticklabels,
            xticklabels=xticklabels, square=True, linewidths=.5, cbar_kws={"shrink": .5})

plt.savefig('Correlation Map with Mean value of three runs_Network 2')
#plt.show()
#plt.close()

# Generate a mask for the upper triangle
mask = np.triu(np.ones_like(mean_mat, dtype=bool))

# Set up the matplotlib figure
f, ax = plt.subplots(figsize=(8,8))

# Generate a custom diverging colormap
cmap = sns.diverging_palette(230, 20, as_cmap=True)

# Draw the heatmap with the mask and correct aspect ratio
xticklabels = ["PD-L1", "miR-200","SNAI1","let-7","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7","AR"]
yticklabels = ["PD-L1", "miR-200","SNAI1","let-7","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7","AR"]
sns.heatmap(stdev_mat, mask=mask, cmap=cmap, vmax=.3, center=0,annot=True,yticklabels=yticklabels,
            xticklabels=xticklabels,square=True, linewidths=.5, cbar_kws={"shrink": .5})
plt.savefig('Standard Deviation Map 1_network2')
#plt.show()
#plt.close()

df_mean_mat = pd.DataFrame(mean_mat)
#print(df_mean_mat)

df_stdev_mat = pd.DataFrame(stdev_mat)

p_values = corr_sig(df_mean_mat)                     # get p-Value
#mask = np.invert(np.tril(p_values<0.05))    # mask - only get significant corr
f, ax = plt.subplots(figsize=(8,8))
cmap = sns.diverging_palette(230, 20, as_cmap=True)
xticklabels = ["PD-L1", "miR-200","SNAI1","let-7","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7","AR"]
yticklabels = ["PD-L1", "miR-200","SNAI1","let-7","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7","AR"]
sns.heatmap(p_values, mask=mask, cmap=cmap, vmax=.3, center=0,annot=True,yticklabels=yticklabels,
            xticklabels=xticklabels, square=True, linewidths=.5, cbar_kws={"shrink": .5})
plt.savefig('p-value of the three runs 1_Network2')
#plt.show()
#plt.close()

#print(type(p_values))



for i in range (0,10):
    for j in range (0,10):
        if p_values[i,j] > 0.05:
            var = annot_lis[i][j]
            annot_lis[i][j] = str(var) + '**'
        else:
            var = annot_lis[i][j]
            annot_lis[i][j] = str(var)
            

#print(annot_lis)
some = np.array(annot_lis)
#print(some)
#print(type(some))
def plot_cor_matrix(corr, mask=None):
    f, ax = plt.subplots(figsize=(11, 9))
    sns.heatmap(corr, ax=ax,
                mask=mask,
                # cosmetics
                annot=True, vmax=0.3, center=0,yticklabels=yticklabels,
                xticklabels=xticklabels,cmap=cmap, linewidths=.5, cbar_kws={"shrink": .5})


#corr = df.corr()                            # get correlation
#p_values = corr_sig(df)                     # get p-Value
# mask - only get significant corr
fig, ax = plt.subplots(figsize=(9,9))
sns.heatmap(df_mean_mat, mask=mask,annot=annot_lis,fmt = '', vmax=1,vmin=-1, center=0,yticklabels=yticklabels,
                xticklabels=xticklabels,cmap=cmap,square=True, linewidths=.5, cbar_kws={"shrink": .7},annot_kws={"size": 15})
#plot_cor_matrix(df_mean_mat,mask)

plt.yticks(rotation=45)
plt.xticks(rotation=45)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
#plt.title('Correlation Plot with significantly important data')
plt.savefig('Correlation Plot with significantly important data 1_Network2')
#plt.show()



#mask = np.invert(np.tril(p_values<0.05))    # mask - only get significant corr
f, ax = plt.subplots(figsize=(11, 9))
sns.heatmap(df_stdev_mat, ax=ax,mask=mask,annot=annot_lis,fmt = '', vmax=0.3, center=0,yticklabels=yticklabels,
                xticklabels=xticklabels,cmap=cmap, linewidths=.5, cbar_kws={"shrink": .5})
plt.yticks(rotation=45)
plt.xticks(rotation=45)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
#plt.title('Standard Deviation for the correlation plot ')
plt.savefig('Standard Deviation for the correlation plot_Network2 ')
#plt.show()



for i in range (0,10):
    for j in range (0,10):
        if p_values[i,j] > 0.05:
            var = annot_lis[i][j]
            annot_lis[i][j] = 'X'
        else:
            var = annot_lis[i][j]
            annot_lis[i][j] = ' '

xticklabels = ["PD-L1", "miR-200","SNAI1","let-7","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7",""]
yticklabels = [" ", "miR-200","SNAI1","let-7","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7","AR"]

f, ax = plt.subplots(figsize=(9,9))
sns.heatmap(df_mean_mat, mask=mask,annot=annot_lis,annot_kws={'size': 22},fmt = '', vmax=1,vmin=-1, center=0,yticklabels=yticklabels,
                xticklabels=xticklabels,cmap=cmap,square=True, linewidths=.5, cbar_kws={"shrink": .7})
#plot_cor_matrix(df_mean_mat,mask)
plt.yticks(rotation=45)
plt.xticks(rotation=45)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
#plt.title('Correlation Plot with significantly important data')
plt.savefig('Correlation Plot with significantly important data 2_Network2')
#plt.show()
