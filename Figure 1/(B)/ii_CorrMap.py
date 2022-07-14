from string import ascii_letters
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import statistics
import fontstyle


'''df_n are the dataframes of n different runs, where the columns are in the order of "Model No.","Stable States","let-7" , "miR-200","SNAI1", "AR","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7"'''

sns.set_theme(style="white")
df_1 = pd.read_table ('TS_RunA_10_4.dat', names = ["Model No.","Stable States","let-7" , "miR-200","SNAI1", "AR","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7"], usecols=["let-7", "miR-200","SNAI1", "AR","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7"])
#print(df_1)
corr1 = df_1.corr()
#print(corr1)

df_2 = pd.read_table ('TS_RunB_10_4.dat', names = ["Model No.","Stable States","let-7" , "miR-200","SNAI1", "AR","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7"], usecols=["let-7", "miR-200","SNAI1", "AR","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7"])
corr2 = df_2.corr()
#print(corr2)

df_3 = pd.read_table ('TS_RunC_10_4.dat', names = ["Model No.","Stable States","let-7" , "miR-200","SNAI1", "AR","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7"], usecols=["let-7", "miR-200","SNAI1", "AR","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7"])
corr3 = df_3.corr()
#print(corr3)

corrmat1 = corr1.to_numpy()
corrmat2 = corr2.to_numpy()
corrmat3 = corr3.to_numpy()

#print(corrmat3)
#print(corrmat3[2,2])

a = (corrmat1[2,2], corrmat2[2,2], corrmat3[2,2])
b= statistics.mean(a)
c =statistics.stdev(a)
#print(c)

mean_mat = np.zeros(shape=(9, 9))
stdev_mat = np.zeros(shape=(9, 9))


for i in range (0,9):
    for j in range (0,9):
        a = (corrmat1[i,j], corrmat2[i,j], corrmat3[i,j])
        b= statistics.mean(a)
        mean_mat[i,j] = b
        c =statistics.stdev(a)
        stdev_mat[i,j] = c
        
#print(mean_mat)

annot_mat = mean_mat
annot_matrix = np.round_(annot_mat, decimals = 2)
#print (annot_matrix[0,8])
annot_lis = annot_matrix.tolist()
#print(annot_lis)
var = annot_lis[0][8]
annot_lis[0][8] = str(var) + '*'
#print(annot_lis[0][8])


'''
xlsx = pd.ExcelFile('D:/systems biology/Summer Project/Review/AR work begins/Initial/datafornoyellowARandSNAIsi.xlsx')
call = pd.read_excel(xlsx, 'abcd')'''


def corr_sig(df=None):
    p_matrix = np.zeros(shape=(df.shape[1],df.shape[1]))
    for col in df.columns:
        for col2 in df.drop(col,axis=1).columns:
            _ , p = stats.pearsonr(df[col],df[col2])
            p_matrix[df.columns.to_list().index(col),df.columns.to_list().index(col2)] = p
    return p_matrix

#p_values = corr_sig(df)
#mask = np.invert(np.tril(p_values<0.05))
# note seaborn will hide correlation were the boolean value is True in the mask

# Compute the correlation matrix


# Generate a mask for the upper triangle
mask = np.triu(np.ones_like(mean_mat))

# Set up the matplotlib figure
f, ax = plt.subplots(figsize=(8,8))

# Generate a custom diverging colormap
cmap = sns.diverging_palette(230, 20, as_cmap=True)

# Draw the heatmap with the mask and correct aspect ratio
xticklabels = ["let-7", "miR-200","SNAI1", "AR","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7"]
yticklabels = ["let-7", "miR-200","SNAI1", "AR","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7"]
plt.yticks(fontsize=14)
sns.heatmap(mean_mat, mask=mask, cmap=cmap, vmax=.3, center=0,annot=True,yticklabels=yticklabels,
            xticklabels=xticklabels, square=True, linewidths=.5, cbar_kws={"shrink": .5})

#plt.savefig('Correlation Map with Mean value of three runs')
#plt.show()
#plt.close()

# Generate a mask for the upper triangle
mask = np.triu(np.ones_like(mean_mat))

# Set up the matplotlib figure
f, ax = plt.subplots(figsize=(8,8))

# Generate a custom diverging colormap
cmap = sns.diverging_palette(230, 20, as_cmap=True)

# Draw the heatmap with the mask and correct aspect ratio
xticklabels = ["let-7", "miR-200","SNAI1", "AR","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7"]
yticklabels = ["let-7", "miR-200","SNAI1", "AR","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7"]
sns.heatmap(stdev_mat, mask=mask, cmap=cmap, vmax=.3, center=0,annot=True,yticklabels=yticklabels,
            xticklabels=xticklabels,square=True, linewidths=.5, cbar_kws={"shrink": .5})
#plt.savefig('Standard Deviation Map 1')
#plt.show()
#plt.close()

df_mean_mat = pd.DataFrame(mean_mat)
#print(df_mean_mat)

df_stdev_mat = pd.DataFrame(stdev_mat)

p_values = corr_sig(df_mean_mat)                     # get p-Value
#mask = np.invert(np.tril(p_values<0.05))    # mask - only get significant corr
f, ax = plt.subplots(figsize=(8,8))
cmap = sns.diverging_palette(230, 20, as_cmap=True)
xticklabels = ["let-7", "miR-200","SNAI1", "AR","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7"]
yticklabels = ["let-7", "miR-200","SNAI1", "AR","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7"]
sns.heatmap(p_values, mask=mask, cmap=cmap, vmax=.3, center=0,annot=True,yticklabels=yticklabels,
            xticklabels=xticklabels, square=True, linewidths=.5, cbar_kws={"shrink": .5})
#plt.savefig('p-value of the three runs 1')
#plt.show()
#plt.close()

#print(type(p_values))



for i in range (0,9):
    for j in range (0,9):
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
#plt.savefig('Correlation Plot with significantly important data 1')
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
#plt.savefig('Standard Deviation for the correlation plot ')
#plt.show()



for i in range (0,9):
    for j in range (0,9):
        if p_values[i,j] > 0.05:
            var = annot_lis[i][j]
            annot_lis[i][j] = 'X'
        else:
            var = annot_lis[i][j]
            annot_lis[i][j] = ' '

xticklabels = ["let-7", "miR-200","SNAI1", "AR","SLUG","ZEB1","LIN28","hnRNPA1"]
yticklabels = [" ","miR-200","SNAI1", "AR","SLUG","ZEB1","LIN28","hnRNPA1","AR-v7"]

f, ax = plt.subplots(figsize=(9,9))
sns.heatmap(df_mean_mat, mask=mask,annot=annot_lis,annot_kws={'size': 22},fmt = '', vmax=1,vmin=-1, center=0,yticklabels=yticklabels,
                xticklabels=xticklabels,cmap=cmap,square=True, linewidths=.5, cbar_kws={"shrink": .7})
#plot_cor_matrix(df_mean_mat,mask)
plt.yticks(rotation=45)
plt.xticks(rotation=45)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
#plt.title('Correlation Plot with significantly important data')
plt.savefig('Correlation Plot with significantly important data 2')
#plt.show()
