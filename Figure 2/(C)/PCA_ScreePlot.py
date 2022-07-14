#import all important packages
from sklearn.decomposition import PCA
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mb 
import pandas as pd
import seaborn as sns; sns.set_theme(color_codes=True)
import scipy
from scipy import stats
from matplotlib.patches import Rectangle

 
#input the run datafile
dfz = pd.read_table ('TS_RunA_10_4.dat', names = ["Model No.","Stable States", "ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"], usecols=["ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"])
df= pd.DataFrame(stats.zscore(dfz), columns = ["ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"])

#print(df)
n_components = 9
 
# Do the PCA.
pca = PCA(n_components=n_components)
reduced = pca.fit_transform(df[["ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"]])

# Append the principle components for each entry to the dataframe
for i in range(0, n_components):
    df['PC' + str(i + 1)] = reduced[:, i]

#display(df.head())

df.head()
# Do a scree plot
ind = np.arange(1, n_components+1)
(fig, ax) = plt.subplots(figsize=(8, 6))
sns.pointplot(x=ind, y=pca.explained_variance_ratio_)
#ax.set_title('Scree plot')
ax.set_xticks(ind)
ax.set_xticklabels(ind)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
ax.set_xlabel('PCA Component Number',fontsize=18)
ax.set_ylabel('Explained Variance',fontsize=18)
plt.savefig('PCA_Scree Plot1.png')

#plt.show()

# Show the points in terms of the first two PCs
'''g = sns.lmplot('PC1',
               'PC2',
               hue='species',data=df,
               fit_reg=False,
               scatter=True,
               height=7)

plt.show()
'''
# Plot a variable factor map for the first two dimensions.
kwargs1 ={'color':'red' }
kwargs2 ={'color':'black'}
kwargs3 ={'color':'orangered' }
kwargs4 ={'color':'peru'}
kwargs5 ={'color':'gold' }
kwargs6 ={'color':'greenyellow'}
kwargs7 ={'color':'teal' }
kwargs8 ={'color':'slategray'}
kwargs0 ={'color':'magenta' }
colour_tag =['red' ,'black','orangered','peru','gold','greenyellow','teal','slategray','magenta' ]
label_tag = ["ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"]

(fig, ax) = plt.subplots(figsize=(8, 8))
for i in range(0, pca.components_.shape[1]):
    print(i)
    '''ax.arrow(0,
             0,  # Start the arrow at the origin
             pca.components_[0, i],  #0 for PC1
             pca.components_[1, i],  #1 for PC2
             head_width=0.1,
                 head_length=0.1)'''
    if(i==0):
        ax.arrow( 0,
             0,  # Start the arrow at the origin
             pca.components_[0, i],  #0 for PC1
             pca.components_[1, i],  #1 for PC2
             head_width=0.1,
             head_length=0.1, **kwargs0 )
        Arrow0 = Rectangle((0, 0), 1, 1, label=label_tag[i], fc= colour_tag[i])
    if(i==1):
        ax.arrow( 0,
             0,  # Start the arrow at the origin
             pca.components_[0, i],  #0 for PC1
             pca.components_[1, i],  #1 for PC2
             head_width=0.1,
             head_length=0.1, **kwargs1 )
        Arrow1 = Rectangle((0, 0), 1, 1, label=label_tag[i], fc= colour_tag[i])
    if(i==2):
        ax.arrow( 0,
             0,  # Start the arrow at the origin
             pca.components_[0, i],  #0 for PC1
             pca.components_[1, i],  #1 for PC2
             head_width=0.1,
             head_length=0.1, **kwargs2 )
        Arrow2 = Rectangle((0, 0), 1, 1, label=label_tag[i], fc= colour_tag[i])
    if(i==3):
        ax.arrow( 0,
             0,  # Start the arrow at the origin
             pca.components_[0, i],  #0 for PC1
             pca.components_[1, i],  #1 for PC2
             head_width=0.1,
             head_length=0.1, **kwargs3 )
        Arrow3 = Rectangle((0, 0), 1, 1, label=label_tag[i], fc=colour_tag[i])
    if(i==4):
        ax.arrow( 0,
             0,  # Start the arrow at the origin
             pca.components_[0, i],  #0 for PC1
             pca.components_[1, i],  #1 for PC2
             head_width=0.1,
             head_length=0.1, **kwargs4 )
        Arrow4 = Rectangle((0, 0), 1, 1, label=label_tag[i], fc=colour_tag[i])
    if(i==5):
        ax.arrow( 0,
             0,  # Start the arrow at the origin
             pca.components_[0, i],  #0 for PC1
             pca.components_[1, i],  #1 for PC2
             head_width=0.1,
             head_length=0.1, **kwargs5 )
        Arrow5 = Rectangle((0, 0), 1, 1, label=label_tag[i], fc=colour_tag[i])
    if(i==6):
        ax.arrow( 0,
             0,  # Start the arrow at the origin
             pca.components_[0, i],  #0 for PC1
             pca.components_[1, i],  #1 for PC2
             head_width=0.1,
             head_length=0.1, **kwargs6 )
        Arrow6 = Rectangle((0, 0), 1, 1, label=label_tag[i], fc=colour_tag[i])
    if(i==7):
        ax.arrow( 0,
             0,  # Start the arrow at the origin
             pca.components_[0, i],  #0 for PC1
             pca.components_[1, i],  #1 for PC2
             head_width=0.1,
             head_length=0.1, **kwargs7 )
        Arrow7 = Rectangle((0, 0), 1, 1, label=label_tag[i], fc=colour_tag[i])
    if(i==8):
        ax.arrow( 0,
             0,  # Start the arrow at the origin
             pca.components_[0, i],  #0 for PC1
             pca.components_[1, i],  #1 for PC2
             head_width=0.1,
             head_length=0.1, **kwargs8 )
        Arrow8 = Rectangle((0, 0), 1, 1, label=label_tag[i], fc=colour_tag[i])
    

    
ax.legend( handles=[Arrow0,Arrow1, Arrow2,Arrow3, Arrow4,Arrow5, Arrow6,Arrow7, Arrow8], loc='upper left' )
an = np.linspace(0, 2 * np.pi, 100)
plt.plot(np.cos(an), np.sin(an))  # Add a unit circle for scale
plt.axis('equal')
ax.set_title('PCA Circle for Network')
#plt.show()
#plt.close()

label_tag = ["ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"]
kwargs2 ={'color':'royalblue'}
kwargs ={'fontsize':'16'}

(fig, ax) = plt.subplots(figsize=(8, 8))
for i in range(0, pca.components_.shape[1]):
    print(i)
    font = {'size': 16}
    ax.arrow(0,
             0,  # Start the arrow at the origin
             pca.components_[0, i],  #0 for PC1
             pca.components_[1, i],  #1 for PC2
             head_width=0.1,
             head_length=0.1,**kwargs2)
    if(i!=7 and i!=0 and i!=1 and i!=2 and i!=5 and i!=3) :
        plt.text(pca.components_[0, i] -0.4,
             pca.components_[1, i] +0.01,
             df.columns.values[i],**kwargs)
    if(i==7):
        plt.text(pca.components_[0, i] -0.5,
             pca.components_[1, i] -0.08,
             df.columns.values[i],**kwargs)
    if(i==0):
         plt.text(pca.components_[0, i] +0.08,
             pca.components_[1, i] +0.01,
             df.columns.values[i],**kwargs)

    if(i==1):
         plt.text(pca.components_[0, i] +0.06,
             pca.components_[1, i] +0.01,
             df.columns.values[i],**kwargs)

    if(i==2):
         plt.text(pca.components_[0, i] +0.06,
             pca.components_[1, i] +0.01,
             df.columns.values[i],**kwargs)
    if(i==5):
         plt.text(pca.components_[0, i] +0.06,
             pca.components_[1, i] +0.01,
             df.columns.values[i],**kwargs)
    if(i==3):
         plt.text(pca.components_[0, i] +0.06,
             pca.components_[1, i] +0.01,
             df.columns.values[i],**kwargs)
         

    
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.xlim([-1,1])
plt.ylim([-1,1])
an = np.linspace(0, 2 * np.pi, 100)
plt.plot(np.cos(an), np.sin(an))  # Add a unit circle for scale

plt.axis('equal')
plt.savefig('PCA Circle for Network.png')
plt.show()
