# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 22:49:34 2021

@author: csb
"""

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.neighbors import KernelDensity
def kde2D(x, y, bandwidth, xbins=100j, ybins=100j, **kwargs): 
    """Build 2D kernel density estimate (KDE)."""

    # create grid of sample locations (default: 100x100)
    xx, yy = np.mgrid[x.min():x.max():xbins, 
                      y.min():y.max():ybins]

    xy_sample = np.vstack([yy.ravel(), xx.ravel()]).T
    xy_train  = np.vstack([y, x]).T

    kde_skl = KernelDensity(bandwidth=bandwidth, **kwargs)
    kde_skl.fit(xy_train)

    # score_samples() returns the log-likelihood of the samples
    z = np.exp(kde_skl.score_samples(xy_sample))
    return xx, yy, np.reshape(z, xx.shape)



df = pd.read_csv('./input_files/100_0.5/118_simDat.csv')
 
# Creating dataset
x = df["EMTScore"].to_numpy()
y = df["Sensitivity"].to_numpy()
xx, yy, zz = kde2D(x, y, 1.0)

#plt.scatter(x, y, s=2, facecolor='white

#%%
# Creating figure
fig = plt.figure(figsize =(14, 9))
ax = plt.axes(projection ='3d')
 
# Creating plot
ax.plot_surface(xx, yy, zz)
 
# show plot
plt.show()