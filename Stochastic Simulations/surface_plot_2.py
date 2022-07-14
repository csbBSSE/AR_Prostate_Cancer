# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 23:11:10 2021

@author: csb
"""


from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.neighbors import KernelDensity
from scipy import stats
import plotly.graph_objects as go
import os

def landscape(dfTraj,nam, pthresh=-3, output_flder="./"): 
    # Retrieving the scores
    m1 = dfTraj["EMTScore"].to_numpy()
    m2 = dfTraj["Sensitivity"].to_numpy()
    
    xmin = m1.min();xmax = m1.max()
    ymin = m2.min();ymax = m2.max()
    ##  Calculating KDE
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([m1, m2])
    kernel = stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)
    Z=np.log10(Z)
    Z[Z<pthresh]=np.nan
    ## getting color axis
    colorAxis =1*Z
    colorAxis[X > 0.3, ] = 1
    colorAxis[X < 0.3, ] =-1
    colorAxis[X == 0.3]  = 0    
    

    
    #%% plotly figure
    fig = go.Figure(data=[go.Surface(z=-Z, x=X, y=Y)])
    
    fig.update_traces(surfacecolor = colorAxis,showscale=False) 
    fig.update_traces(contours_z=dict(show=True, usecolormap=False,highlightcolor="#ff0000", project_z=True))

    fig.update_layout(title=nam, xaxis_title="EMT Score", yaxis_title="Sensitivty",autosize=False,width=1000, height=1000,margin=dict(l=65, r=50, b=65, t=90))    
    fig.show()
    fig.write_image(output_flder+nam+"_landscape.png")
    fig.write_html(output_flder+nam+"_landscape.html")

folder='./input_files/100_0.75/'
filenames = [file for file in os.listdir(folder) if file.endswith('.csv')]

for filename in filenames:
    
    parID=filename.replace('_simDat.csv', '')
    print(parID)
    df = pd.read_csv(folder + filename)
    landscape(df,parID,pthresh=-2.5,output_flder=folder)     



# #%% matpotlib figure:
# # Creating figure
# fig = plt.figure(figsize =(14, 9))
# ax = plt.axes(projection ='3d')
 
# # Creating plot
# ax.plot_surface(X, Y, -Z)
 
# # show plot
# plt.show()
