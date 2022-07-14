#Version 2: 23/06/2021
#By Abheepsa

#this code takes a RACIPE solution file, calculates EMT, Stemness and Resistance metascores
#and plots the histograms of these scores fit with Gaussian distributions


import pandas as pd
import numpy as np
import seaborn as sns
from scipy import stats
import matplotlib.pyplot as plt

df = pd.read_table('TS_runA_10_4.dat', names = ["Model No.","Stable States", "ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"], usecols=["ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"])
df.apply(stats.zscore)

# here we define a metascore for EMTness such that EMT score = ZEB1 - miR200
def emtscoring1 (gene_data):
    zeb1 = np.array(gene_data['ZEB1'])
    miR200= np.array(gene_data['miR-200'])
    emtscore= np.subtract(zeb1, miR200)
    return emtscore

#here we define a metascore for EMTness such that EMT score = (ZEB1+SLUG+SNAI1)/3-miR200
def emtscoring2 (gene_data):
    zeb1 = np.array(gene_data['ZEB1'])
    miR200= np.array(gene_data['miR-200'])
    slug = np.array(gene_data['SLUG'])
    snail = np.array(gene_data['SNAI1'])
    emtscore= np.subtract(np.add(zeb1, slug, snail)/3, miR200)
    return emtscore

#here we define a metascore for EMTness such that EMT score = (ZEB1+SLUG)/2-miR200
def emtscoring3 (gene_data):
    zeb1 = np.array(gene_data['ZEB1'])
    miR200= np.array(gene_data['miR-200'])
    slug = np.array(gene_data['SLUG'])
    emtscore= np.subtract(np.add(zeb1, slug)/2, miR200)
    return emtscore

#here we define a metascore for Stemness such that Stemness Score = Lin28-Let7
def stemnessscoring (gene_data):
    lin28 = np.array(gene_data['LIN28'])
    let7= np.array(gene_data['let-7'])
    stemnessscore= np.subtract(lin28, let7)
    return stemnessscore

#here we define a metascore for resistance such that Resistance Score = AR+ AR-v7
def resistancescoring (gene_data):
    arv7 = np.array(gene_data['AR-v7'])
    ar= np.array(gene_data['AR'])
    resistancescore= np.add(ar, arv7)
    return resistancescore

emtscore1 = stats.zscore(emtscoring1(df))
'''emtscore2 = emtscoring2(df)
emtscore3 = emtscoring3(df)'''
stemnessscore= stemnessscoring(df)
resistancescore= stats.zscore(resistancescoring(df))

np.savetxt('EMTScore1.dat', emtscore1, delimiter=',')
#np.savetxt('EMTScore2.dat', emtscore2, delimiter=',')
np.savetxt('StemnessScore.dat', stemnessscore, delimiter=',')
np.savetxt('ResistanceScore.dat', resistancescore, delimiter=',')

plt.hist(x = emtscore1, bins = 40, density = True)
plt.xticks(font='Arial', fontsize=12)
plt.yticks(font='Arial', fontsize=12)
plt.ylabel('Density', font = 'Arial', fontsize = 16)
plt.xlabel('EM score = ZEB1 - miR200', font = 'Arial', fontsize = 16)
plt.title('Histogram of EM Score', font = 'Arial', fontsize = 19)
plt.savefig('Histogram EM Scores.png')
plt.clf()
'''
plt.hist(x = emtscore2, bins = 40, density=True)
plt.xlabel('EMT score = (ZEB1+SLUG+SNAI1)/3-miR200')
plt.title('Histogram of EMT Score2')
plt.savefig('Histogram EMT Scores 2.png')
plt.clf()

plt.hist(x = emtscore3, bins = 40, density=True)
plt.xlabel('EMT score = (ZEB1+SLUG)/2-miR200')
plt.title('Histogram of EMT Score3')
plt.savefig('AR+ARv7 KD Histogram EMT Scores 3.png')
plt.clf()

plt.hist(x = stemnessscore, bins = 40, density=True)
plt.xlabel('Stemness Score = LIN28 - let7')
plt.title('Histogram of Stemness Score')
plt.savefig('Histogram of Stemness Scores.png')
plt.clf()'''

plt.hist(x = resistancescore, bins = 40, density = True)
plt.xticks(font='Arial', fontsize=14)
plt.yticks(font='Arial', fontsize=14)
plt.xlabel('Resistance Score = AR + AR-v7', font = 'Arial', fontsize = 16)
plt.ylabel('Density', font = 'Arial', fontsize = 16)
plt.title('Histogram of Resistance Score', font = 'Arial', fontsize = 19)
plt.savefig('Histogram of Resistance Scores.png')
plt.clf()