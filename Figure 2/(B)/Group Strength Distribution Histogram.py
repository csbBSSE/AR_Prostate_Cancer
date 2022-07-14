import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import os
from pathlib import Path
from scipy.stats import zscore

def emtscoring (gene_data):
    zeb1 = np.array(gene_data['ZEB1'])
    miR200= np.array(gene_data['miR-200'])
    emtscore= np.subtract(zeb1, miR200)
    
    return emtscore

def resistancescoring (gene_data):
    arv7 = np.array(gene_data['AR-v7'])
    ar= np.array(gene_data['AR'])
    resistancescore= np.add(ar, arv7)
    return resistancescore

def emtvsrescorr (emt, resistance):
    corr, p = pearsonr(emt, resistance)
    return corr, p

def Interactionmatrix(n):
    filename = "TSRand_"+str(n)+".topo"
    topol = pd.read_table(filename)
    genes = pd.unique(topol.loc[:, 'Source'])
    num_genes = genes.size
    length = topol.shape[0]
    intmat = pd.DataFrame(data = np.zeros((num_genes, num_genes), dtype=np.int32), index = genes, columns= genes)
    for i in range (0, length):
        source = topol.loc[i,'Source']
        targ = topol.loc[i, 'Target']
        int_type= topol.loc[i,'Type']
        intmat.loc[source, targ]= 1 if int_type==1 else -1
    return intmat.to_numpy(), genes


def TeamStrength(out, top, nodes, intermat):
    if not os.path.exists(Path(out)):
        os.mkdir(Path(out))
        
    for n in range(11,12):
        N = float(np.count_nonzero(intermat))
        inf = np.copy(intermat)
        M_max = np.copy(intermat)
        M_max[M_max != 0] = 1.0
        for i in range(2,n):
            a = np.linalg.matrix_power(intermat, i).astype(float)
            b = np.linalg.matrix_power(M_max, i).astype(float)
            inf = inf + np.divide(a, b, out=np.zeros_like(a), where=b!=0)

        inf = inf/n

        data = pd.DataFrame(data=inf,index=nodes,columns=nodes)
        data = data.loc[top,top]
        team1 = ['AR', 'AR-v7', 'LIN28', 'SLUG','hnRNPA1','ZEB1']
        team2 = ['SNAI1', 'let-7', 'miR-200']
        datared1 = np.array(data.loc[team1, team1])/(len(team1)*len(team1))
        datared2 = np.array(data.loc[team2, team2])/(len(team2)*len(team2))
        datablue1= np.array(data.loc[team1, team2])/(len(team1)*len(team2))
        datablue2= np.array(data.loc[team2, team1])/(len(team1)*len(team2))
        team_strength = np.sum(datared1)+np.sum(datared2)-np.sum(datablue1)-np.sum(datablue2)
        return team_strength

corrvector = [0]*100
pvector = [0]*100
teamstrength = [0]*100

dfbio = pd.read_table('TS_runA_10_4.dat', names = ["Model No.","Stable States", "ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"], usecols=["ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"])
emtbio = zscore(emtscoring(dfbio))
resbio = zscore(resistancescoring(dfbio))

corrbio, pbio = emtvsrescorr(emtbio, resbio)
print("CORRELATION:", corrbio, 'bla', pbio)

for i in range(1, 101):
    filename = "TSRand_"+str(i)+"_solution.dat"
    df = pd.read_table(filename, names = ["Model No.","Stable States", "ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"], usecols=["ZEB1", "miR-200","SNAI1", "AR","SLUG","let-7","LIN28","hnRNPA1","AR-v7"])
    df.apply(zscore)
    emt = zscore(emtscoring(df))
    res = zscore(resistancescoring(df))
    corrvector[i-1], pvector[i-1]= emtvsrescorr(emt, res)
    intermat, nodes = Interactionmatrix(i)
    top = ['AR', 'AR-v7', 'LIN28', 'SLUG','hnRNPA1','ZEB1', 'SNAI1', 'let-7', 'miR-200']
    teamstrength[i-1] = TeamStrength("Influence",top,nodes,intermat)

print(corrvector, pvector, teamstrength)
plt.scatter(corrvector, teamstrength)
#plt.scatter(corrbio, 19.47053206394525, color = 'red')
print(pearsonr(corrvector, teamstrength))
plt.xlabel('EMT-Resistance Correlation Coefficient')
plt.ylabel('Team Strength')
plt.title('Scatter: Group Strength vs. EMT-Res Correlation')
corr, p = pearsonr(corrvector, teamstrength)
text = "Corr = {:.6f}".format(corr)
plt.text(x = 0.4, y = -0.075, s = text)
plt.savefig('ScatterCorrectGroup')
print(pearsonr(corrvector, teamstrength))
#plt.show()