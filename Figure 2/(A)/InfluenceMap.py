import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import os
from pathlib import Path
import pickle
import matplotlib.ticker as mticker

def Influ(out, top, nodes, intermat):
    '''Added new varible intermat here'''

    if not os.path.exists(Path(out)):
        os.mkdir(Path(out))

    for n in range(11,12):
        #intermat = pickle.load(open('intermat.f', 'rb')) '''change the intermat.f, to TOPO.TS'''
        '''You need to pass intermat here'''
        intermat = intermat
        '''This line is redundant, but I wrote to get you a better idea on what variable is passed'''
        N = float(np.count_nonzero(intermat))
        inf = intermat.copy()
        M_max = intermat.copy()
        M_max[M_max != 0] = 1.0
        for i in range(2,n):
            a = np.linalg.matrix_power(intermat, i).astype(float)
            b = np.linalg.matrix_power(M_max, i).astype(float)
            inf = inf + np.divide(a, b, out=np.zeros_like(a), where=b!=0)

        inf = inf/n

        data = pd.DataFrame(data=inf,index=nodes,columns=nodes)
        data = data.loc[top,top]
        inf = np.array(data)
        print(inf)
        if n == 11:
            with open("bool_influ.data",'wb') as f:
                pickle.dump(inf,f)

        fig = plt.figure(figsize=(8,8))
        ax1 = fig.add_subplot(111)
        plt.imshow(inf, cmap='seismic', interpolation='nearest')
        #plt.clim(-1,1)
        cb = plt.colorbar()
        
        ax1.set_yticks(np.arange(len(top)))
        ax1.set_xticks(np.arange(-0.5,len(top),1), minor=True);
        ax1.set_yticks(np.arange(-0.5,len(top),1), minor=True);
        ax1.set_yticklabels(labels=top, minor=False, font='Arial', fontsize = 14 )
        ticks_loc = ax1.get_xticks().tolist()
        ticks_loc = ticks_loc[1:10]
        print(ticks_loc)
        ax1.xaxis.set_major_locator(mticker.FixedLocator(ticks_loc))
        ax1.set_xticklabels(labels = top, rotation = 90, font='Arial', fontsize = 14)

        ax1.grid(which='minor', color='w', linestyle='-', linewidth=0.2)
        #plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.savefig((Path(out,'influence_n={}.jpeg'.format(n))), dpi=300)
        plt.close()

        #################################################################################
        '''You don't need this part now'''
        '''STATES = pickle.load(open('states.f', 'rb'))
        freq = STATES[:,-1]
        stable_vects = STATES[:,:-1]
        stable_vects[stable_vects < 0] = 0
        stable_vects = stable_vects.astype('int')

        INTERMAT = pd.DataFrame(inf.T, index=top, columns=top)

        arr = np.array([[0]*len(top)]*len(top))
        for run, vects in enumerate(stable_vects):
            for num_i, i in enumerate(top):
                for num_j, j in enumerate(top):
                    if INTERMAT[i][j] > 0:
                        arr[num_i][num_j] += 1*freq[run] if vects[nodes.index(i)] == vects[nodes.index(j)] else -1*freq[run]
                    if INTERMAT[i][j] < 0:
                        arr[num_i][num_j] += 1*freq[run] if vects[nodes.index(i)] != vects[nodes.index(j)] else -1*freq[run]

        arr = arr/sum(freq)
        fig = plt.figure(figsize=(8,8))
        ax1 = fig.add_subplot(111)
        plt.imshow(arr, cmap='seismic', interpolation='nearest')
        plt.clim(-1,1)
        plt.colorbar()
        ax1.set_yticks(np.arange(len(top)))
        ax1.set_xticks(np.arange(-0.5,len(top),1), minor=True);
        ax1.set_yticks(np.arange(-0.5,len(top),1), minor=True);
        ax1.set_yticklabels(labels=np.arange(1,len(top)+1), minor=False)
        ax1.grid(which='minor', color='w', linestyle='-', linewidth=0.2)
        plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.savefig((Path(out,'influ_with_states_n={}.jpeg'.format(n))), dpi=300)
        plt.close()'''

if __name__ == '__main__':

    '''Pass nodes and interaction matrix'''
    '''Here I wrote a function which reads .ids and .topo file'''
    '''Write your own functions'''
    '''I am writing here'''

    with open("TS.topo", "r") as f:
        lines  = f.readlines()[1:]
        nodes = []
        for line in lines:
            strings = line.strip("\n").split("\t")
            nodes.append(strings[0])
            nodes.append(strings[1])

        #Take unique and sorted nodes
        nodes = sorted(list(set(nodes)))

        #Initiate an array
        intermat = np.zeros((len(nodes), len(nodes)))
        for line in lines:
            res = line.strip("\n").split("\t")
            if res[2] == '1': #activating
                intermat[nodes.index(res[1])][nodes.index(res[0])] = 1
            if res[2] == '2':
                intermat[nodes.index(res[1])][nodes.index(res[0])] = -1

    print(nodes)
    # After running the code once, I hab a look at influence matrix and wrote top accordingly
    top = ['AR', 'AR-v7', 'LIN28', 'SLUG', 'hnRNPA1', 'ZEB1', 'SNAI1', 'let-7', 'miR-200']
    Influ("Influence",top,nodes,intermat)
