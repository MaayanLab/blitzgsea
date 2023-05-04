import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import blitzgsea as blitz

def running_sum(signature, geneset, library, result=None, compact=False):
    """
    Plot the running sum for a given geneset and signature.

    Parameters:
    signature (array-like): The gene expression signature to analyze.
    geneset (str): The name of the gene set for a gene set in the library.
    library (array-like): The gene set library to use for enrichment analysis.
    result (array-like, optional): A precomputed enrichment result. Default is None.
    compact (bool, optional): If True, return a compact representation of the running sum plot for better readability in small plots. Default is False.

    Returns:
    figure: The running sum plot for the given geneset and signature.
    """
    plt.ioff()
    signature = signature.sort_values(1, ascending=False).set_index(0)
    signature = signature[~signature.index.duplicated(keep='first')]
    
    signature_map = {}
    for i,h in enumerate(signature.index):
        signature_map[h] = i
    
    gs = set(library[geneset])
    hits = [i for i,x in enumerate(signature.index) if x in gs]
    
    running_sum, es = blitz.enrichment_score(np.array(np.abs(signature.iloc[:,0])), signature_map, gs)
    running_sum = list(running_sum)
    fig = plt.figure(figsize=(7,5))
    
    if compact:
        ax = fig.add_gridspec(5, 11, wspace=0, hspace=0)
        ax1 = fig.add_subplot(ax[0:4, 0:11])
    else:
        ax = fig.add_gridspec(12, 11, wspace=0, hspace=0)
        ax1 = fig.add_subplot(ax[0:7, 0:11])
  
    if compact:
        ax1.plot(list(running_sum), color=(0,1,0), lw=5)
        ax1.tick_params(labelsize=24)
    else:
        ax1.plot(list(running_sum), color=(0,1,0), lw=3)
        ax1.tick_params(labelsize=16)
    plt.xlim([0, len(running_sum)])
    
    nn = np.where(np.abs(running_sum)==np.max(np.abs(running_sum)))[0][0]
    ax1.vlines(x=nn, ymin=np.min(running_sum), ymax=np.max(running_sum),linestyle = ':', color="red")
    
    if not result is None:
        if es > 0:
            if compact:
                ax1.text(len(running_sum)/30, 0, "NES="+"{:.3f}".format(result.loc[geneset,"nes"]), size=25, bbox={'facecolor':'white','alpha':0.8,'edgecolor':'none','pad':1}, ha='left', va='bottom', zorder=100)
            else:
                ax1.text(len(running_sum)/30, 0, "NES="+"{:.3f}".format(result.loc[geneset,"nes"]), size=20, bbox={'facecolor':'white','alpha':0.8,'edgecolor':'none','pad':1}, ha='left', va='bottom', zorder=100)
        else:
            if compact:
                ax1.text(len(running_sum)/30, 0, "NES="+"{:.3f}".format(result.loc[geneset,"nes"]), size=25, bbox={'facecolor':'white','alpha':0.8,'edgecolor':'none','pad':1}, ha='left', va='top', zorder=100)
            else:
                ax1.text(len(running_sum)/30, 0, "NES="+"{:.3f}".format(result.loc[geneset,"nes"]), size=20, bbox={'facecolor':'white','alpha':0.8,'edgecolor':'none','pad':1}, ha='left', va='top', zorder=100)
    else:
        if es > 0:
            if compact:
                ax1.text(len(running_sum)/30, 0, "ES="+"{:.3f}".format(running_sum[nn]), size=25, bbox={'facecolor':'white','alpha':0.8,'edgecolor':'none','pad':1}, ha='left', va='bottom', zorder=100)
            else:
                ax1.text(len(running_sum)/30, 0, "ES="+"{:.3f}".format(running_sum[nn]), size=20, bbox={'facecolor':'white','alpha':0.8,'edgecolor':'none','pad':1}, ha='left', va='bottom', zorder=100)
        else:
            if compact:
                ax1.text(len(running_sum)/30, 0, "ES="+"{:.3f}".format(running_sum[nn]), size=25, bbox={'facecolor':'white','alpha':0.8,'edgecolor':'none','pad':1}, ha='left', va='top', zorder=100)
            else:
                ax1.text(len(running_sum)/30, 0, "ES="+"{:.3f}".format(running_sum[nn]), size=20, bbox={'facecolor':'white','alpha':0.8,'edgecolor':'none','pad':1}, ha='left', va='top', zorder=100)

    ax1.grid(True, which='both')
    ax1.set(xticks=[])
    plt.title(geneset, fontsize=18)
    
    if compact:
        plt.ylabel("ES", fontsize=24)
    else:
        plt.ylabel("Enrichment Score (ES)", fontsize=16)
    
    if compact:
        ax1 = fig.add_subplot(ax[4:, 0:11])
        rank_vec = signature[1]
        ax1.vlines(x=hits, ymin=0, ymax=1, color=(0,0,0,1), lw=1.5)
        plt.xlim([0, len(running_sum)])
        plt.ylim([0, 1])
        ax1.set(yticks=[])
        ax1.set(xticks=[])
        x = np.arange(0.0, len(rank_vec), 20).astype("int")
        x = np.append(x, signature.shape[0]-1)
        plt.xlim([0, len(running_sum)])
        plt.ylim([0, 1])
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        plt.xlabel("Rank", fontsize=24)
        
        posv = np.percentile(range(len(rank_vec[rank_vec > 0])), np.linspace(0,100,10))
        for i in range(9):
            plt.gca().add_patch(Rectangle((posv[i],0),posv[i+1]-posv[i],0.5,linewidth=0,facecolor='red', alpha=0.6*(1-i*0.1)))
        negv = np.percentile(range(len(rank_vec[rank_vec <= 0])), np.linspace(0,100,10))
        for i in range(9):
            plt.gca().add_patch(Rectangle((posv[-1]+negv[i],0),negv[i+1]-negv[i],0.5,linewidth=0,facecolor='blue', alpha=0.6*(0.1+i*0.1)))
    else:
        ax1 = fig.add_subplot(ax[7:8, 0:11])
        ax1.vlines(x=hits, ymin=-1, ymax=1, color=(0,0,0,1), lw=0.5)
        plt.xlim([0, len(running_sum)])
        plt.ylim([-1, 1])
        ax1.set(yticks=[])
        ax1.set(xticks=[])
        ax1 = fig.add_subplot(ax[8:, 0:11])
        rank_vec = signature[1]
        x = np.arange(0.0, len(rank_vec), 20).astype("int")
        x = np.append(x, signature.shape[0]-1)
        ax1.fill_between(x, np.array(rank_vec)[x], color="lightgrey")
        ax1.plot(x, np.array(rank_vec)[x], color=(0.2,0.2,0.2), lw=1)
        ax1.hlines(y=0, xmin=0, xmax=len(rank_vec), color="black", zorder=100, lw=0.6)
        plt.xlim([0, len(running_sum)])
        plt.ylim([np.min(rank_vec), np.max(rank_vec)])
        minabs = np.min(np.abs(rank_vec))
        zero_cross = int(np.where(np.abs(rank_vec)==minabs)[0][0])
        ax1.vlines(x=zero_cross, ymin=np.min(rank_vec), ymax=np.max(rank_vec),linestyle = ':',)
        ax1.text(zero_cross, np.max(rank_vec)/3, "Zero crosses at "+str(zero_cross), bbox={'facecolor':'white','alpha':0.5,'edgecolor':'none','pad':1}, ha='center', va='center')
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        plt.xlabel("Rank in Ordered Dataset", fontsize=16)
        plt.ylabel("Ranked list metric", fontsize=16)
        ax1.tick_params(labelsize=16)
    plt.ion()
    fig.patch.set_facecolor('white')
    return fig

def top_table(signature, library, result, n=10):
    """
    Plot a table to enrichment results for top N enriched gene sets for a given geneset and signature.

    Parameters:
    signature (array-like): The gene expression signature to analyze.
    library (array-like): The gene set library to use for enrichment analysis.
    result (array-like, optional): A precomputed enrichment result. Default is None.
    n (integer): number of top enriched gene sets to be plotted
    
    Returns:
    figure: The running sum plot for the given geneset and signature.
    """
    sig = signature.sort_values(1, ascending=False).set_index(0)
    sig = sig[~sig.index.duplicated(keep='first')]

    plt.ioff()
    fig = plt.figure(figsize=(5,0.5*n), frameon=False)
    ax = fig.add_subplot(111)
    fig.patch.set_visible(False)
    plt.axis('off')

    ax.vlines(x=[0.2,0.8], ymin=-0.1, ymax=1, color="black")
    ln = np.linspace(-0.1,1,n+1)[::-1]
    ax.hlines(y=ln, xmin=0, xmax=1, color="black")

    ax.text(0.03, 1.03, "NES", fontsize=16)
    ax.text(0.84, 1.03, "SET", fontsize=16)

    for i in range(n):
        ax.text(0.03, (ln[i]+ln[i+1])/2, "{:.3f}".format(result.iloc[i, 1]), verticalalignment='center')
        ax.text(0.84, (ln[i]+ln[i+1])/2, result.index[i], verticalalignment='center')
        
        gs = set(library[result.index[i]])
        hits = np.array([i for i,x in enumerate(sig.index) if x in gs])
        hits = (hits/len(sig.index))*0.6+0.2
        
        if result.iloc[i, 1] > 0:
            ax.vlines(hits, ymax=ln[i], ymin=ln[i+1], color="red", lw=0.5, alpha=0.3)
        else:
            ax.vlines(hits, ymax=ln[i], ymin=ln[i+1], color="blue", lw=0.5, alpha=0.3)
    fig.patch.set_facecolor('white')
    plt.ion()
    return fig
