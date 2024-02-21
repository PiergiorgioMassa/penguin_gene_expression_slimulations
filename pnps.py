#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import statistics as stats

#import data
mod=["PSMnWF","GEMnWF","PSMWF","GEMWF"]
psmod=["PSMnWF","PSMWF"] #pop size models
psmod_ext=["Pop Size Effect nonWF Model","Pop Size Effect WF Model"]
gemod=["GEMnWF","GEMWF"] #gene expression models
gemod_ext=["Gene Expression Effect nonWF Model","Gene Expression Effect WF Model"]
K=[1000,10000,100000]
rep=[1,2,3]
ds = {}
counter=0
for i in psmod: #import genome stats data from PS models
    for j in K:
        for k in rep:
            ds[i + "_K{0}".format(j) + "_rep{0}".format(k)] = pd.read_csv(i + "_K{0}".format(j) + "_seed{0}".format(k) + '.csv')
            ds[i + "_K{0}".format(j) + "_rep{0}".format(k)]["mod"]=i
            ds[i + "_K{0}".format(j) + "_rep{0}".format(k)]["N"]=j
            ds[i + "_K{0}".format(j) + "_rep{0}".format(k)]["rep"]="rep{0}".format(k)
            counter=counter+1
        ds[i + "_K{0}".format(j)]=pd.concat([ds[i + "_K{0}".format(j) + "_rep{0}".format(rep[0])],ds[i + "_K{0}".format(j) + "_rep{0}".format(rep[1])],ds[i + "_K{0}".format(j) + "_rep{0}".format(rep[2])]])
counter=0
for i in gemod: #import gene stats data from GE models
    for j in K:
        for k in rep:
            ds[i + "_K{0}".format(j) + "_rep{0}".format(k)] = pd.read_csv(i + "_geneStats_K{0}".format(j) + '_seed{0}'.format(k) + '.csv')
            ds[i + "_K{0}".format(j) + "_rep{0}".format(k)]["mod"]=i
            ds[i + "_K{0}".format(j) + "_rep{0}".format(k)]["N"]=j
            ds[i + "_K{0}".format(j) + "_rep{0}".format(k)]["Rep"]="rep{0}".format(k)
            counter=counter+1
        ds[i + "_K{0}".format(j)]=pd.concat([ds[i + "_K{0}".format(j) + "_rep{0}".format(rep[0])],ds[i + "_K{0}".format(j) + "_rep{0}".format(rep[1])],ds[i + "_K{0}".format(j) + "_rep{0}".format(rep[2])]])

#gene stats dataset correction: remove genes with pnps==0/0 and pnps==n/0, replace genes with pnps==0/s with 0
for i in gemod:
    for j in K:
        ds[i + "_K{0}".format(j)]=ds[i + "_K{0}".format(j)][ds[i + "_K{0}".format(j)]["pnps"]!="0/0"]
        ds[i + "_K{0}".format(j)].pnps[ds[i + "_K{0}".format(j)].pnps=="0/s"] = 0.0
        ds[i + "_K{0}".format(j)]=ds[i + "_K{0}".format(j)][ds[i + "_K{0}".format(j)]["pnps"]!="n/0"]
        ds[i + "_K{0}".format(j)]=ds[i + "_K{0}".format(j)].reset_index()
        ds[i + "_K{0}".format(j)]["pnps"] = ds[i + "_K{0}".format(j)]["pnps"].astype(float)

#create a single dataframe for each model
for i in mod:
    ds[i]=pd.concat([ds[i + "_K{0}".format(K[0])],ds[i + "_K{0}".format(K[1])],ds[i + "_K{0}".format(K[2])]])

#figure 4
fig, axes = plt.subplots(2,2,sharey=False)
fig.set_size_inches(21,14)
counter=0
for i in psmod:
    sns.boxplot(data=ds[i], x="N", y="pnps", linewidth=1.5, showfliers=False, palette="Greys", ax=axes[0,counter]).set_title(psmod_ext[counter], size=23, y=1.02)
    axes[0,counter].set_xlabel('N', size=15)
    axes[0,counter].set_ylabel('pN/pS', size=15)
    axes[0,counter].set_ylim(0.0,2.2)
    counter=counter+1
counter=0
for i in gemod:
    sns.boxplot(data=ds[i + "_K100000"], x="selCoeff", y="pnps", order=[-0.001,-0.01,-0.1], linewidth=1.5, showfliers=False, palette="Greys", ax=axes[1,counter]).set_title(gemod_ext[counter], size=23, y=1.02)
    axes[1,counter].set_xlabel('Selection coefficient', size=15)
    axes[1,counter].set_ylabel('pN/pS', size=15)
    axes[1,counter].set_ylim(0,2.2)
    counter=counter+1
plt.savefig('Fig4.pdf') 

#supplemetary figure 8 
suppfig_plottitl=["nonWF","WF"]
fig, ax = plt.subplots(2,2,sharey=True)
fig.set_size_inches(21,14)
counter=0
for i in gemod:
    sns.boxplot(data=ds[i], x="N", y="pnps", hue="selCoeff", hue_order=[-0.001,-0.01,-0.1], showfliers=False, linewidth=1.5, palette="Greys", ax=ax[0,counter]).set_title(suppfig_plottitl[counter], size=23, y=1.02)
    ax[0,counter].set_xlabel('Population size', size=15)
    ax[0,counter].set_ylim(0,4.0)
    sns.boxplot(data=ds[i], x="selCoeff", y="pnps", hue="N", order=[-0.001,-0.01,-0.1], showfliers=False, palette="Greys", ax=ax[1,counter])
    ax[1,counter].set_xlabel('Selection coefficient', size=15)
    ax[1,counter].set_ylim(0,4.0)
    counter=counter+1
ax[0,0].set_ylabel('pN/pS', size=15)
ax[0,1].set_ylabel('')
ax[1,0].set_ylabel('pN/pS', size=15)
ax[1,1].set_ylabel('')
ax[0,0].legend().remove()
ax[0,1].legend(loc="upper right", fontsize=12, title="Selection coefficient", title_fontsize=15)
ax[1,0].legend().remove()
ax[1,1].legend(loc="upper right", fontsize=12, title="Population size", title_fontsize=15)
plt.savefig('SuppFig8.pdf') 
