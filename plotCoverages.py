
import numpy as np
import pandas as pd
import matplotlib.pylab as pl
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
# from sklearn.cluster import KMeans
import copy, os, pdb, random, shutil, subprocess, time
# import cPickle as pk
import smooth as sm
import detect_peaks as peak

import pdb
################################################################################
# plotCoverages.py
#
# Calculates and plots the coverage score for each chromosome for the given dataset
#
################################################################################

################################################################################
# main
################################################################################
pl.ioff()
saveDir = 'Coverage'
if not os.path.exists(saveDir):
    os.makedirs(saveDir)
sns.set_style('white')

def main():
    wn=5
    chrLens = pd.read_csv('../data/chromosomeLengths.txt',delimiter='\t',header=None)
    chrLens.rename(columns={0:'chrNo',1:'len'},inplace=True)
    if False:
        sampleList = ['AsyncWT','AFac','simData','WT1','fkh1','fkh2']
        coverages = dict()
        for smp in sampleList:
            intra = pd.read_csv('../data/'+smp+'_intra.txt',delimiter='\t',header= None)
            intra = intra.rename(columns={0:'chr1',1:'frg1',2:'chr2',3:'frg2',4:'freq'})
            covVec = [None]*16
            for ch in range(1,17):
                binMat = getBinMat(intra,chrLens,ch=ch)
                covVec[ch-1],_ = getCoverage(binMat,wn)

            coverages[smp] = np.concatenate(covVec)

    if False:
        ##### Plot coverage correlations as scatter

        pp = PdfPages(saveDir+'/'+'CoverageScatter.pdf')
        for i in range(len(sampleList)-1):
            for j in range(i+1,len(sampleList)):
                g = sns.jointplot(coverages[sampleList[i]],coverages[sampleList[j]],kind="reg",
                          xlim=(0, 1), ylim=(0, 1), color="r", size=7)
                g.set_axis_labels(sampleList[i],sampleList[j])
                pp.savefig()
        pp.close()


    ##### Plot coverage profiles

    for wn in [8,10]:
        print('performing for: '+ str(10*wn))
        cols_ = ['#d9544d','#3b5b92','#6f828a','red']
        # cols_ = ['blue','red','grey']
        # pp = PdfPages(saveDir+'/'+'CoverageProfiles_AsyncWT_AFac_Simulation_v2.pdf')
        # pp = PdfPages(saveDir+'/'+'CoverageProfiles_WT1_fkh1_fkh2_50kbSmooth.pdf')
        for ch in range(1,17):
        #     for smp in ['WT1','fkh1','fkh2']:
        #     fig=pl.figure()
            qq=0
            # for smp in ['AsyncWT','AFac','simData']:
            for smp in ['WT1']:
                intra = pd.read_csv('../data/'+smp+'_intra.txt',delimiter='\t',header= None)
                intra = intra.rename(columns={0:'chr1',1:'frg1',2:'chr2',3:'frg2',4:'freq'})
        #         intra=intra[intra[5]<5e-2]
                binMat = getBinMat(intra,chrLens,ch=ch)
                covVec,xran = getCoverage(binMat,wn)
                # pl.plot(xran*10,covVec,color=cols_[qq],ls='--',lw=0.5)
                covVec = sm.smooth(covVec,6)[5:-5]
                # pl.plot(xran*10,covVec,label=smp,color=cols_[qq])
                vls = peak.detect_peaks(np.r_[-1,covVec,-1],valley=True,show=False)
                pks = peak.detect_peaks(np.r_[-1,covVec,-1],valley=False,show=False)
                if len(vls)>0:
                    if max(pks)< max(vls):
                        pks = np.r_[pks,len(covVec)-1]
                        covVec[len(covVec)-1] = covVec[max(vls)]
                    if min(pks)> min(vls):
                        pks = np.r_[0,pks]
                        covVec[0]=covVec[min(vls)]
                    if len(pks)==(len(vls)+2):
                        pks = pks[1:]
                    mean_trophs = (covVec[pks[:-1]-1]-covVec[vls-1] + covVec[pks[1:]-1]-covVec[vls-1])/2
                    assert all(mean_trophs > 0), 'Peaks should have higher values than valleys. check the implementation'
                    tf = mean_trophs > 0.05
                    tf = tf&((covVec[pks[:-1]]-covVec[vls])>0.05)
                    borders = vls[tf]
                    ff = open(saveDir+'/'+smp+'_chr'+str(ch)+'_borders_'+str(wn*10)+'kb.txt','w')
                    print >> ff,0
                    for x in borders:
                        print >> ff,x
                        # if smp=='AFac':
                            # pl.axvline(10*(x),ls='--',color='black',lw=1)
                    print >> ff,len(covVec)+1
                    ff.close()

                qq+=1
        #     pl.title('Chromosome '+str(ch))
        #     pl.xlabel('Position (kbp)')
        #     pl.ylabel('Relative Coverage Score (a.u.)')
        #     pl.legend(loc=2, bbox_to_anchor=(.8, 0.3),
        #                fancybox=True, shadow=True)
        #     pp.savefig()
        #     pl.close(fig)
        # pp.close()

def getBinMat(intra,chrLens,ch=4):

    binMat = np.zeros((int(np.floor(chrLens.iloc[ch-1,1]/10000)+1),int(np.floor(chrLens.iloc[ch-1,1]/10000)+1)),dtype=int)

    dat = intra[intra.chr1=='chr'+str(ch)]

    for q in range(len(dat)):
        frg1 = dat.frg1.iloc[q]
        frg2 = dat.frg2.iloc[q]
        ix1 = int(np.floor(frg1/10000))
        ix2 = int(np.floor(frg2/10000))

        freq = dat.freq.iloc[q]
        if (ix1<binMat.shape[0]) and (ix2<binMat.shape[0]):
            binMat[ix1,ix2] = freq
    binMat = binMat + binMat.T
    return binMat

def getCoverage(binMat,wn = 5,norm=True):
    covVec1 = [binMat[0:ix,ix:(ix+ix)].mean() for ix in range(1,wn)]
    covVec2 = [binMat[(binMat.shape[0]-ix):ix,ix:binMat.shape[0]].mean() for ix in range((binMat.shape[0]-wn),binMat.shape[0]-1)]
    covVec = [binMat[(ix-wn):ix,ix:(ix+wn)].mean() for ix in range(wn,binMat.shape[0]-wn)]
    covVec = np.r_[covVec1,covVec,covVec2]
    if norm:
        covVec-=min(covVec)
        covVec = covVec.astype(float)
        covVec/=max(covVec)
    return covVec,np.arange(1,binMat.shape[0]-1)

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
