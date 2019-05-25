#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 10:29:11 2018

@author: yahya
"""

import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt;
import seaborn as sns;
import pandas as pd;
import numpy as np;
import pickle
import pandas as pd;
from pandas import DataFrame;
import scipy.stats;
import scipy;
from collections import defaultdict;
import statsmodels.stats.multitest;
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

'''

Fig2f: n=600(300/300)
Fig2f: Atopobium_vaginae FDR-adjusted p-value:0.0062832683299410145
Fig2f: Clostridiales_BVAB2 FDR-adjusted p-value:0.04879595353645941
Fig2f: Dialister_cluster51 FDR-adjusted p-value:0.04879595353645941
Fig2f: Dialister_micraerophilus FDR-adjusted p-value:0.011643523563466218
Fig2f: Gardnerella_vaginalis FDR-adjusted p-value:0.0016785634248519816
Fig2f: Lactobacillus_crispatus_cluster FDR-adjusted p-value:0.03080785167133271
Fig2f: Lactobacillus_iners FDR-adjusted p-value:0.00893978100322039
Fig2f: Prevotella_bivia FDR-adjusted p-value:0.0023654296305745193
Fig2f: Prevotella_cluster2 FDR-adjusted p-value:0.0023654296305745193
Fig2f: Sneathia_amnii FDR-adjusted p-value:0.005173991060656701
Fig2f: Sneathia_sanguinegens FDR-adjusted p-value:0.03385959206111903

Fig2g: n=312(156/156)
Fig2g: Atopobium_vaginae FDR-adjusted p-value:0.010795962736652393
Fig2g: Clostridiales_BVAB2
Fig2g: Dialister_cluster51
Fig2g: Dialister_micraerophilus
Fig2g: Gardnerella_vaginalis FDR-adjusted p-value:0.029979130084145884
Fig2g: Lactobacillus_crispatus_cluster
Fig2g: Lactobacillus_iners
Fig2g: Prevotella_bivia FDR-adjusted p-value:0.011046625028185025
Fig2g: Prevotella_cluster2 FDR-adjusted p-value:0.029979130084145884
Fig2g: Sneathia_amnii FDR-adjusted p-value:0.0207524729707748
Fig2g: Sneathia_sanguinegens

Fig2h: n=122(61/61)
Fig2h: Atopobium_vaginae
Fig2h: Clostridiales_BVAB2
Fig2h: Dialister_cluster51
Fig2h: Dialister_micraerophilus
Fig2h: Gardnerella_vaginalis
Fig2h: Lactobacillus_crispatus_cluster FDR-adjusted p-value:0.03513569801232595
Fig2h: Lactobacillus_iners
Fig2h: Prevotella_bivia
Fig2h: Prevotella_cluster2
Fig2h: Sneathia_amnii
Fig2h: Sneathia_sanguinegens

Fig2i: n=166(83/83)
Fig2i: Atopobium_vaginae
Fig2i: Clostridiales_BVAB2
Fig2i: Dialister_cluster51
Fig2i: Dialister_micraerophilus
Fig2i: Gardnerella_vaginalis FDR-adjusted p-value:0.0017991165737597527
Fig2i: Lactobacillus_crispatus_cluster
Fig2i: Lactobacillus_iners
Fig2i: Prevotella_bivia
Fig2i: Prevotella_cluster2 FDR-adjusted p-value:0.0388409186211142
Fig2i: Sneathia_amnii FDR-adjusted p-value:0.03569263943844617
Fig2i: Sneathia_sanguinegens

'''

subfigure='g';
dfAll = pickle.load(open('metadata.pickle','rb'))
df_data = pickle.load(open('stirrups_profiles.pickle','rb'))
(df1,df2,df3,df4,df5,df2s,df3s,df4s,df5s,df5v)=dfAll

if (subfigure=='g'):
    df1.drop(df1[df1.ethnicity_preg != 'african_american'].index, inplace=True)
    myDir = 'g_AA_'
if (subfigure=='h'):
    df1.drop(df1[df1.ethnicity_preg != 'caucasian'].index, inplace=True)
    myDir = 'h_C_'
if (subfigure=='i'):
    df1.drop(df1[df1.ethnicity_preg != 'hispanic_or_latino'].index, inplace=True)
    myDir = 'i_HL_'
if (subfigure=='f'):
    myDir = 'f_ALL_'

if (subfigure!='f'):

    goodTaxaALL=['Atopobium_vaginae', 
    'Clostridiales_BVAB2', 
    'Dialister_cluster51', 
    'Dialister_micraerophilus', 
    'Gardnerella_vaginalis', 
    'Lactobacillus_crispatus_cluster', 
    'Lactobacillus_iners', 
    'Prevotella_bivia', 
    'Prevotella_cluster2', 
    'Sneathia_amnii', 
    'Sneathia_sanguinegens'];
    
    
    goodTaxaALL_VT=['VT-Atopobium_vaginae', 
    'VT-Clostridiales_BVAB2', 
    'VT-Dialister_cluster51', 
    'VT-Dialister_micraerophilus', 
    'VT-Gardnerella_vaginalis', 
    'VT-Lactobacillus_crispatus_cluster', 
    'VT-Lactobacillus_iners', 
    'VT-Prevotella_bivia', 
    'VT-Prevotella_cluster2', 
    'VT-Sneathia_amnii', 
    'VT-Sneathia_sanguinegens'];
else:
    goodTaxaALL=[];
    goodTaxaALL_VT=[];    

f = open ('Abbrev-Species.csv','r')

#Create dictionary for abbreviations

abbrev = {}

for line in iter(f):

    ln = line.split(",")

    abbrev[ln[0].strip()] = ln[1].strip()

f.close()

def boxPlotAll(x_all,isSignif,y,colNames,pAdj,log10='log10-4',L='N',R='P'):
    doYLabel=False;
    if (subfigure=='f'):
        doYLabel=True;
        
    letters='abcdefghijklmnop';
    
    yp=y>0;
    yn=y<=0;

    fCnt=isSignif.shape[0];
    
    sigCnt=0;
    for cI in range(0,fCnt):
        txName=colNames[cI][3:];
        if (isSignif[cI]>0 or txName in goodTaxaALL):
            sigCnt=sigCnt+1;
            
    if (sigCnt<1):
        return None;
    sigCntX=sigCnt;
    if (sigCntX<2):
        sigCntX=2;

    f, ax = plt.subplots(nrows=1,ncols=sigCntX, sharex=False, sharey=False)
    
    #color of median
    medianprops = dict(color='white') 
    whiskerprops=dict(color='black',linestyle='-')
    boxprops=dict(color='black')
    
    sI=0;
    
    bplots=[];
    
    print("Fig2%s"%subfigure)
    nn=x_all.shape[0];
    nP=np.count_nonzero(yp);
    nN=np.count_nonzero(yn);
    print("Blue: %s, Red: %s"%(L,R))
    print("On the plot, from left to right:")
    print("Fig2"+subfigure+": "+"n="+str(nn)+"("+str(nP)+"/"+str(nN)+")");
    for cI in range(0,fCnt):
        txName=colNames[cI][3:];
        #print(txName)
        if (isSignif[cI]>0 or (txName in goodTaxaALL)):
            
            if (isSignif[cI]>0):
                print("Fig2"+subfigure+": "+colNames[cI][3:]+" FDR-adjusted p-value:"+str(pAdj[cI]))
            else:
                print("Fig2"+subfigure+": "+colNames[cI][3:])
        
            x=x_all[:,cI];

            xp=x[yp];
            xn=x[yn];
            

            nCnt=xn.shape[0];
            pCnt=xp.shape[0];
            
            axC=ax[sI];
            
            medianprops = dict(color='white') 
            whiskerprops=dict(color='black',linestyle='-')
            boxprops=dict(color='black');
            
            yl=xn;
            yr=xp;
            bplt=axC.boxplot((yl,yr), whis='range',labels=[L,R], sym='',manage_xticks=True,positions=[-0.1,0.6],widths=[0.35,0.35],medianprops=medianprops,notch=False,patch_artist=True,whiskerprops=whiskerprops,boxprops=boxprops)  #whis=[25, 75]


    
            if (isSignif[cI]>0):
                bplt['boxes'][0].set_facecolor('#4393c3')
                bplt['boxes'][1].set_facecolor('#d6604d')
            else:
                bplt['boxes'][0].set_facecolor('#999999')
                bplt['boxes'][1].set_facecolor('#999999')
            if (sI>0):
                axC.spines['left'].set_visible(False)
            axC.spines['right'].set_visible(False)
            axC.spines['top'].set_visible(False)
            axC.spines['bottom'].set_visible(False)
            axC.tick_params(labeltop='off', labelright='off',top='off', right='off')
            if (sI==0):
                if (log10=='log10-4'):
                    axC.set_ylim(ymin=-0.1,ymax=4.1)
                    axC.set_yticks([0,1,2,3,4])
                    axC.set_yticklabels(['0',"$10^{-3}$", "$10^{-2}$", "$10^{-1}$", "1"])
            else:
                if (log10=='log10-4'):
                    axC.set_ylim(ymin=-0.1,ymax=4.1)
                    axC.set_yticks(())

            myabb=abbrev[colNames[cI][3:]];
            axC.set_xlabel(myabb, rotation=45) 
            sI=sI+1;
#        elif (subfigure!='f'):
#            print("NOPLOT: "+colNames[cI][3:])
            

    if (doYLabel):
        ax[0].set_ylabel('vaginal 16S rRNA abundance')
    
    plt.show();
        
    return f;




def getMWUPs(x_all,x_all2,y,log10=False):
    fCnt=x_all.shape[1];
    fCnt2=x_all2.shape[1];

    yp=y>0;
    yn=y<=0;

    cols=list(x_all.columns);
    cols2=list(x_all2.columns);
    

    print((len(cols),len(cols2)))

    xp=x_all[yp];
    xn=x_all[yn];
    p=np.ones((fCnt,))
    fStats=np.zeros((fCnt,100));

    for i in range(0,fCnt):
        try:    
            xpf=xp.iloc[:,i];
            xnf=xn.iloc[:,i];
            _,p[i]=scipy.stats.mannwhitneyu(xpf,xnf,alternative='two-sided');
            fStats[i,0]=np.mean(xp)
        except ValueError:
            pass

    hAdj, pAdj, alphacSidak, alphacBonf = statsmodels.stats.multitest.multipletests(p,alpha=0.05,method='fdr_bh')
    
    p2=np.ones((fCnt2,));
    pAdj2=np.ones((fCnt2,))
    
    for i in range(0,fCnt):
        co=cols[i];
        co2=cols2.index(co);
        p2[co2]=p[i];
        pAdj2[co2]=pAdj[i];
        
    
    return p, pAdj, p2, pAdj2, alphacBonf;


def normalizeX(x,nt):
    x=x.copy();
    if (nt=='log10-4'):
        x=x-0.0001;
        x[x<0]=0;
        return np.log10((x+0.0001)/0.0001);


def makePlots(df,df2,y,commonIndex, columns,expName,normType2,L,R):
    x_all_base=df.as_matrix().astype(dtype='float32');
    x_all2=normalizeX(x_all_base,normType2);

    x_all_baseALL=df2.as_matrix().astype(dtype='float32');
    x_all2ALL=normalizeX(x_all_baseALL,normType2);

    ### Mann–Whitney 
    pRaw, pAdj ,pRaw2, pAdj2 , alphacBonf = getMWUPs(df,df2,y.flatten(),log10=False);
    fCnt=pAdj2.shape[0];

    isSignif=np.zeros((fCnt,))
    
    for cI in range(0,fCnt):
        if (pAdj2[cI]<=0.05):
            isSignif[cI]=1;

    boxPlotFileName='paperFigs/Fig2%s.pdf'%subfigure;

    f=boxPlotAll(x_all2ALL,isSignif,y.flatten(),list(df2.columns),pAdj2,log10=normType2,L=L,R=R)
    if (f):
        f.savefig(boxPlotFileName,dpi=600,bbox_inches='tight')
    
    correctedP_FDR = pAdj[pAdj<0.05]
    correctedP_Bonf = pRaw[pRaw<alphacBonf]
    

    featuresOfInterestFDR =  df.columns[pAdj<0.05]
    featuresOfInterestBONF =  df.columns[pRaw<alphacBonf]
    
    resultDf = pd.DataFrame(None, index=commonIndex, columns=columns )
    resultDf.loc[featuresOfInterestFDR, columns[0]] = correctedP_FDR
    resultDf.loc[featuresOfInterestBONF, columns[1]] = correctedP_Bonf
   
    return resultDf
    
def paperFilter(df,subfigure):
    ## Criteria 1
    
    subX = df.copy()>=0.01
    subX_5Percent = (subX.copy().sum(0)/subX.shape[0])>0.05
    subX_5PercentIndex = subX_5Percent[subX_5Percent==True]

    
    # Criteria 2
    subX = df.copy()>=0.001
    subX_15Percent = (subX.sum(0)/subX.shape[0])>0.15
    subX_15PercentIndex = subX_15Percent[subX_15Percent==True]
    
    ## merging the two criteria
    ind1= subX_5PercentIndex.index
    ind2 = subX_15PercentIndex.index
    significant = ind1.union(ind2)
    df_threshod = df.copy().loc[:, significant]  #new

    if (subfigure!='f'):
        significant2=significant.union(goodTaxaALL_VT)
        df_threshod2 = df.copy().loc[:, significant2]  #new
    else:
        df_threshod2 = df.copy().loc[:, significant]  #new
   
    return df_threshod,df_threshod2

def runPipeline(df_filterd,df_filterd2, y, name,L,R):
    df_FilteredThreshod0001=df_filterd.copy().fillna(0.0);
    df_FilteredThreshod0001[df_FilteredThreshod0001<0.0001]=0.0
    nameN=name+str(df_FilteredThreshod0001.shape[0]);

    df_FilteredThreshod0001_2=df_filterd2.copy().fillna(0.0);
    df_FilteredThreshod0001_2[df_FilteredThreshod0001_2<0.0001]=0.0


    resultsFT0001df = makePlots(df_FilteredThreshod0001,df_FilteredThreshod0001_2,y,UnifiedIndex, columns=['pAdjusted_F_0001', 'pBonferroni_F_0001'],expName=nameN+"_thresh_0001",normType2='log10-4',L=L,R=R)
    
###########


######################################################################################################
# Pregnant vs nonPregnant
######################################################################################################
### Extracting the classes  and the data
class1Preg = df1['sid_preg']  # class 1 in sid
class0NotPreg = df1['sid_notpreg']  # class 0 in sid

mydict = pd.Series(df5.vid.values,index=df5.sid).to_dict() # sid vid dictionary 
class1Preg = class1Preg.copy().apply(lambda x: mydict[x])  #converting sid to vid c1
class0NotPreg = class0NotPreg.copy().apply(lambda x: mydict[x])  #converting sid to vid c0
y1 = df_data.index.isin(class1Preg)  #extracting indices of c1 from the big df
y2 = df_data.index.isin(class0NotPreg)  #extracting indices of c0 from the big df
y_temp = y1+y2  # merging the indices of c1 and c0 in the big dataframe

data_2groups_P_N = df_data.copy()[y_temp]
y = data_2groups_P_N.index.isin(class1Preg)  # class 1 is pregnant class 0 is not pregnent 
######################################################################################################
# all data without Fitering or Threshold
######################################################################################################
class1Preg = df1['sid_preg']  # class 1 in sid
class0NotPreg = df1['sid_notpreg']  # class 0 in sid

mydict_VtoS = pd.Series(df5.sid.values,index=df5.vid).to_dict() # vid sid dictionary 
mydict_StoV = pd.Series(df5.vid.values,index=df5.sid).to_dict() # sid vid dictionary 

class1Preg = class1Preg.copy().apply(lambda x: mydict_StoV[x])  #converting sid to vid c1
class0NotPreg = class0NotPreg.copy().apply(lambda x: mydict_StoV[x])  #converting sid to vid c0
y1 = df_data.index.isin(class1Preg)  #extracting indices of c1 from the big df
y2 = df_data.index.isin(class0NotPreg)  #extracting indices of c0 from the big df
y_temp = y1+y2  # merging the indices of c1 and c0 in the big dataframe

data_2groups = df_data.copy()[y_temp]
y = data_2groups.index.isin(class1Preg)  # class 1 is pregnant class 0 is not pregnent 


#########################################################################################################
### all data analysis
########################################################################################################
##results_noT_noF= levelOfanalysis(data_2groups,y,columns=['pAdjusted', 'pBonferroni'] )
#
#
#######################################################################################################
##  data with deleting features that is below paper criteria
#######################################################################################################
data_2groups_P_N_Filtered,data_2groups_P_N_Filtered2 = paperFilter(data_2groups_P_N,subfigure) 


AllFeatures =data_2groups_P_N_Filtered.columns

meanOfFilterData = data_2groups[AllFeatures].describe().loc['mean',:]  # AA has more features after filtering
meanOfFilterData.sort_values(ascending=False,inplace=True)
UnifiedIndex = meanOfFilterData.index
#######################################################################################################
#### filter the thresholded data by set small quantity of reads 0.0001 to 0   
#######################################################################################################\


runPipeline(data_2groups_P_N_Filtered, data_2groups_P_N_Filtered2, y, myDir+'P_N_n=','N','P')

