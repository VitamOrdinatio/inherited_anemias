## Use this script to generate a dataframe that relates 199 anemia-enriched loci to 619 different HPO terms.
## Goal is to plot as a single gene-dot in 619th dimensional space, where each axis serves 1 of 619 HPO terms.
## TDA may reveal how these genes are interacting in 619th dimensional space.
## Hypothesis: genes with gene products involved in similar processes (Fanconi anemia translesion repair or Diamond-Blackfan anemia ribosomal ptn subunit functions) should cluster together


## How I generated all 619 HPO listVars
    ## 1) I uploaded a list of 199 loci to gprofiler and resolved discrepancies
    ## 2) I downloaded the gprofiler HPO intersections as CSV
    ## 3) I saved the CSV as XLSX
    ## 4) In Excel, I performed several string operations or manual character replacements to assemble the genPYscript column
    ## 5) The genPYscript column represents python code (that is pasted below)


### VENV pocket for vPCA:
### BASH code to setup VENV pocket for vPCA:
# python3 -m venv vPCA
# source vPCA/bin/activate
# python3 -m pip install --upgrade pip
# python3 -m pip install matplotlib # will also install numpy
# python3 -m pip install seaborn  # will also install pandas
# python3 -m pip install numpy
# python3 -m pip install pandas
# python3 -m pip install -U scikit-learn
# python3 -m pip install kneed
# deactivate



##### BASH
source gotovenv.sh
source vPCA/bin/activate
cd anemia/gProfiler/loci_199/gprofiler_HPO_TDA/prepData
python


##### PYTHON

# LIBRARY IMPORTS
import pandas as pd
import numpy as np
# import seaborn as sns
# import matplotlib.pyplot as plt
# from sklearn.manifold import TSNE

## A list of 199 genes (by geneName)
genes = ['ABCA1','ABCG8','ACVRL1','ADA2','ALG8','ALPL','ANKS6','ATP6V1B1','ATP7B','BAAT','BRCA1','BRCA2','BRIP1','BTK','C3','CD40LG','CD46','CDAN1','CDIN1','CEP164','CEP83','CFB','CFH','CFHR1','CFHR3','CFI','CLCN7','COG1','COL3A1','COL4A1','COL7A1','COQ2','CP','CTC1','DCDC2','DDX41','DKC1','DNAJC21','ELANE','ENG','EPB42','ERCC4','ETV6','F8','F9','FAH','FAM111A','FANCA','FANCB','FANCC','FANCD2','FANCE','FANCF','FANCG','FANCI','FANCL','FARS2','FMO3','FOXP3','G6PC1','G6PC3','GATA1','GBA1','GLA','GLIS2','HAMP','HBA1','HBA2','HBB','HCFC1','HLA-DQA1','HLA-DQB1','IFIH1','INVS','KCNE1','KCNQ1','LIPA','LMBRD1','LPIN2','LYST','MAD2L2','MEN1','MMAA','MMAB','MMACHC','MMADHC','MMP1','MTR','MTRR','MTTP','MUC1','NAF1','NEK8','NPHP1','NPHP3','NPHP4','PALB2','PCCA','PCCB','PEPD','PLEC','PNPO','PRDX1','RAD51','RAD51C','RBM8A','REN','RFWD3','RMRP','RPL11','RPL15','RPL18','RPL26','RPL35','RPL35A','RPL5','RPS10','RPS15A','RPS17','RPS19','RPS24','RPS26','RPS28','RPS29','RPS7','SAMD9','SAMD9L','SBDS','SLC19A2','SLC25A13','SLC2A1','SLC37A4','SLC46A1','SLC4A1','SLC7A7','SLX4','SMAD4','SMARCAL1','SMPD1','SRP54','STK11','SURF1','TBX1','TERC','TERT','TFR2','TGFB1','THBD','TINF2','TMEM67','TREX1','TSR2','UBE2T','UROS','VWF','WAS','WDR19','XK','XPNPEP3','XRCC2','MT-ATP8','MT-ATP6','MT-CO1','MT-CO2','MT-CO3','MT-CYB','MT-ND1','MT-ND2','MT-ND3','MT-ND4L','MT-ND4','MT-ND5','MT-ND6','MT-TA','MT-TR','MT-TN','MT-TD','MT-TC','MT-TE','MT-TQ','MT-TG','MT-TH','MT-TI','MT-TL1','MT-TL2','MT-TK','MT-TM','MT-TF','MT-TP','MT-TS1','MT-TS2','MT-TT','MT-TW','MT-TY','MT-TV','MT-RNR1','MT-RNR2','POLG','POLG2']


## A dictionary that relates all 199 geneNames to all 199 ENSEMBL geneIDs
dictGenes = {'ABCA1': 'ENSG00000165029','ABCG8': 'ENSG00000143921','ACVRL1': 'ENSG00000139567','ADA2': 'ENSG00000093072','ALG8': 'ENSG00000159063','ALPL': 'ENSG00000162551','ANKS6': 'ENSG00000165138','ATP6V1B1': 'ENSG00000116039','ATP7B': 'ENSG00000123191','BAAT': 'ENSG00000136881','BRCA1': 'ENSG00000012048','BRCA2': 'ENSG00000139618','BRIP1': 'ENSG00000136492','BTK': 'ENSG00000010671','C3': 'ENSG00000125730','CD40LG': 'ENSG00000102245','CD46': 'ENSG00000117335','CDAN1': 'ENSG00000140326','CDIN1': 'ENSG00000186073','CEP164': 'ENSG00000110274','CEP83': 'ENSG00000173588','CFB': 'ENSG00000243649','CFH': 'ENSG00000000971','CFHR1': 'ENSG00000244414','CFHR3': 'ENSG00000116785','CFI': 'ENSG00000205403','CLCN7': 'ENSG00000103249','COG1': 'ENSG00000166685','COL3A1': 'ENSG00000168542','COL4A1': 'ENSG00000187498','COL7A1': 'ENSG00000114270','COQ2': 'ENSG00000173085','CP': 'ENSG00000047457','CTC1': 'ENSG00000178971','DCDC2': 'ENSG00000146038','DDX41': 'ENSG00000183258','DKC1': 'ENSG00000130826','DNAJC21': 'ENSG00000168724','ELANE': 'ENSG00000197561','ENG': 'ENSG00000106991','EPB42': 'ENSG00000166947','ERCC4': 'ENSG00000175595','ETV6': 'ENSG00000139083','F8': 'ENSG00000185010','F9': 'ENSG00000101981','FAH': 'ENSG00000103876','FAM111A': 'ENSG00000166801','FANCA': 'ENSG00000187741','FANCB': 'ENSG00000181544','FANCC': 'ENSG00000158169','FANCD2': 'ENSG00000144554','FANCE': 'ENSG00000112039','FANCF': 'ENSG00000183161','FANCG': 'ENSG00000221829','FANCI': 'ENSG00000140525','FANCL': 'ENSG00000115392','FARS2': 'ENSG00000145982','FMO3': 'ENSG00000007933','FOXP3': 'ENSG00000049768','G6PC1': 'ENSG00000131482','G6PC3': 'ENSG00000141349','GATA1': 'ENSG00000102145','GBA1': 'ENSG00000177628','GLA': 'ENSG00000102393','GLIS2': 'ENSG00000126603','HAMP': 'ENSG00000105697','HBA1': 'ENSG00000206172','HBA2': 'ENSG00000188536','HBB': 'ENSG00000244734','HCFC1': 'ENSG00000172534','HLA-DQA1': 'ENSG00000196735','HLA-DQB1': 'ENSG00000179344','IFIH1': 'ENSG00000115267','INVS': 'ENSG00000119509','KCNE1': 'ENSG00000180509','KCNQ1': 'ENSG00000053918','LIPA': 'ENSG00000107798','LMBRD1': 'ENSG00000168216','LPIN2': 'ENSG00000101577','LYST': 'ENSG00000143669','MAD2L2': 'ENSG00000116670','MEN1': 'ENSG00000133895','MMAA': 'ENSG00000151611','MMAB': 'ENSG00000139428','MMACHC': 'ENSG00000132763','MMADHC': 'ENSG00000168288','MMP1': 'ENSG00000196611','MTR': 'ENSG00000116984','MTRR': 'ENSG00000124275','MTTP': 'ENSG00000138823','MUC1': 'ENSG00000185499','NAF1': 'ENSG00000145414','NEK8': 'ENSG00000160602','NPHP1': 'ENSG00000144061','NPHP3': 'ENSG00000113971','NPHP4': 'ENSG00000131697','PALB2': 'ENSG00000083093','PCCA': 'ENSG00000175198','PCCB': 'ENSG00000114054','PEPD': 'ENSG00000124299','PLEC': 'ENSG00000178209','PNPO': 'ENSG00000108439','PRDX1': 'ENSG00000117450','RAD51': 'ENSG00000051180','RAD51C': 'ENSG00000108384','RBM8A': 'ENSG00000265241','REN': 'ENSG00000143839','RFWD3': 'ENSG00000168411','RMRP': 'ENSG00000277027','RPL11': 'ENSG00000142676','RPL15': 'ENSG00000174748','RPL18': 'ENSG00000063177','RPL26': 'ENSG00000161970','RPL35': 'ENSG00000136942','RPL35A': 'ENSG00000182899','RPL5': 'ENSG00000122406','RPS10': 'ENSG00000124614','RPS15A': 'ENSG00000134419','RPS17': 'ENSG00000182774','RPS19': 'ENSG00000105372','RPS24': 'ENSG00000138326','RPS26': 'ENSG00000197728','RPS28': 'ENSG00000233927','RPS29': 'ENSG00000213741','RPS7': 'ENSG00000171863','SAMD9': 'ENSG00000205413','SAMD9L': 'ENSG00000177409','SBDS': 'ENSG00000126524','SLC19A2': 'ENSG00000117479','SLC25A13': 'ENSG00000004864','SLC2A1': 'ENSG00000117394','SLC37A4': 'ENSG00000137700','SLC46A1': 'ENSG00000076351','SLC4A1': 'ENSG00000004939','SLC7A7': 'ENSG00000155465','SLX4': 'ENSG00000188827','SMAD4': 'ENSG00000141646','SMARCAL1': 'ENSG00000138375','SMPD1': 'ENSG00000166311','SRP54': 'ENSG00000100883','STK11': 'ENSG00000118046','SURF1': 'ENSG00000148290','TBX1': 'ENSG00000184058','TERC': 'ENSG00000270141','TERT': 'ENSG00000164362','TFR2': 'ENSG00000106327','TGFB1': 'ENSG00000105329','THBD': 'ENSG00000178726','TINF2': 'ENSG00000092330','TMEM67': 'ENSG00000164953','TREX1': 'ENSG00000213689','TSR2': 'ENSG00000158526','UBE2T': 'ENSG00000077152','UROS': 'ENSG00000188690','VWF': 'ENSG00000110799','WAS': 'ENSG00000015285','WDR19': 'ENSG00000157796','XK': 'ENSG00000047597','XPNPEP3': 'ENSG00000196236','XRCC2': 'ENSG00000196584','MT-ATP8': 'ENSG00000228253','MT-ATP6': 'ENSG00000198899','MT-CO1': 'ENSG00000198804','MT-CO2': 'ENSG00000198712','MT-CO3': 'ENSG00000198938','MT-CYB': 'ENSG00000198727','MT-ND1': 'ENSG00000198888','MT-ND2': 'ENSG00000198763','MT-ND3': 'ENSG00000198840','MT-ND4L': 'ENSG00000212907','MT-ND4': 'ENSG00000198886','MT-ND5': 'ENSG00000198786','MT-ND6': 'ENSG00000198695','MT-TA': 'ENSG00000210127','MT-TR': 'ENSG00000210174','MT-TN': 'ENSG00000210135','MT-TD': 'ENSG00000210154','MT-TC': 'ENSG00000210140','MT-TE': 'ENSG00000210194','MT-TQ': 'ENSG00000210107','MT-TG': 'ENSG00000210164','MT-TH': 'ENSG00000210176','MT-TI': 'ENSG00000210100','MT-TL1': 'ENSG00000209082','MT-TL2': 'ENSG00000210191','MT-TK': 'ENSG00000210156','MT-TM': 'ENSG00000210112','MT-TF': 'ENSG00000210049','MT-TP': 'ENSG00000210196','MT-TS1': 'ENSG00000210151','MT-TS2': 'ENSG00000210184','MT-TT': 'ENSG00000210195','MT-TW': 'ENSG00000210117','MT-TY': 'ENSG00000210144','MT-TV': 'ENSG00000210077','MT-RNR1': 'ENSG00000211459','MT-RNR2': 'ENSG00000210082','POLG': 'ENSG00000140521','POLG2': 'ENSG00000256525'}



# Test dictionary access:
dictGenes['SURF1']
# 'ENSG00000148290'


dfAnemiaHPO = pd.DataFrame()
dfAnemiaHPO = pd.read_csv('199loci_gProfiler_hsapiens_HPO_GO_only.csv')


# Pull out 12 columns.  These represent different facets of the HPO dataset
hpo_source = dfAnemiaHPO['source'].tolist()  
hpo_term_VarName = dfAnemiaHPO['HPO_term_VarName'].tolist()  
hpo_term_fullname = dfAnemiaHPO['HPO_term_name'].tolist()  
hpo_term_id = dfAnemiaHPO['HPO_term_id'].tolist()  
hpo_highlighted = dfAnemiaHPO['highlighted'].tolist()  
hpo_adjusted_p_value = dfAnemiaHPO['adjusted_p_value'].tolist()  
hpo_negative_log10_of_adjusted_p_value = dfAnemiaHPO['negative_log10_of_adjusted_p_value'].tolist()  
hpo_term_size = dfAnemiaHPO['term_size'].tolist()  
hpo_query_size = dfAnemiaHPO['query_size'].tolist()  
hpo_intersection_size = dfAnemiaHPO['intersection_size'].tolist()  
hpo_effective_domain_size = dfAnemiaHPO['effective_domain_size'].tolist()  
hpo_intersections = dfAnemiaHPO['intersections'].tolist()  


### Demo code to retrieve all intersecting genes of the first HPO element (Anemia):
# x = hpo_intersections[0]
# x is a string variable that characters which represents a gene list 'ABCA1,ABCG8,ACVRL1,ADA2,ALG8,ALPL,ATP7B,BRCA1,BRCA2,BRIP1,BTK,C3,CD40LG,CD46,CDAN1,CDIN1,CFB,CFH,CFHR1,CFHR3,CFI,CLCN7,COG1,COL3A1,COL4A1,COL7A1,COQ2,CP,CTC1,DCDC2,DDX41,DKC1,DNAJC21,ELANE,ENG,EPB42,ERCC4,ETV6,F8,FAH,FAM111A,FANCA,FANCB,FANCC,FANCD2,FANCE,FANCF,FANCG,FANCI,FANCL,FARS2,FMO3,FOXP3,G6PC3,GATA1,GBA1,GLA,HAMP,HBA1,HBA2,HBB,HLA-DQB1,IFIH1,KCNE1,KCNQ1,LIPA,LMBRD1,LPIN2,LYST,MAD2L2,MMAA,MMAB,MMACHC,MMADHC,MMP1,MTR,MTRR,MTTP,MUC1,NAF1,NPHP1,NPHP4,PALB2,PCCA,PCCB,PEPD,PLEC,PNPO,PRDX1,RAD51,RAD51C,RBM8A,REN,RFWD3,RMRP,RPL11,RPL15,RPL18,RPL26,RPL35,RPL35A,RPL5,RPS10,RPS15A,RPS17,RPS19,RPS24,RPS26,RPS28,RPS29,RPS7,SAMD9,SAMD9L,SBDS,SLC19A2,SLC25A13,SLC2A1,SLC37A4,SLC46A1,SLC4A1,SLC7A7,SLX4,SMAD4,SMARCAL1,SMPD1,SRP54,STK11,SURF1,TBX1,TERT,TFR2,TGFB1,THBD,TINF2,TMEM67,TREX1,TSR2,UBE2T,UROS,WAS,XRCC2'

## Demo code on how to search using string class find methods:
# x.find('ABCA1') # 0    which indicates that the ABCA1 gene starts at index 0
# x.find('abca1') # -1   which indicates that the query is not found
# x.find('TMEM67') # 774 which indicates that the TMEM67 gene starts at index 774
# x.find('tmem67') # -1 which indicates that the query is not found


## PSEUDOCODE
    # 1) Fetch an HPO (1 of 539), track using counter i 
    # 2) Store this fetched HPO varName in a temp string variable.  Need this for writing a pd column.
    # 3) Store intersections of current fetched HPO in a string variable
    # 4) Flush an accumulator listVar, iterated by j for genes.  Increment j counter when moving over gene list.
    # 5) The accumulator listVar will store either 0 or 1 values, depending on gene intersection status
    # 6) Write the accumulator listVar to a pandas dataframe column
    # 7) Increment i counter, and work on next HPO


# Initalize i counter, which tracks hpo items:
i = 0

newdf = pd.DataFrame()
newdf['gene'] = genes

for i in range(len(hpo_term_VarName)):
    # flush accumulator listVar
    geneHits = []
    # get current hpo term
    hpoColName = hpo_term_VarName[i]
    hpoHits = hpo_intersections[i]
    # initialize j counter which tracks genes:
    j = 0
    for j in range(len(genes)):
        # get a gene
        gene = genes[j]
        query = hpoHits.find(gene)
        # Test query variable.  If -1, then gene wsa not found on HPO list.
        if query == -1:
            geneHits.append(0)
        elif query >= 0:
            geneHits.append(1)
    newdf[hpoColName] = geneHits
    

newdf.to_csv('anemia_199_loci_HPO_TDAprep.csv')