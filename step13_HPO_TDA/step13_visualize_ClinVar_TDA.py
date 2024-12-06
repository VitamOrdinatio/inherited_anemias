## The goal of this script is to to generate SVG circles that are of relative proportions for the purposes of TDA node visualization.

## In its current implementation, the script simply performs a log2-transformation of the geneNumber size of each TDA node.
## These log2(geneNumber) values are then manually shunted to BioRender to draw approximate circles with relative size differences.

## Additionally, this script outputs a CSV file that serves as a digestion of tracks, nodes, and gene relations for HPO/TDA kmapper returns.


## In future implementations, I'd like to implement the entire graphical output in the following pipeline:    
## Pseudocode: 
    # 1) Set the node number (taken from the kmapper HTML output)
    # 2) Create a dataframe that tracks the number of genes per nodeID
    # 3) Calculate the log2-transformation of that the numGenes per nodeID
    # 4) Plot a series of circles with diameters set to the nodeID's log2(numGenes)  (large nodes have larger diameters)
    # 5) Save as SVG file to disk (assemble later via BioRender)
    
#   For each gene, raw number of each of 4 ClinVar categories: B, LB, LP, P.  Follow this color schema:
        #   CATEGORY: Benign                (dark blue)     #29386F
        #   CATEGORY: Likely benign         (light blue)    #DFEDFA
        #   CATEGORY: Likely pathogenic     (light red)     #FCE5EA
        #   CATEGORY: Pathogenic            (dark red)      #802A2A      


### Here is an ASCII abstraction of the TDA output of the 199 loci + 619 HPOterms:
    ## A node is a cluster.  Nodes/clusters have guild members.  
    ## A guild member is a gene.
    ## TDA components often resemble subway lines.  A TDA component is known as a track.


    
    ## In our analysis, we have 3 distinct TDA tracks (components) with varying numbers of nodes (clusters) and each node contains a varying number of guildmembers (genes).
    
    ## Track 1     track::structure                       track::nodeIDs                track::geneNUM
    ##             0-0-0-0-0-0                            1-2-3-4-5-6                   42-64-54-31-3-6
    ##                   |                                      |                                |
    ##                   0                                      7                                6

    ## Track 2
    ##             0-0-0                                  1-2-3                         13-17-3

    ## Track 3
    ##             0-0                                    1-2                           13-11




### 1SHOT: BASH code to setup VENV pocket for vPCA:
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

### After VENV vPCA construction, upon WSL bootup:
source gotovenv.sh
source vPCA/bin/activate
cd /mnt/c/wslshare/github/inherited_anemias/step13_HPO_TDA
python

############################################################ IMPORTS
import pandas as pd
import numpy as np
import math
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import math


## In our analysis, we have 3 distinct TDA tracks (components) with varying numbers of nodes (clusters) and each node contains a varying number of guildmembers (genes).

## Track 1     track::structure                       track::nodeIDs                track::geneNUM
##             0-0-0-0-0-0                            1-2-3-4-5-6                   42-64-54-31-3-6
##                   |                                      |                                |
##                   0                                      7                                6

## Track 2
##             0-0-0                                  1-2-3                         13-17-3

## Track 3
##             0-0                                    1-2                           13-11


# Set the total number of TDA components (which I'm calling tracks akin to subway lines).
numTracks = 3

# Set the total node number (not index number) for each track.  Get this by manually opening the kmapper HTML output, and counting the nodes.
track1_numNodes = 7
track2_numNodes = 3
track3_numNodes = 2

# Set each nodes' number of genes (guild membership).  Get this by manually examining the kmapper HTML output, and recording the guild membership for each node.  Guilds are TDA clusters, and guild members are genes.
numGenes_PerNode_Track1 = [42,64,54,31,3,6,6]
numGenes_PerNode_Track2 = [13,17,3]
numGenes_PerNode_Track3 = [13,11]

# Based on above information, manually assert listVariables for dataframe iteration. Calculate log2-tranformation of geneNumList
trackList = [1,1,1,1,1,1,1,2,2,2,3,3]
nodeList = [1,2,3,4,5,6,7,1,2,3,1,2]
geneNumList = [42,64,54,31,3,6,6, 13,17,3, 13,11]

# Convert trackList to a listVar of string elements
strTrackList = []
i = 0
for i in range(len(trackList)):
    trackLabel = 't' + str(trackList[i])
    strTrackList.append(trackLabel)

# View strTrackList
strTrackList    



# Convert nodeList to a listVar of string elements    
strNodeList = []
i = 0
for i in range(len(nodeList)):
    nodeLabel = 'n' + str(nodeList[i])
    strNodeList.append(nodeLabel)

# View listVar
strNodeList



# Calculate each nodes' log2-transformed geneNumber:
i = 0
log2numGenesPerNode = []
for i in range(len(geneNumList)):
    # Fetch a node's gene number:
    geneNUM = geneNumList[i]
    result = math.log2(geneNUM)
    log2numGenesPerNode.append(result)

# View log2 transformed geneNumber per node:
log2numGenesPerNode
# [6.741466986401147, 5.882643049361842, 5.169925001442312, 4.643856189774724, 3.807354922057604, 2.0]


## Write CSV to disk
newdf = pd.DataFrame()
newdf['trackID'] = strTrackList
newdf['nodeID'] = strNodeList
newdf['numGenes'] = geneNumList
newdf['log2(numGenes)'] = log2numGenesPerNode
newdf.to_csv('./out/anemia_199loci_HPO_TDA_3tracks_nodeSizes.csv')





#### Next, we need to manually annotate the keplermapper output file: 'Anemia199HPOBands.csv'
#### MANUAL CHANGES (all in EXCEL):
    # 1) Save AS XLSX file
    # 2) Add column names: cluster, band, gene
    # 3) Appended two new columns before 'cluster'
    #       a) 'track'
    #       b) 'node'
    # 4) Appended two new columns after 'gene'
    #       a) 'disorder'
    #       b) 'alias'
    
##### Manual data entry for 'track' and 'node' columns.  Consulting HTML file as a reference for all 12 nodes across 3 tracks.

    # 5) Sorted sheet (A->Z) by track, then node, then cluster, then band.
    #       a) I have n = 263 genes (not unique) across 3 tracks comprising 12 TDA cluster nodes.
    # 6) Created a new EXCEL tab called listVar
##### Using XLSX saved settings to reassert a pythonic CSV file:
    
hpoTDAtrackList = ['1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3']

hpoTDAnodeList = ['1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','1','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','3','4','4','4','4','4','4','4','4','4','4','4','4','4','4','4','4','4','4','4','4','4','4','4','4','4','4','4','4','4','4','4','5','5','5','6','6','6','6','6','6','7','7','7','7','7','7','1','1','1','1','1','1','1','1','1','1','1','1','1','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','3','3','3','1','1','1','1','1','1','1','1','1','1','1','1','1','2','2','2','2','2','2','2','2','2','2','2']

hpoTDAclusterList = ['Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 0','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 1','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 2','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 3','Cluster 7','Cluster 7','Cluster 7','Cluster 10','Cluster 10','Cluster 10','Cluster 10','Cluster 10','Cluster 10','Cluster 5','Cluster 5','Cluster 5','Cluster 5','Cluster 5','Cluster 5','Cluster 4','Cluster 4','Cluster 4','Cluster 4','Cluster 4','Cluster 4','Cluster 4','Cluster 4','Cluster 4','Cluster 4','Cluster 4','Cluster 4','Cluster 4','Cluster 6','Cluster 6','Cluster 6','Cluster 6','Cluster 6','Cluster 6','Cluster 6','Cluster 6','Cluster 6','Cluster 6','Cluster 6','Cluster 6','Cluster 6','Cluster 6','Cluster 6','Cluster 6','Cluster 6','Cluster 8','Cluster 8','Cluster 8','Cluster 9','Cluster 9','Cluster 9','Cluster 9','Cluster 9','Cluster 9','Cluster 9','Cluster 9','Cluster 9','Cluster 9','Cluster 9','Cluster 9','Cluster 9','Cluster 11','Cluster 11','Cluster 11','Cluster 11','Cluster 11','Cluster 11','Cluster 11','Cluster 11','Cluster 11','Cluster 11','Cluster 11']

hpoTDAbandList = ['Band 0','Band 1','Band 10','Band 11','Band 12','Band 13','Band 14','Band 15','Band 16','Band 17','Band 18','Band 19','Band 2','Band 20','Band 21','Band 22','Band 23','Band 24','Band 25','Band 26','Band 27','Band 28','Band 29','Band 3','Band 30','Band 31','Band 32','Band 33','Band 34','Band 35','Band 36','Band 37','Band 38','Band 39','Band 4','Band 40','Band 41','Band 5','Band 6','Band 7','Band 8','Band 9','Band 0','Band 1','Band 10','Band 11','Band 12','Band 13','Band 14','Band 15','Band 16','Band 17','Band 18','Band 19','Band 2','Band 20','Band 21','Band 22','Band 23','Band 24','Band 25','Band 26','Band 27','Band 28','Band 29','Band 3','Band 30','Band 31','Band 32','Band 33','Band 34','Band 35','Band 36','Band 37','Band 38','Band 39','Band 4','Band 40','Band 41','Band 42','Band 43','Band 44','Band 45','Band 46','Band 47','Band 48','Band 49','Band 5','Band 50','Band 51','Band 52','Band 53','Band 54','Band 55','Band 56','Band 57','Band 58','Band 59','Band 6','Band 60','Band 61','Band 62','Band 63','Band 7','Band 8','Band 9','Band 0','Band 1','Band 10','Band 11','Band 12','Band 13','Band 14','Band 15','Band 16','Band 17','Band 18','Band 19','Band 2','Band 20','Band 21','Band 22','Band 23','Band 24','Band 25','Band 26','Band 27','Band 28','Band 29','Band 3','Band 30','Band 31','Band 32','Band 33','Band 34','Band 35','Band 36','Band 37','Band 38','Band 39','Band 4','Band 40','Band 41','Band 42','Band 43','Band 44','Band 45','Band 46','Band 47','Band 48','Band 49','Band 5','Band 50','Band 51','Band 52','Band 53','Band 6','Band 7','Band 8','Band 9','Band 0','Band 1','Band 10','Band 11','Band 12','Band 13','Band 14','Band 15','Band 16','Band 17','Band 18','Band 19','Band 2','Band 20','Band 21','Band 22','Band 23','Band 24','Band 25','Band 26','Band 27','Band 28','Band 29','Band 3','Band 30','Band 4','Band 5','Band 6','Band 7','Band 8','Band 9','Band 0','Band 1','Band 2','Band 0','Band 1','Band 2','Band 3','Band 4','Band 5','Band 0','Band 1','Band 2','Band 3','Band 4','Band 5','Band 0','Band 1','Band 10','Band 11','Band 12','Band 2','Band 3','Band 4','Band 5','Band 6','Band 7','Band 8','Band 9','Band 0','Band 1','Band 10','Band 11','Band 12','Band 13','Band 14','Band 15','Band 16','Band 2','Band 3','Band 4','Band 5','Band 6','Band 7','Band 8','Band 9','Band 0','Band 1','Band 2','Band 0','Band 1','Band 10','Band 11','Band 12','Band 2','Band 3','Band 4','Band 5','Band 6','Band 7','Band 8','Band 9','Band 0','Band 1','Band 10','Band 2','Band 3','Band 4','Band 5','Band 6','Band 7','Band 8','Band 9']


hpoTDAgeneList = ['ABCA1','ANKS6','GLIS2','HLA-DQA1','KCNE1','NAF1','PNPO','TERC','VWF','XK','MT-TA','MT-TR','ATP6V1B1','MT-TN','MT-TD','MT-TC','MT-TE','MT-TQ','MT-TG','MT-TH','MT-TI','MT-TL1','MT-TL2','BAAT','MT-TK','MT-TM','MT-TF','MT-TP','MT-TS1','MT-TS2','MT-TT','MT-TW','MT-TY','MT-TV','CDIN1','MT-RNR1','MT-RNR2','CP','DDX41','F9','FARS2','FMO3','ABCA1','ABCG8','CDAN1','CDIN1','CEP164','CEP83','CFB','COQ2','CP','DCDC2','DDX41','ELANE','ACVRL1','EPB42','ETV6','F8','FAM111A','FARS2','FMO3','FOXP3','G6PC1','G6PC3','GLA','ALG8','HAMP','HBA1','HBA2','HCFC1','INVS','KCNQ1','LIPA','LPIN2','MMAA','MMAB','ALPL','MMP1','MTTP','MUC1','NAF1','NEK8','NPHP4','PCCA','PCCB','PEPD','PNPO','ANKS6','REN','SAMD9L','SLC19A2','SLC2A1','SLC46A1','SMPD1','SURF1','TFR2','THBD','VWF','BAAT','WAS','XPNPEP3','MT-ATP8','MT-ND4L','BTK','CD40LG','CD46','ACVRL1','ALG8','COL4A1','COL7A1','CTC1','DCDC2','DKC1','DNAJC21','ENG','FAH','FOXP3','G6PC3','ATP7B','GLA','HBA1','HLA-DQB1','LYST','MEN1','MMAA','MMAB','MMADHC','MMP1','MTR','CD46','MTRR','NPHP3','PCCA','PCCB','PEPD','PLEC','PRDX1','RMRP','SAMD9','SBDS','CFHR1','SLC25A13','SLC2A1','SLC4A1','SLC7A7','SMARCAL1','SRP54','STK11','TGFB1','TINF2','UROS','CFHR3','WAS','WDR19','MT-ND3','POLG2','CFI','CLCN7','COG1','COL3A1','CFH','CFI','MMACHC','MTR','MTRR','NPHP1','PRDX1','RBM8A','RMRP','SAMD9','SLC37A4','SLC7A7','COL3A1','SMARCAL1','TGFB1','TINF2','TMEM67','TREX1','WDR19','MT-ATP6','MT-CO2','MT-CYB','MT-ND2','CTC1','POLG','DKC1','DNAJC21','GBA1','HBB','IFIH1','LMBRD1','MT-ATP6','MT-CO1','MT-CO2','MT-CO1','MT-CO3','MT-ND1','MT-ND4','MT-ND5','MT-ND6','C3','GBA1','SLC37A4','SMAD4','TERT','TREX1','RPL15','RPL18','RPS29','RPS7','TSR2','RPL35','RPL35A','RPS10','RPS15A','RPS17','RPS24','RPS26','RPS28','RPL11','RPL15','RPS19','RPS24','RPS26','RPS28','RPS29','RPS7','TSR2','RPL18','RPL26','RPL35','RPL35A','RPL5','RPS10','RPS15A','RPS17','ADA2','GATA1','RPL26','BRIP1','FANCB','SLX4','UBE2T','XRCC2','FANCF','FANCG','FANCI','FANCL','MAD2L2','RAD51','RAD51C','RFWD3','BRCA1','BRCA2','RAD51C','ERCC4','FANCA','FANCB','FANCC','FANCD2','FANCE','PALB2','RAD51']

hpoTDAcondnList = ['Tangier disease','renal','renal','Celiac disease, susceptibility to, 1','Jervell and Lange-Nielsen syndrome 2','Pulmonary fibrosis and/or bone marrow failure syndrome, telomere-related, 7','B6','Dyskeratosis congenita, autosomal dominant 1','clot','McLeod neuroacanthocytosis syndrome','mtDNA','mtDNA','renal','mtDNA','mtDNA','mtDNA','mtDNA','mtDNA','mtDNA','mtDNA','mtDNA','mtDNA','mtDNA','Bile acid conjugation defect 1','mtDNA','mtDNA','mtDNA','mtDNA','mtDNA','mtDNA','mtDNA','mtDNA','mtDNA','mtDNA','CDA','mtDNA','mtDNA','Fe','DDX41-related hematologic malignancy predisposition syndrome','clot','OXPHOS','Trimethylaminuria','Tangier disease','Sitosterolemia 1','CDA','CDA','renal','renal','Atypical hemolytic-uremic syndrome with B factor anomaly','CoQ10','Fe','renal','DDX41-related hematologic malignancy predisposition syndrome','Neutropenia, severe congenital, 1, autosomal dominant','Telangiectasia, hereditary hemorrhagic, type 2','HS','clot','clot','Autosomal dominant Kenny-Caffey syndrome','OXPHOS','Trimethylaminuria','Insulin-dependent diabetes mellitus secretory diarrhea syndrome','G6P','G6P','Fabry disease','ALG8 congenital disorder of glycosylation','Hemochromatosis type 2B','thalassemia','thalassemia','B12','renal','Jervell and Lange-Nielsen syndrome 1','Cholesteryl ester storage disease','Majeed syndrome','B12','B12','Infantile hypophosphatasia','ECM','Abetalipoproteinaemia','renal','Pulmonary fibrosis and/or bone marrow failure syndrome, telomere-related, 7','renal','renal','Propionic acidemia','Propionic acidemia','Prolidase deficiency','B6','renal','Familial juvenile hyperuricemic nephropathy type 2','Ataxia-pancytopenia syndrome','Megaloblastic anemia, thiamine-responsive, with diabetes mellitus and sensorineural deafness','Childhood onset GLUT1 deficiency syndrome 2','B9','Niemann-Pick disease, type B','OXPHOS','Hemochromatosis type 3','Atypical hemolytic-uremic syndrome with thrombomodulin anomaly','clot','Bile acid conjugation defect 1','clot','renal','mtDNA','mtDNA','X-linked agammaglobulinemia','Hyper-IgM syndrome type 1','Atypical hemolytic-uremic syndrome with MCP/CD46 anomaly','Telangiectasia, hereditary hemorrhagic, type 2','ALG8 congenital disorder of glycosylation','Brain small vessel disease 1 with or without ocular anomalies','ECM','Cerebroretinal microangiopathy with calcifications and cysts 1','renal','Dyskeratosis congenita, X-linked','Shwachman-Diamond syndrome 1','Telangiectasia, hereditary hemorrhagic, type 1','Tyrosinemia type I','Insulin-dependent diabetes mellitus secretory diarrhea syndrome','G6P','Wilson disease','Fabry disease','thalassemia','Celiac disease, susceptibility to, 1','Chediak-Higashi syndrome','Multiple endocrine neoplasia, type 1','B12','B12','B12','ECM','B12','Atypical hemolytic-uremic syndrome with MCP/CD46 anomaly','B12','renal','Propionic acidemia','Propionic acidemia','Prolidase deficiency','ECM','B12','Metaphyseal dysplasia','MIRAGE syndrome','Shwachman-Diamond syndrome 1','Hemolytic uremic syndrome, atypical, susceptibility to, 1','Neonatal intrahepatic cholestasis due to citrin deficiency','Childhood onset GLUT1 deficiency syndrome 2','renal','Lysinuric protein intolerance','Schimke immuno-osseous dysplasia','Shwachman-Diamond syndrome 1','Peutz-Jeghers syndrome','Diaphyseal dysplasia','Dyskeratosis congenita','Cutaneous porphyria','Hemolytic uremic syndrome, atypical, susceptibility to, 1','clot','renal','mtDNA','mtDNA depletion syndrome','Atypical hemolytic-uremic syndrome with I factor anomaly','Osteopetrosis','COG1 congenital disorder of glycosylation','ECM','Hemolytic uremic syndrome, atypical, susceptibility to, 1','Atypical hemolytic-uremic syndrome with I factor anomaly','B12','B12','B12','renal','B12','clot','Metaphyseal dysplasia','MIRAGE syndrome','G6P','Lysinuric protein intolerance','ECM','Schimke immuno-osseous dysplasia','Diaphyseal dysplasia','Dyskeratosis congenita','renal','renal','renal','mtDNA','mtDNA','mtDNA','mtDNA','Cerebroretinal microangiopathy with calcifications and cysts 1','mtDNA depletion syndrome','Dyskeratosis congenita, X-linked','Shwachman-Diamond syndrome 1','Gaucher disease','thalassemia','Aicardi-Goutieres syndrome 7','B12','mtDNA','mtDNA','mtDNA','mtDNA','mtDNA','mtDNA','mtDNA','mtDNA','mtDNA','Atypical hemolytic-uremic syndrome with C3 anomaly','Gaucher disease','G6P','Juvenile polyposis/hereditary hemorrhagic telangiectasia syndrome','Pulmonary fibrosis and/or bone marrow failure, Telomere-related, 1','renal','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','DBA','Vasculitis due to ADA2 deficiency','clot','DBA','FA','FA','FA','FA','FA','FA','FA','FA','FA','FA','FA','FA','FA','FA','FA','FA','FA','FA','FA','FA','FA','FA','FA','FA']


dfHPOtda = pd.DataFrame()
dfHPOtda['track'] = hpoTDAtrackList
dfHPOtda['node'] = hpoTDAnodeList
dfHPOtda['band'] = hpoTDAbandList
dfHPOtda['gene'] = hpoTDAgeneList
dfHPOtda['disorder_tagName'] = hpoTDAcondnList
dfHPOtda.to_csv('./out/anemia_199loci_HPO_TDA_3tracks_geneHOOKS.csv')