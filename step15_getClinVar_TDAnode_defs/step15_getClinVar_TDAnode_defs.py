## The following python script generates TDA node definition plots by averaging allelic categorical frequencies for all genes found in a particular TDA node.

## The script requires the following two files in the input folder (./in)
    # 1) "Anemia_199_genes_AlleleCount_4_ClinVar_classifiers.csv"
                # This file contains the raw allelic counts (n = 112,534) for all 199 loci
                    # Column Names:
                        # gene_name     <- this column has UNIQUE gene names
                        # B
                        # LB
                        # LP
                        # P
            # NOTA BENE: The CSV file must contain the original allele counts (that are not normalized frequencies).
    # 2) "Anemia199ClinVarBands_colNames.csv"
                # This file contains the TDA output where nodes (aka clusters) contain genes
                    # Column Names:
                        # cluster (aka nodes)
                        # gene_name     <- this column does NOT have unique gene names

# Pseudocode:
    # 1) Load both CSV files into two pandas dataframe variables
    # 2) Perform the equivalent of an SQL inner join using pd.merge
            # Merge on the 'gene_name' column
            # Make sure that the TDA-ClinVar df is loaded first in the .merge(args)
            # Perform a LEFT JOIN
                # This retains all columns of the TDA-ClinVar dataframe
                # Next, the ClinVar-allele-counts data frame is appended to the TDA-ClinVar dataframe.
    # 3) Output mean TDA node definitions for all TDA clusters (aka nodes) as barplots in PNG and SVG formats
                


#   CATEGORY: Benign                (dark blue)     #29386F
#   CATEGORY: Likely benign         (light blue)    #DFEDFA
#   CATEGORY: Likely pathogenic     (light red)     #FCE5EA
#   CATEGORY: Pathogenic            (dark red)      #802A2A       

### CUSTOM COLOR PALETTE:
# COLOR	        FILL
# black	        #363636
# grey	        #E5E5E5
# blue	        #29386F
# light blue	#DFEDFA
# red	        #802A2A
# light red	    #FCE5EA
# green	        #44A043
# light green	#E0F4DA


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
cd /mnt/c/wslshare/github/inherited_anemias/step15_getClinVar_TDAnode_defs
python

### Start PYTHON codeblock:

############################################################ IMPORTS
import pandas as pd
import numpy as np
import math
import seaborn as sns
import matplotlib.pyplot as plt


### read csv data for TDA
tdaDF = pd.read_csv('./in/Anemia199ClinVarBands_colNames.csv') 
# View tdaDF
tdaDF
# [245 rows x 3 columns]

### read csv data for raw ClinVar allelic counts
alleleCountsDF = pd.read_csv('./in/Anemia_199_genes_AlleleCount_4_ClinVar_classifiers.csv') 
# View alleleCountsDF
alleleCountsDF
# [199 rows x 5 columns]

## Append a rowSUM to the loaded alleleCountsDF:
# Use pd.sum function to calculate rowSUMS
alleleCountsDF['alleleSUM'] = alleleCountsDF.sum(axis=1, numeric_only=True)
# View alleleCountsDF
alleleCountsDF
# [199 rows x 6 columns]

## Convert raw allele counts into frequencies:
genes = alleleCountsDF['gene_name'].tolist() 
rawBcounts = alleleCountsDF['B'].tolist() 
rawLBcounts = alleleCountsDF['LB'].tolist() 
rawLPcounts = alleleCountsDF['LP'].tolist() 
rawPcounts = alleleCountsDF['P'].tolist()
alleleSUMS = alleleCountsDF['alleleSUM'].tolist()

# All listVars should currently have 199 elements (each element represents 1 gene)
len(genes)
len(rawBcounts)
len(rawLBcounts)
len(rawLPcounts)
len(rawPcounts)
len(alleleSUMS)

####### The following codeblock is where I convert raw alleleCounts into normalized alleleFreqs
#### This is necessary because our TDA approach operates on the alleleFreqs to prevent raw alleleCounts from dominating the locus behavior.  Thus to get the mean TDA cluster (aka node) defintion, I will need to average the clinvar categorical allele frequencies.

len(genes)  # 199 loci
i = 0

## Accumulator listvars to store a list of float elements.
# We first calculate the relative allelic frequencies PER genetic locus (row)
freqBalleles = []         # floats
freqLBalleles = []        # floats
freqLPalleles = []        # floats
freqPalleles = []         # floats

for i in range(len(genes)):
    # get first gene in list
    gene = genes[i]
    # get raw allele counts (numerators)
    rawB = rawBcounts[i]
    rawLB = rawLBcounts[i]
    rawLP = rawLPcounts[i]
    rawP = rawPcounts[i]
    # get denominator (total allele counts across 4 ClinVar categories: B, LB, LP, and P)
    denom = alleleSUMS[i]
    # Append normalized frequences as float elements to accumulator listvars
    freqBalleles.append(rawB/denom)
    freqLBalleles.append(rawLB/denom)
    freqLPalleles.append(rawLP/denom)
    freqPalleles.append(rawP/denom)

alleleFreqsDF = pd.DataFrame()
alleleFreqsDF['gene_name'] = genes
alleleFreqsDF['B_freq'] = freqBalleles
alleleFreqsDF['LB_freq'] = freqLBalleles
alleleFreqsDF['LP_freq'] = freqLPalleles
alleleFreqsDF['P_freq'] = freqPalleles

alleleFreqsDF
# [199 rows x 5 columns]
# Use pd.sum function to calculate rowSUMS
alleleFreqsDF['alleleFreqSUM'] = alleleFreqsDF.sum(axis=1, numeric_only=True)

## View the per-locus normalized allele frequencies (4 ClinVar categories)
alleleFreqsDF
    # gene_name    B_freq   LB_freq   LP_freq    P_freq  alleleFreqSUM
# 0       ABCA1  0.272035  0.631121  0.019587  0.077258            1.0
# 1       ABCG8  0.230114  0.573864  0.068182  0.127841            1.0
# 2      ACVRL1  0.113801  0.231235  0.203390  0.451574            1.0
# 3        ADA2  0.127168  0.433526  0.089595  0.349711            1.0
# 4        ALG8  0.339713  0.411483  0.095694  0.153110            1.0
# ..        ...       ...       ...       ...       ...            ...
# 194     MT-TV  0.294118  0.294118  0.235294  0.176471            1.0
# 195   MT-RNR1  0.411765  0.500000  0.044118  0.044118            1.0
# 196   MT-RNR2  0.000000  0.200000  0.200000  0.600000            1.0
# 197      POLG  0.078125  0.655625  0.112500  0.153750            1.0
# 198     POLG2  0.132275  0.682540  0.037037  0.148148            1.0
# [199 rows x 6 columns]



## Perform a LEFT JOIN (merge two dataframes, hinge is gene_name column)   
mergeDF = pd.merge(tdaDF,  
                   alleleFreqsDF,  
                   on ='gene_name',  
                   how ='left') 
# View mergeDF:
mergeDF
       # cluster     band gene_name    B_freq   LB_freq   LP_freq    P_freq  alleleFreqSUM
# 0    Cluster 0   Band 0     ABCG8  0.230114  0.573864  0.068182  0.127841            1.0
# 1    Cluster 0   Band 1      ADA2  0.127168  0.433526  0.089595  0.349711            1.0
# 2    Cluster 0   Band 2      ALG8  0.339713  0.411483  0.095694  0.153110            1.0
# 3    Cluster 0   Band 3      ALPL  0.086705  0.431599  0.271676  0.210019            1.0
# 4    Cluster 0   Band 4     ANKS6  0.197842  0.597122  0.039568  0.165468            1.0
# ..         ...      ...       ...       ...       ...       ...       ...            ...
# 240  Cluster 4  Band 13   MT-RNR2  0.000000  0.200000  0.200000  0.600000            1.0
# 241  Cluster 5   Band 0        F8  0.042932  0.073298  0.227225  0.656545            1.0
# 242  Cluster 5   Band 1      RMRP  0.042755  0.040380  0.327791  0.589074            1.0
# 243  Cluster 5   Band 2      TERC  0.061728  0.061728  0.246914  0.629630            1.0
# 244  Cluster 5   Band 3      TSR2  0.082418  0.164835  0.049451  0.703297            1.0
#
# [245 rows x 8 columns]

#### For reference, this is what the LEFT JOIN would have looked like if we joined tdaDF directly to the raw alleleCountsDF.  
#### For instance, ABCG8 has 81 Benign of 352 total alleles.  81/352 = 0.230114 = mergeDF.B_freq[0]
       # cluster     band gene_name   B   LB   LP    P  alleleSUM
# 0    Cluster 0   Band 0     ABCG8  81  202   24   45        352
# 1    Cluster 0   Band 1      ADA2  44  150   31  121        346
# 2    Cluster 0   Band 2      ALG8  71   86   20   32        209
# 3    Cluster 0   Band 3      ALPL  90  448  282  218       1038
# 4    Cluster 0   Band 4     ANKS6  55  166   11   46        278
# ..         ...      ...       ...  ..  ...  ...  ...        ...
# 240  Cluster 4  Band 13   MT-RNR2   0    1    1    3          5
# 241  Cluster 5   Band 0        F8  41   70  217  627        955
# 242  Cluster 5   Band 1      RMRP  18   17  138  248        421
# 243  Cluster 5   Band 2      TERC   5    5   20   51         81
# 244  Cluster 5   Band 3      TSR2  15   30    9  128        182

# [245 rows x 8 columns]

# Extract the LEFT JOIN mergeDF columns into listVars
tdaGenes = mergeDF['gene_name'].tolist() 
normBfreqs = mergeDF['B_freq'].tolist() 
normLBfreqs = mergeDF['LB_freq'].tolist() 
normLPfreqs = mergeDF['LP_freq'].tolist() 
normPfreqs = mergeDF['P_freq'].tolist()
clusters = mergeDF['cluster'].tolist()
bands = mergeDF['band'].tolist()
alleleFreqSUMS = mergeDF['alleleFreqSUM'].tolist()

# Verify that all listVar lengths are 245 elements.
len(tdaGenes)
len(normBfreqs)
len(normLBfreqs)
len(normLPfreqs)
len(normPfreqs)
len(clusters)
len(bands)
len(alleleFreqSUMS)

type(tdaGenes[0])              # STR elements are gene names
type(normBfreqs[0])     # FLOAT
type(normLBfreqs[0])    # FLOAT
type(normLPfreqs[0])    # FLOAT    
type(normPfreqs[0])     # FLOAT
type(clusters[0])           # STR elements are TDA cluster names.  'Cluster 0' ... 'Cluster 5'
type(bands[0])              # STR elements are guildmembers to TDA clusters. In this case, a band is a gene
type(alleleFreqSUMS[0])     # FLOAT

### In memory space, these are the current useful listVars.  All are reference the LEFT JOIN mergeDF:
# tdaGenes
# normBfreqs
# normLBfreqs
# normLPfreqs
# normPfreqs
# clusters
# bands
# alleleFreqSUMS

### In memory space, these are the current useful dfVars.  Note that allelCountsDF has an appended column.
# tdaDF                # [245 rows x 3 columns]
# alleleCountsDF       # [199 rows x 6 columns]   <- the last column is 'alleleSUM' appended in this script
# mergeDF              # [245 rows x 8 columns]



# Programmaticaly create a new listVar that tracks TDA cluster names as string elements.
# Specify how many TDA clusters (in the case of the ClinVar anemia dataset, there are six clusters.
numClusters = 6
### In the future, you can check the mergeDF.cluster column for unique cluster names.
# counter i tracks a cluster iteration
i = 0
clusterNames = []
for i in range(numClusters):
    # print("Cluster " + str(i))
    clusterNames.append('Cluster ' + str(i))

# View listVar of clusterName handles
clusterNames
# ['Cluster 0', 'Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5']
len(clusterNames)
# 6


# Next, calculate average column values for B_freq, LB_freq, LP_freq, P_freq on a per TDA cluster basis

# Reset i counter
i = 0

# Flush accumulator listVars
clusterAvgBfreq = []
clusterAvgLBfreq = []
clusterAvgLPfreq = []
clusterAvgPfreq = []
clusterAvgAlleleFreqSUM = []
clusterNumGenes = []

for i in range(len(clusterNames)):
    # Retrive a cluster name handle
    clusterName = clusterNames[i]
    # Subset the master mergeDF for only genes that belong to the current TDA cluster
    clusterDF = mergeDF[ mergeDF['cluster'] == clusterName ]
    # Create a non_numeric mask listvar
    non_numeric = ['cluster', 'band', 'gene_name']
    # Generate a numeric only version of clusterDF
    numeric_clusterDF = clusterDF.drop(non_numeric, axis=1)
    # Calculate col_avgs for numeric_clusterDF
    col_avgs = numeric_clusterDF.mean(axis=0)
    # Convert pandas series to a list of int elements:
    col_avgs = col_avgs.tolist()
    # Append current cluster categorical allele average frequencies to accumulator listVars:
    clusterAvgBfreq.append(col_avgs[0])
    clusterAvgLBfreq.append(col_avgs[1])
    clusterAvgLPfreq.append(col_avgs[2])
    clusterAvgPfreq.append(col_avgs[3])
    clusterAvgAlleleFreqSUM.append(col_avgs[4])
    # Store the number of genes found in current TDA cluster into an accumulator listVar:
    clusterNumGenes.append( len(numeric_clusterDF) )

# LEN = 6
len(clusterAvgBfreq)            
len(clusterAvgLBfreq)
len(clusterAvgLPfreq)
len(clusterAvgPfreq)
len(clusterAvgAlleleFreqSUM)
len(clusterNumGenes)

# TYPE = lists of floats
type(clusterAvgBfreq[0])
type(clusterAvgLBfreq[0])
type(clusterAvgLPfreq[0])
type(clusterAvgPfreq[0])
type(clusterAvgAlleleFreqSUM[0])
# TYPE = list of INT
type(clusterNumGenes[0])

# VIEW
clusterAvgBfreq
clusterAvgLBfreq
clusterAvgLPfreq
clusterAvgPfreq
clusterAvgAlleleFreqSUM
clusterNumGenes



############# The current cluster listVars are parsed based on listVars tracking features across clusters.
### To graph, I need to rearray the listVars to collate features that are cluster-specific

#   CATEGORY: Benign                (dark blue)     #29386F
#   CATEGORY: Likely benign         (light blue)    #DFEDFA
#   CATEGORY: Likely pathogenic     (light red)     #FCE5EA
#   CATEGORY: Pathogenic            (dark red)      #802A2A   

# Define a custom SNS palette
customPalette = ['#29386F', '#DFEDFA', '#FCE5EA', '#802A2A']

# flush counter i which will track the current cluster
i = 0
# Reuse user-defined TDA cluster number variable (initialized earlier)
numClusters     # INT
# Resuse clusterNames listVar (initialized earlier)
clusterNames    # LIST of string elements


for i in range(len(clusterNames)):
    # Fetch current clusterName as a string
    clusterName = clusterNames[i]
    # Fetch number of genes found in current cluster
    clusterGeneNumber = clusterNumGenes[i]
    # Flush newDF
    newDF = pd.DataFrame()
    # Populate current dataframe with features specific to current TDA cluster 
    newDF['avgB_freq'] = [ clusterAvgBfreq[i] ]
    newDF['avgLB_freq'] = [ clusterAvgLBfreq[i] ]
    newDF['avgLP_freq'] = [ clusterAvgLPfreq[i] ]
    newDF['avgP_freq'] = [ clusterAvgPfreq[i] ]
    # Make a simple Seaborn bar graph
    fileTypes = ['png', 'svg']
    # Reset j counter, which tracks the file type output
    j = 0
    for j in range(len(fileTypes)):
        # Build current filename
        filename = 'TDA_node' + str(i) + '_avgAlleleFreqs.' + fileTypes[j]
        ax = sns.barplot(data=newDF, palette = customPalette, linewidth=0.8, edgecolor='#29386F')
        plt.ylim(top=1.0) 
        plt.xlabel('ClinVar allelic categories')
        plt.ylabel('ClinVar allelic frequencies')
        plt.xticks(fontsize=8)  
        plt.yticks(fontsize=8)
        plt.legend( title = clusterName + ': ' + str(clusterGeneNumber) + ' genes' )
        plt.title(clusterName)
        # plt.show()
        plt.savefig('./out/' + filename)
        ## Clear plot space
        plt.clf()
        plt.cla()
        plt.close()