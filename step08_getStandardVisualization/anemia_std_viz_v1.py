## The following python script performs a standard visualizatin pipeline when given a ClinVar raw allele counts table as a CSV file in this format:
    # n rows where n = total gene number
    # 5 columns with labels: gene_name   B   LB   LP    P
        # B = benign
        # LB = likely benign
        # LP = likely pathogenic
        # P = pathogenic
    # NOTA BENE: The CSV file must contain the original allele counts (that are not normalized frequencies).

# The output of this script follows this general standard analysis pipeline.  CSV outs + PNG/SVG plots are written to working directory:
    # 1. A raw allele counts (CSV) file that contains the following columns: 
        #   # (an index column)
        #   Gene name
        #   CATEGORY: Benign                (dark blue)     #29386F
        #   CATEGORY: Likely benign         (light blue)    #DFEDFA
        #   CATEGORY: Likely pathogenic     (light red)     #FCE5EA
        #   CATEGORY: Pathogenic            (dark red)      #802A2A
        #   Allele SUM
        #   Log2 of Allele SUM

    # 2. A frequency-normalized allele counts (CSV) file that contains the following columns.  Allele frequencies are calculated on a per-locus basis.
        #   # (an index column)
        #   Gene name
        #   B.freq                          (dark blue)     #29386F
        #   LB.freq                         (light blue)    #DFEDFA
        #   LP.freq                         (light red)     #FCE5EA
        #   P.freq                          (dark red)      #802A2A
        #   SUM of categorical allele freqs

    # 3. Bar plot of log2-transformed allele sums (B+LB+LP+P) vs. genes (n = 199 loci)
    
    # 4. Histogram (10% bin sizes) of total allele counts (B+LB+LP+P)
    
    # 5. PIE chart of genes sorted by total allele number (with top20 genes as a 10% slice)
    
    # 6. Horizontal bar chart showing the sorted list of genes by top allele counts (descending sort).
        #   For each gene, raw number of each of 4 ClinVar categories: B, LB, LP, P.  Follow this color schema:
                #   CATEGORY: Benign                (dark blue)     #29386F
                #   CATEGORY: Likely benign         (light blue)    #DFEDFA
                #   CATEGORY: Likely pathogenic     (light red)     #FCE5EA
                #   CATEGORY: Pathogenic            (dark red)      #802A2A       

    # 7. PIE chart of genes sorted by normalized PROBLEMATIC allele categorical frequencies (with top20 genes as a 10% slice)
    
    # 8. Horizontal bar chart showing the sorted list of genes by top allele frequencies (descending sort of P then LP then LB, then B).
        #   For each gene, raw number of each of 4 ClinVar categories: B, LB, LP, P.  Follow this color schema:
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

### After VENV vPCA construciton, upon WSL bootup:
source gotovenv.sh
source vPCA/bin/activate
cd /mnt/c/wslshare/github/inherited_anemias/step08_getStandardVisualization
python

############################################################ IMPORTS
import pandas as pd
import numpy as np
import math
import seaborn as sns
import matplotlib.pyplot as plt
# from sklearn.manifold import TSNE

# LOAD DATA
alleles_df = pd.DataFrame()                             # Note: alleles_df will never be sorted
alleles_df = pd.read_csv('Anemia_199_genes_AlleleCount_4_ClinVar_classifiers.csv')
# Extract columns into list variables:
genes = alleles_df['gene_name'].tolist() 
rawBalleles = alleles_df['B'].tolist() 
rawLBalleles = alleles_df['LB'].tolist() 
rawLPalleles = alleles_df['LP'].tolist() 
rawPalleles = alleles_df['P'].tolist() 

# Pull out numeric only columns (remove the 'gene_name' column)
# non_numeric = ['gene_name']
# numeric_alleles_df = alleles_df.drop(non_numeric, axis=1)

########################################################### step 1: Table of raw allele counts

# Create raw allele counts file (CSV):
        # 1. A raw allele counts (CSV) file that contains the following columns: 
        #   # (an index column)
        #   Gene name
        #   CATEGORY: Benign                (dark blue)     #29386F
        #   CATEGORY: Likely benign         (light blue)    #DFEDFA
        #   CATEGORY: Likely pathogenic     (light red)     #FCE5EA
        #   CATEGORY: Pathogenic            (dark red)      #802A2A
        #   Allele SUM
        #   Log2 of Allele SUM


# Use pd.sum function to calculate rowSUMS
alleles_df['alleleSUM'] = alleles_df.sum(axis=1, numeric_only=True)
# Generate a listvar that tracks allele SUMS per locus:
alleleSUMS = alleles_df['alleleSUM'].tolist()       # INT type list
# Preview df
alleles_df
    # gene_name    B    LB   LP    P  alleleSUM
# 0       ABCA1  250   580   18   71        919
# 1       ABCG8   81   202   24   45        352
# 2      ACVRL1   94   191  168  373        826
# 3        ADA2   44   150   31  121        346
# 4        ALG8   71    86   20   32        209
# ..        ...  ...   ...  ...  ...        ...
# 194     MT-TV    5     5    4    3         17
# 195   MT-RNR1   28    34    3    3         68
# 196   MT-RNR2    0     1    1    3          5
# 197      POLG  125  1049  180  246       1600
# 198     POLG2   25   129    7   28        189
# [199 rows x 6 columns]

# Write to disk (CSV file)
alleles_df = alleles_df.set_index('gene_name')
alleles_df.to_csv('anemia_ClinVar_B_LB_LP_P_raw_allele_counts.csv')



########################################################### step 2: Table of categorical allele frequencies (normalized per locus)

# Create transformed categorical allele frequencies file (CSV):
    # 2. A frequency-normalized allele counts (CSV) file that contains the following columns.  Allele frequencies are calculated on a per-locus basis.
        #   # (an index column)
        #   Gene name
        #   B.freq                          (dark blue)     #29386F
        #   LB.freq                         (light blue)    #DFEDFA
        #   LP.freq                         (light red)     #FCE5EA
        #   P.freq                          (dark red)      #802A2A
        #   SUM of categorical allele freqs

## Loop over list of genes
    ## Fetch a gene
        ## Fetch raw counts of B, LB, LP, and P alleles
        ## Fetch the alleleSUM for current gene
        ## Normalize B, LB, LP, and P to arrive at categorical allele frequencies

## Useful listvars:
genes           # str elems
rawBalleles         # INT
rawLBalleles        # INT
rawLPalleles        # INT
rawPalleles         # INT
alleleSUMS          # INT

len(genes)  # 199 loci
i = 0

## Accumulator listvars to store a list of float elements:
freqBalleles = []         # floats
freqLBalleles = []        # floats
freqLPalleles = []        # floats
freqPalleles = []         # floats

for i in range(len(genes)):
    # get first gene in list
    gene = genes[i]
    # get raw allele counts (numerators)
    rawB = rawBalleles[i]
    rawLB = rawLBalleles[i]
    rawLP = rawLPalleles[i]
    rawP = rawPalleles[i]
    # get denominator (total allele counts across 4 ClinVar categories: B, LB, LP, and P)
    denom = alleleSUMS[i]
    # Append normalized frequences as float elements to accumulator listvars
    freqBalleles.append(rawB/denom)
    freqLBalleles.append(rawLB/denom)
    freqLPalleles.append(rawLP/denom)
    freqPalleles.append(rawP/denom)

# Create a new dataframe to capture categorical allele frequencies    
dfAlleleFreqs = pd.DataFrame()
# Write columns to pandas dataframe:
dfAlleleFreqs['gene_name'] = genes
dfAlleleFreqs['B.freq'] = freqBalleles
dfAlleleFreqs['LB.freq'] = freqLBalleles
dfAlleleFreqs['LP.freq'] = freqLPalleles
dfAlleleFreqs['P.freq'] = freqPalleles
# QC step to ensure that allele frequency SUMS each equal to 1.0
# Use pd.sum function to calculate rowSUMS
dfAlleleFreqs['alleleFreqSUM'] = dfAlleleFreqs.sum(axis=1, numeric_only=True)
# Generate a listvar that tracks allele frequency SUMS per locus:
alleleFreqSUMS = dfAlleleFreqs['alleleFreqSUM'].tolist()       # float type list
# Set index to gene_name
dfAlleleFreqs = dfAlleleFreqs.set_index('gene_name')
# Write pd to disk
dfAlleleFreqs.to_csv('anemia_ClinVar_B_LB_LP_P_categorical_allele_freqs.csv')



########################################################### step 3: bar plots of log2-transformed total allele counts (B+LB+LP+P)
    
    # 3. Bar plot of log2-transformed allele sums (B+LB+LP+P) vs. genes (n = 199 loci)

#### Useful variables in memory:
# labels
genes               # str 
# raw allele counts:
rawBalleles         # INT
rawLBalleles        # INT
rawLPalleles        # INT
rawPalleles         # INT
alleleSUMS          # INT
# normalized allele frequencies:
freqBalleles        # floats
freqLBalleles       # floats
freqLPalleles       # floats
freqPalleles        # floats
alleleFreqSUMS      # floats
# pd dataframes:
alleles_df          # 199 (genes) x 5 cols (raw allele counts: B, LB, LP, P, SUM).  gene_name is set to df index
dfAlleleFreqs       # 199 (genes) x 5 cols (norm allele counts: B.freq, LB.freq, LP.freq, P.freq, FreqSUM).  gene_name is set to df index

# View a summary of raw allele counts dataframe:
alleles_df.describe()
                # B           LB          LP            P     alleleSUM
# count  199.000000   199.000000  199.000000   199.000000    199.000000
# mean    62.055276   289.552764   54.994975   158.894472    565.497487
# std     97.449712   481.069102   94.265770   474.646486   1061.639084
# min      0.000000     1.000000    0.000000     1.000000      5.000000
# 25%     18.000000    30.500000    5.000000    18.000000     95.500000
# 50%     38.000000   130.000000   18.000000    48.000000    235.000000
# 75%     69.000000   360.000000   62.500000   139.000000    673.000000
# max    845.000000  3752.000000  638.000000  5005.000000  10240.000000

## Across all 199 loci, there is thus a large spread (variance) for total allele counts (alleleSUM), from 5-10K.
## To view this together, we can perform a log2 transformation of the total allele counts per genetic locus.

# Use math.log(x, base)
math.log(8,2)   # returns 3.0, which means that 2^3.0 = 8

# Reset counter
i = 0
# Set an accumulator listvar
log2_alleleSUMS = []    

for i in range(len(genes)):
    # fetch the current gene's alleleSUM
    alleleSum = alleleSUMS[i]
    # calc the log2-value of the current locus' total allele count
    log2_TFN = math.log(alleleSum, 2)
    log2_alleleSUMS.append(log2_TFN)

len(log2_alleleSUMS)    # 199 loci

# Create a copy of alleles_df
log2df = alleles_df
# Append the log2-transformed alleleSUM column
log2df['log2_alleleSUM'] = log2_alleleSUMS
# Sort the log2df by log2_alleleSUM, in descending order
log2df_sorted_byLOG2 = log2df.sort_values(by='log2_alleleSUM', ascending=False)


# preview log2df
log2df
             # B    LB   LP    P  alleleSUM  log2_alleleSUM
# gene_name
# ABCA1      250   580   18   71        919        9.843921
# ABCG8       81   202   24   45        352        8.459432
# ACVRL1      94   191  168  373        826        9.689998
# ADA2        44   150   31  121        346        8.434628
# ALG8        71    86   20   32        209        7.707359
# ...        ...   ...  ...  ...        ...             ...
# MT-TV        5     5    4    3         17        4.087463
# MT-RNR1     28    34    3    3         68        6.087463
# MT-RNR2      0     1    1    3          5        2.321928
# POLG       125  1049  180  246       1600       10.643856
# POLG2       25   129    7   28        189        7.562242
# [199 rows x 6 columns]


# preview log2df_sorted_byLOG2
log2df_sorted_byLOG2
             # B    LB   LP     P  alleleSUM  log2_alleleSUM
# gene_name
# BRCA2      845  3752  638  5005      10240       13.321928
# BRCA1      805  2554  462  4065       7886       12.945078
# COL7A1     128  2895  287   645       3955       11.949462
# FANCA      251  2013  430   933       3627       11.824561
# PALB2      119  1305  390  1246       3060       11.579316
# ...        ...   ...  ...   ...        ...             ...
# HLA-DQA1     4     4    1     4         13        3.700440
# MT-TM        6     1    3     3         13        3.700440
# MT-TP        6     2    2     2         12        3.584963
# HLA-DQB1     2     4    1     4         11        3.459432
# MT-RNR2      0     1    1     3          5        2.321928
# [199 rows x 6 columns]


## Barplot (PNG)
ax = sns.barplot(data=log2df_sorted_byLOG2, x='gene_name', y='log2_alleleSUM', color='#29386F')
ax.set_xticklabels(ax.get_xticklabels(), ha="right")
plt.xlabel('gene name\n(n = 199 loci)')
plt.ylabel('log2-transformation of total allele counts\n(B+LB+LP+P)')
plt.xticks(fontsize=8, rotation=45)  # Rotate labels by 45 degrees, set font to 8
plt.xticks(range(0, len(genes), 6))  # Show every 6th label (POLG visible)
plt.yticks(fontsize=8)
plt.subplots_adjust(bottom=0.2, left=0.114, right=0.945, top=0.802)  # Adjust the value as needed
plt.savefig('log2barplot.png')
# plt.show()
## Clear plot space
plt.clf()
plt.cla()
plt.close()

## Barplot (SVG)
ax = sns.barplot(data=log2df_sorted_byLOG2, x='gene_name', y='log2_alleleSUM', color='#29386F')
ax.set_xticklabels(ax.get_xticklabels(), ha="right")
plt.xlabel('gene name\n(n = 199 loci)')
plt.ylabel('log2-transformation of total allele counts\n(B+LB+LP+P)')
plt.xticks(fontsize=8, rotation=45)  # Rotate labels by 45 degrees, set font to 8
plt.xticks(range(0, len(genes), 6))  # Show every 6th label (POLG visible)
plt.yticks(fontsize=8)
plt.subplots_adjust(bottom=0.2, left=0.114, right=0.945, top=0.802)  # Adjust the value as needed
plt.savefig('log2barplot.svg')
# plt.show()
## Clear plot space
plt.clf()
plt.cla()
plt.close()




########################################################### step 4: histogram plots of log2-transformed total allele counts (B+LB+LP+P)

    # 4. Histogram (10% bin sizes) of total allele counts (B+LB+LP+P)

## Histogram (PNG)
ax = sns.histplot(data=log2df_sorted_byLOG2, x="log2_alleleSUM", bins = 10, color='#802A2A', alpha=0.92)
plt.xlabel('ten percentile bins of log2-transformed allele sums', fontsize=10)
plt.ylabel('gene count', fontsize=10)
plt.savefig('log2histogram.png')
# plt.show()
plt.clf()
plt.cla()
plt.close()

## Histogram (SVG)
ax = sns.histplot(data=log2df_sorted_byLOG2, x="log2_alleleSUM", bins = 10, color='#802A2A', alpha=0.92)
plt.xlabel('ten percentile bins of log2-transformed allele sums', fontsize=10)
plt.ylabel('gene count', fontsize=10)
plt.savefig('log2histogram.svg')
# plt.show()
plt.clf()
plt.cla()
plt.close()


########################################################### step 5: PIE chart showing top20 genes for total allele counts (B+LB+LP+P)

    # 5. PIE chart of genes sorted by total allele number (with top20 genes as a 10% slice)

## Of 199 anemia-enriched loci, let's display the top20 genes sorted by total allele counts (B+LB+LP+P).  Top20 = ~10% of the 199 loci.

# Let's create a dictionary relating categorical labels to %
data = {'category': ['top 20 genes', ''], 'values': [10, 90]}
# Initialize a pd df using this simple dictionary
df = pd.DataFrame(data)
# View the pd df
df
  # category  values
# 0        A      10
# 1        B      90
# set colors
colorNames = ['black',   'grey',    'blue',    'light blue',  'red',      'light red',  'green',   'light green']
colors = [    '#363636', '#E5E5E5', '#29386F', '#DFEDFA',     '#802A2A',  '#FCE5EA',    '#44A043',  '#E0F4DA']
# declaring exploding pie 
explode = [0.15, 0]

## Get PNG
# Set Seaborn style
sns.set_theme()
sns.set_style("whitegrid")
# Create pie chart
plt.figure(figsize=(9, 6))
plt.pie(df['values'], labels=df['category'], explode=explode, autopct='%1.1f%%', startangle = 345, pctdistance=1.15, labeldistance=1.5, colors=colors)
plt.title('Genes sorted by total allele number')
plt.tight_layout()
# plt.show()
plt.savefig('pieTop20rawCounts.png')
# Clear plot space
plt.clf()
plt.cla()
plt.close()

## Get SVG
# Set Seaborn style
sns.set_theme()
sns.set_style("whitegrid")
# Create pie chart
plt.figure(figsize=(9, 6))
plt.pie(df['values'], labels=df['category'], explode=explode, autopct='%1.1f%%', startangle = 345, pctdistance=1.15, labeldistance=1.5, colors=colors)
plt.title('Genes sorted by total allele number')
plt.tight_layout()
# plt.show()
plt.savefig('pieTop20rawCounts.svg')
# Clear plot space
plt.clf()
plt.cla()
plt.close()


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

########################################################### step 6: Horizontal bar chart showing  sorted list of genes by top allele counts (descending sort).

    # 6. Horizontal bar chart showing the sorted list of genes by top allele counts (descending sort).
        #   For each gene, raw number of each of 4 ClinVar categories: B, LB, LP, P.  Follow this color schema:
                #   CATEGORY: Benign                (dark blue)     #29386F
                #   CATEGORY: Likely benign         (light blue)    #DFEDFA
                #   CATEGORY: Likely pathogenic     (light red)     #FCE5EA
                #   CATEGORY: Pathogenic            (dark red)      #802A2A     


# Slice the top20 loci from the sorted log2df, and store in a pd variable
top20byAlleleCount = log2df_sorted_byLOG2[0:20]
topBvals = top20byAlleleCount['B'].tolist()
topLBvals = top20byAlleleCount['LB'].tolist()
topLPvals = top20byAlleleCount['LP'].tolist()
topPvals = top20byAlleleCount['P'].tolist()
topGeneNames = top20byAlleleCount.index        # index = gene names

df = pd.DataFrame()
df['gene_name'] = topGeneNames
df['B'] = topBvals
df['LB'] = topLBvals
df['LP'] = topLPvals
df['P'] = topPvals


############################################################# Left off here!



####################################################################################################################################################################################
# https://stackoverflow.com/questions/61740434/how-to-make-horizontal-bar-chart-using-seaborn-to-be-stacked-to-a-100-by-hue-pa

df['pct'] = df['Count'] / sum(df['Count'])

df.pivot(index=topGeneNames,columns='Count',
         values='pct').plot(kind='barh',stacked=True,colormap='Blues',rot=90) 


####################################################################################################################################################################################

## Barplot (PNG)
ax = sns.barplot(data=df, y='gene_name', color='#29386F')
ax.set_xticklabels(ax.get_xticklabels(), ha="right")
plt.xlabel('gene name\n(n = 199 loci)')
plt.ylabel('log2-transformation of total allele counts\n(B+LB+LP+P)')
plt.xticks(fontsize=8, rotation=45)  # Rotate labels by 45 degrees, set font to 8
plt.xticks(range(0, len(genes), 6))  # Show every 6th label (POLG visible)
plt.yticks(fontsize=8)
plt.subplots_adjust(bottom=0.2, left=0.114, right=0.945, top=0.802)  # Adjust the value as needed
# plt.savefig('log2barplot.png')
plt.show()
## Clear plot space
# plt.clf()
# plt.cla()
# plt.close()





########################################################### step 7: PIE chart showing top20 genes for allele freqs (LP+P)

    # 7. PIE chart of genes sorted by problematic allele frequencies (with top20 genes as a 10% slice)

## Of 199 anemia-enriched loci, let's display the top20 genes sorted by total allele counts (LP+P).  Top20 = ~10% of the 199 loci.

# Let's create a dictionary relating categorical labels to %
data = {'category': ['top 20 genes', ''], 'values': [10, 90]}
# Initialize a pd df using this simple dictionary
df = pd.DataFrame(data)
# View the pd df
df
  # category  values
# 0        A      10
# 1        B      90

# set colors
colorNames = ['grey',   'black',    'blue',    'light blue',  'red',      'light red',  'green',   'light green']
colors = [    '#E5E5E5', '#363636', '#29386F', '#DFEDFA',     '#802A2A',  '#FCE5EA',    '#44A043',  '#E0F4DA']
# declaring exploding pie 
explode = [0.15, 0]

## Get PNG
# Set Seaborn style
sns.set_theme()
sns.set_style("whitegrid")
# Create pie chart
plt.figure(figsize=(9, 6))
plt.pie(df['values'], labels=df['category'], explode=explode, autopct='%1.1f%%', startangle = 165, pctdistance=1.15, labeldistance=1.5, colors=colors)
plt.title('Genes sorted by problematic allele frequencies')
plt.tight_layout()
# plt.show()
plt.savefig('pieTop20normFreqs.png')
# Clear plot space
plt.clf()
plt.cla()
plt.close()

## Get SVG
# Set Seaborn style
sns.set_theme()
sns.set_style("whitegrid")
# Create pie chart
plt.figure(figsize=(9, 6))
plt.pie(df['values'], labels=df['category'], explode=explode, autopct='%1.1f%%', startangle = 165, pctdistance=1.15, labeldistance=1.5, colors=colors)
plt.title('Genes sorted by problematic allele frequencies')
plt.tight_layout()
# plt.show()
plt.savefig('pieTop20normFreqs.svg')
# Clear plot space
plt.clf()
plt.cla()
plt.close()

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
















############################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################




alleles_df.shape
# (199, 7)                    # returns a tuple which means the df has 199 rows and 7 columns

alleles_df.head()              # view just the first five genes of 199 genes
  # gene_name   CC    B   LB   US   LP    P
# 0     ABCA1  113  250  580  467   18   71
# 1     ABCG8   85   81  202  358   24   45
# 2    ACVRL1   52   94  191  208  168  373
# 3      ADA2   14   44  150  242   31  121
# 4      ALG8   18   71   86  127   20   32

alleles_df.describe()
# returns simple stats analysis:
                # CC           B           LB           US          LP            P
# count   199.000000  199.000000   199.000000   199.000000  199.000000   199.000000
# mean     79.407035   62.055276   289.552764   321.407035   54.994975   158.894472
# std     411.702661   97.449712   481.069102   528.190890   94.265770   474.646486
# min       0.000000    0.000000     1.000000     4.000000    0.000000     1.000000
# 25%       1.500000   18.000000    30.500000    50.500000    5.000000    18.000000
# 50%      16.000000   38.000000   130.000000   156.000000   18.000000    48.000000
# 75%      47.000000   69.000000   360.000000   330.000000   62.500000   139.000000
# max    5092.000000  845.000000  3752.000000  3970.000000  638.000000  5005.000000

alleles_df.loc[0]   # gets first of 199 genes
# gene_name    ABCA1
# CC             113
# B              250
# LB             580
# US             467
# LP              18
# P               71
# Name: 0, dtype: object

alleles_df.loc[1]   # gets second of 199 genes
# gene_name    ABCG8
# CC              85
# B               81
# LB             202
# US             358
# LP              24
# P               45
# Name: 1, dtype: object

alleles_df.loc[1][1]
# np.int64(85)                  # 85 alleles of this categorical type
alleles_df.loc[1][2]
# np.int64(81)                  # 81 alleles of this categorical type
alleles_df.loc[1][3]
# np.int64(202)                 # 202 alleles of this categorical type
alleles_df.loc[1][4]
# np.int64(358)                 # 358 alleles of this categorical type
alleles_df.loc[1][5]
# np.int64(24)                  # 24 alleles of this categorical type
alleles_df.loc[1][6]
# np.int64(45)                  # 45 alleles of this categorical type

############################################################ NUMERIC DF INITALIZIATION
# Pull out numeric only columns (remove the 'gene_name' column)
# non_numeric = ['gene_name']
# numeric_alleles_df = alleles_df.drop(non_numeric, axis=1)

numeric_alleles_df.shape
# (199, 6)
numeric_alleles_df.head()
    # CC    B   LB   US   LP    P
# 0  113  250  580  467   18   71
# 1   85   81  202  358   24   45
# 2   52   94  191  208  168  373
# 3   14   44  150  242   31  121
# 4   18   71   86  127   20   32

numeric_alleles_df.describe()
                # CC           B           LB           US          LP            P
# count   199.000000  199.000000   199.000000   199.000000  199.000000   199.000000
# mean     79.407035   62.055276   289.552764   321.407035   54.994975   158.894472
# std     411.702661   97.449712   481.069102   528.190890   94.265770   474.646486
# min       0.000000    0.000000     1.000000     4.000000    0.000000     1.000000
# 25%       1.500000   18.000000    30.500000    50.500000    5.000000    18.000000
# 50%      16.000000   38.000000   130.000000   156.000000   18.000000    48.000000
# 75%      47.000000   69.000000   360.000000   330.000000   62.500000   139.000000
# max    5092.000000  845.000000  3752.000000  3970.000000  638.000000  5005.000000

############################################################ SNS PAIRPLOT
## Run seaborn's pairplot analysis
# sns.pairplot(numeric_alleles_df, diag_kind ='hist')
## Make sure on Win11/WSL2 to first run XLaunch, and ensure that 'disable access controls' is toggled ON in last menu dialog
# plt.show()
## View the plot from memory

# Run seaborn's pairplot analysis
sns.pairplot(numeric_alleles_df, diag_kind ='hist')
# Make sure on Win11/WSL2 to first run XLaunch, and ensure that 'disable access controls' is toggled ON in last menu dialog
# SAVE the pairplot to file system
plt.savefig('raw_pairplot.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
# get a PNG version
sns.pairplot(numeric_alleles_df, diag_kind ='hist')
plt.savefig('raw_pairplot.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############################################################ TSNE
# t-SNE
m = TSNE(learning_rate=50)
tsne_features = m.fit_transform(numeric_alleles_df)
len(tsne_features)
# 164
tsne_features[1:4,:]
# array([[ 2.0416634 ,  2.980683  ],
       # [ 0.77361345, -3.9859142 ],
       # [-0.25464657,  1.362818  ]], dtype=float32)
# Initialize tsne df
tsne_df = pd.DataFrame()
# Transfer tsne output into tsne df
tsne_df['x'] = tsne_features[:,0]
tsne_df['y'] = tsne_features[:,1]

# Setup an sns scatterplot to visualize tsne analysis
sns.scatterplot( x="x", y="y", data=tsne_df )
# sns.scatterplot( x="x", y="y", hue="y", palette=sns.color_palette("hls", 10), data=tsne_df, legend="full", alpha=0.3 )
plt.savefig('TSNE_scatterplot.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

# get a PNG version
sns.scatterplot( x="x", y="y", data=tsne_df )
plt.savefig('TSNE_scatterplot.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

tsne_df.var()
# x    72.292007
# y     7.238809
# dtype: float32
# Boxplot
sns.boxplot(tsne_df)
# plt.show()
plt.savefig('TSNE_boxplot.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
# get a PNG version
sns.boxplot(tsne_df)
plt.savefig('TSNE_boxplot.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############################################################ PAIRWISE CORRELATION
# The correlation coefficient (r) measures the strength of correlation between values of different dimensions (columns)
# r is measured from -1 to 0 to +1
# if r = -1, then the two dimensions are negatively correlated
# if r = 0, then there is no correlation between the two dimensions
# if r = +1, then the two dimensions are positively correlated
## As a rule, dimensions of the data that exhibit either r = -1 or r = +1 show perfect correlation and thus do not add any new variance to the entire dataset. They do not increase the complexity of the dataset, and including both of these redundant dimensions risks lowering the accuracy of test/train ML and thus may lead to overfitting the data.
## Overfit is measured by the accuracy difference between TRAIN and TEST runs.  If the accuracy differential is too high then the fit either favors (too precisely, too specific) the original TRAINING dataset, or it favors (too broadly, too generalized) the test dataset.  Ideally, a good statistical fitted model means the TRAIN accuracy = the TEST accuracy score.
## Goal: remove dimensions that exhibit close to -1 or +1 correlation coefficients (r).  This reduces the risk of overfitting a model, while minimizing loss of complexity to the original dataset.

numeric_alleles_df.corr()
          # CC         B        LB        US        LP         P
# CC  1.000000  0.819468  0.717500  0.676335  0.626177  0.947301
# B   0.819468  1.000000  0.792313  0.718777  0.673274  0.834417
# LB  0.717500  0.792313  1.000000  0.866518  0.797543  0.762214
# US  0.676335  0.718777  0.866518  1.000000  0.736070  0.697099
# LP  0.626177  0.673274  0.797543  0.736070  1.000000  0.756673
# P   0.947301  0.834417  0.762214  0.697099  0.756673  1.000000
cmap = sns.diverging_palette ( h_neg=10, h_pos=240, as_cmap=True)
sns.heatmap(numeric_alleles_df.corr(), center=0, cmap=cmap, linewidths=1, annot=True, fmt=".2f")
plt.savefig('Correlation_heatmap.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

# get a PNG version
numeric_alleles_df.corr()
          # CC         B        LB        US        LP         P
# CC  1.000000  0.819468  0.717500  0.676335  0.626177  0.947301
# B   0.819468  1.000000  0.792313  0.718777  0.673274  0.834417
# LB  0.717500  0.792313  1.000000  0.866518  0.797543  0.762214
# US  0.676335  0.718777  0.866518  1.000000  0.736070  0.697099
# LP  0.626177  0.673274  0.797543  0.736070  1.000000  0.756673
# P   0.947301  0.834417  0.762214  0.697099  0.756673  1.000000
cmap = sns.diverging_palette ( h_neg=10, h_pos=240, as_cmap=True)
sns.heatmap(numeric_alleles_df.corr(), center=0, cmap=cmap, linewidths=1, annot=True, fmt=".2f")
plt.savefig('Correlation_heatmap.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

# create a MASK to remove the upper triangle
corr = numeric_alleles_df.corr()
mask = np.triu( np.ones_like(corr, dtype = bool) )
cmap = sns.diverging_palette ( h_neg=10, h_pos=240, as_cmap=True)
sns.heatmap(numeric_alleles_df.corr(), mask=mask, center=0, cmap=cmap, linewidths=1, annot=True, fmt=".2f")
plt.savefig('Correlation_heatmap_LT.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

# get a PNG version
# create a MASK to remove the upper triangle
corr = numeric_alleles_df.corr()
mask = np.triu( np.ones_like(corr, dtype = bool) )
cmap = sns.diverging_palette ( h_neg=10, h_pos=240, as_cmap=True)
sns.heatmap(numeric_alleles_df.corr(), mask=mask, center=0, cmap=cmap, linewidths=1, annot=True, fmt=".2f")
plt.savefig('Correlation_heatmap_LT.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()


################ PRINCIPAL COMPONENT ANALYSIS (PCA) codeblock ###########################                                                    

# Using a StandardScaler approach:
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
numeric_alleles_df_std = pd.DataFrame(scaler.fit_transform(numeric_alleles_df),columns=numeric_alleles_df.columns)

# View the fit/transformation using the elected scaler:
numeric_alleles_df_std.head()
         # CC         B        LB        US        LP         P
# 0  0.081801  1.933497  0.605276  0.276340 -0.393444 -0.185646
# 1  0.013619  0.194895 -0.182455  0.069455 -0.329633 -0.240562
# 2 -0.066738  0.328634 -0.205379 -0.215250  1.201815  0.452222
# 3 -0.159270 -0.185745 -0.290820 -0.150717 -0.255188 -0.080039
# 4 -0.149530  0.092020 -0.424193 -0.368990 -0.372174 -0.268019

## compare with original (pre-fit/transform), which represent raw allele counts:
numeric_alleles_df.head()
    # CC    B   LB   US   LP    P
# 0  113  250  580  467   18   71
# 1   85   81  202  358   24   45
# 2   52   94  191  208  168  373
# 3   14   44  150  242   31  121
# 4   18   71   86  127   20   32

# Bring in PCA tools:
from sklearn.decomposition import PCA

# Declare a pca object of the PCA class
pca = PCA()

# Fit the standardized dataset to perform PCA:
pca.fit(numeric_alleles_df_std)

# Get the PCA explained variance ratio
print(pca.explained_variance_ratio_)
# This tells me individual PC contributions to the overall PCA.  Usually the first 2-3 components contain most of dataset variance
# [0.8018952  0.09095647 0.05025481 0.03302067 0.01863209 0.00524076]

# Get the accumulated PCA explaind variance ratio
print(pca.explained_variance_ratio_.cumsum())
# [0.8018952  0.89285166 0.94310647 0.97612715 0.99475924 1.        ]
## This tells me that 94.2% of the data is explained by the first 3 PCs of PCA

# Get the PCA component ruleset:
print(pca.components_)
## This lists all principal components:
#       CC              B           LB          US          LP          P           # VARIANCE EXPLANATION BY PRINCIPAL COMPONENTS:
# PC1 [ 0.40726418  0.41126905  0.41883331  0.39776639  0.3883725   0.42487764]     # All classes contribute relatively equally to the variances
# PC2 [ 0.54177701  0.24427407 -0.35049018 -0.4475752  -0.4155178   0.38856901]     # CC / P / B are contributing together, inversely to other three classifiers
# PC3 [-0.02301058 -0.28973272 -0.22238972 -0.47216185  0.75324028  0.27524633]     # LP / P contribute together, inversely to other four classifiers
# PC4 [-0.3794475   0.77255619  0.11524398 -0.42995786  0.12742098 -0.21164945]     # B / LB / LP contribute together with B as a major driver, inversely to other 3 classifiers
# PC5 [ 0.08916003 -0.29989775  0.79908674 -0.47638702 -0.18814181  0.03507524]
# PC6 [-0.62302019 -0.02147484  0.02138479  0.08193396 -0.23980697  0.73939625]

# Perform PCA in a pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline

# Define a pipe object of the Pipeline class:
pipe = Pipeline( [
               ('scaler', StandardScaler() ),
               ('reducer', PCA() ) ] )

# Fit to the defined pipe object
pc = pipe.fit_transform(numeric_alleles_df)
# View the fit
print(pc[ :, :2 ])
# Initialize a new pandas dataframe to capture the PC categories
pc_categories = pd.DataFrame()
# Populate the pandas dataframe with all six principal components:
pc_categories['PC1'] = pc[:,0]
pc_categories['PC2'] = pc[:,1]
pc_categories['PC3'] = pc[:,2]
pc_categories['PC4'] = pc[:,3]
pc_categories['PC5'] = pc[:,4]
pc_categories['PC6'] = pc[:,5]

############## GET ALL PAIRWISE PCA PC IDENTITY PLOTS    as PNGs           # CONSIDER A LOOP STRUCTURE
sns.scatterplot(data=pc_categories,y='PC1',x='PC1',alpha=0.4)
# plt.show()
plt.savefig('PC1_PC1_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC2',x='PC2',alpha=0.4)
# plt.show()
plt.savefig('PC2_PC2_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC3',x='PC3',alpha=0.4)
# plt.show()
plt.savefig('PC3_PC3_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC4',x='PC4',alpha=0.4)
# plt.show()
plt.savefig('PC4_PC4_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC5',x='PC5',alpha=0.4)
# plt.show()
plt.savefig('PC5_PC5_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC6',x='PC6',alpha=0.4)
# plt.show()
plt.savefig('PC6_PC6_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############## GET ALL PAIRWISE PCA PC IDENTITY PLOTS    as SVGs           # CONSIDER A LOOP STRUCTURE
sns.scatterplot(data=pc_categories,y='PC1',x='PC1',alpha=0.4)
# plt.show()
plt.savefig('PC1_PC1_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC2',x='PC2',alpha=0.4)
# plt.show()
plt.savefig('PC2_PC2_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC3',x='PC3',alpha=0.4)
# plt.show()
plt.savefig('PC3_PC3_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC4',x='PC4',alpha=0.4)
# plt.show()
plt.savefig('PC4_PC4_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC5',x='PC5',alpha=0.4)
# plt.show()
plt.savefig('PC5_PC5_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC6',x='PC6',alpha=0.4)
# plt.show()
plt.savefig('PC6_PC6_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############## GET ALL PAIRWISE PCA PC vs. PC GRAPHS                # CONSIDER A LOOP STRUCTURE

############## PNG filetype
############## PC1-GROUP
sns.scatterplot(data=pc_categories,y='PC1',x='PC2',alpha=0.4)
# plt.show()
plt.savefig('PC1_PC2_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC1',x='PC3',alpha=0.4)
# plt.show()
plt.savefig('PC1_PC3_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC1',x='PC4',alpha=0.4)
# plt.show()
plt.savefig('PC1_PC4_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC1',x='PC5',alpha=0.4)
# plt.show()
plt.savefig('PC1_PC5_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC1',x='PC6',alpha=0.4)
# plt.show()
plt.savefig('PC1_PC6_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############## PC2-GROUP
sns.scatterplot(data=pc_categories,y='PC2',x='PC3',alpha=0.4)
# plt.show()
plt.savefig('PC2_PC3_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC2',x='PC4',alpha=0.4)
# plt.show()
plt.savefig('PC2_PC4_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC2',x='PC5',alpha=0.4)
# plt.show()
plt.savefig('PC2_PC5_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC2',x='PC6',alpha=0.4)
# plt.show()
plt.savefig('PC2_PC6_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############## PC3-GROUP
sns.scatterplot(data=pc_categories,y='PC3',x='PC4',alpha=0.4)
# plt.show()
plt.savefig('PC3_PC4_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC3',x='PC5',alpha=0.4)
# plt.show()
plt.savefig('PC3_PC5_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC3',x='PC6',alpha=0.4)
# plt.show()
plt.savefig('PC3_PC6_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############## PC4-GROUP
sns.scatterplot(data=pc_categories,y='PC4',x='PC5',alpha=0.4)
# plt.show()
plt.savefig('PC4_PC5_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC4',x='PC6',alpha=0.4)
# plt.show()
plt.savefig('PC4_PC6_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############## PC5-GROUP
sns.scatterplot(data=pc_categories,y='PC5',x='PC6',alpha=0.4)
# plt.show()
plt.savefig('PC5_PC6_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############## SVG filetype
############## PC1-GROUP
sns.scatterplot(data=pc_categories,y='PC1',x='PC2',alpha=0.4)
# plt.show()
plt.savefig('PC1_PC2_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC1',x='PC3',alpha=0.4)
# plt.show()
plt.savefig('PC1_PC3_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC1',x='PC4',alpha=0.4)
# plt.show()
plt.savefig('PC1_PC4_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC1',x='PC5',alpha=0.4)
# plt.show()
plt.savefig('PC1_PC5_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC1',x='PC6',alpha=0.4)
# plt.show()
plt.savefig('PC1_PC6_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############## PC2-GROUP
sns.scatterplot(data=pc_categories,y='PC2',x='PC3',alpha=0.4)
# plt.show()
plt.savefig('PC2_PC3_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC2',x='PC4',alpha=0.4)
# plt.show()
plt.savefig('PC2_PC4_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC2',x='PC5',alpha=0.4)
# plt.show()
plt.savefig('PC2_PC5_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC2',x='PC6',alpha=0.4)
# plt.show()
plt.savefig('PC2_PC6_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############## PC3-GROUP
sns.scatterplot(data=pc_categories,y='PC3',x='PC4',alpha=0.4)
# plt.show()
plt.savefig('PC3_PC4_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC3',x='PC5',alpha=0.4)
# plt.show()
plt.savefig('PC3_PC5_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC3',x='PC6',alpha=0.4)
# plt.show()
plt.savefig('PC3_PC6_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############## PC4-GROUP
sns.scatterplot(data=pc_categories,y='PC4',x='PC5',alpha=0.4)
# plt.show()
plt.savefig('PC4_PC5_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC4',x='PC6',alpha=0.4)
# plt.show()
plt.savefig('PC4_PC6_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############## PC5-GROUP
sns.scatterplot(data=pc_categories,y='PC5',x='PC6',alpha=0.4)
# plt.show()
plt.savefig('PC5_PC6_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

########################################################
# make the lower triangle plots (for a PC1, PC2, PC3 square 3x3 grid)
########################################################
#### SVGs output ####
sns.scatterplot(data=pc_categories,y='PC2',x='PC1',alpha=0.4)
# plt.show()
plt.savefig('PC2_PC1_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC3',x='PC1',alpha=0.4)
# plt.show()
plt.savefig('PC3_PC1_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC3',x='PC2',alpha=0.4)
# plt.show()
plt.savefig('PC3_PC2_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
#### PNGs output ####
sns.scatterplot(data=pc_categories,y='PC2',x='PC1',alpha=0.4)
# plt.show()
plt.savefig('PC2_PC1_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC3',x='PC1',alpha=0.4)
# plt.show()
plt.savefig('PC3_PC1_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC3',x='PC2',alpha=0.4)
# plt.show()
plt.savefig('PC3_PC2_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()


####################################################################################################
### Execute a k-means / PCA pipe via python's pipeline library
####################################################################################################

import tarfile
import urllib

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import LabelEncoder, MinMaxScaler

# Copy the raw allele counts to a new pandas dataframe called data
data = numeric_alleles_df

# Set the preprocessor parameters for PCA operations
# Set PCA number of components (nc) to 2, since the first two PCs explain 88.8% of the variance in the anemia allele dataset
preprocessor = Pipeline(
    [
        ("scaler", MinMaxScaler()),
        ("pca", PCA(n_components=2)),
    ]
)

# Set the clusterer parameters for k-means operations
## Setting cluster number to match TDA node number to align 
clusterer = Pipeline(
   [
       (
           "kmeans",
           KMeans(
               n_clusters=6,
               init="k-means++",
               n_init=50,
               max_iter=500,
           ),
       ),
   ]
)

# Set the pipe strategy, which sequentially triggers preprocessor followed by clusterer.  This effectively executes 2-component PCA followed by k-means clustering.
pipe = Pipeline(
    [
        ("preprocessor", preprocessor),
        ("clusterer", clusterer)
    ]
)

# Calling .fit() with data as the argument performs all the pipeline steps on the data:
pipe.fit(data)

# Apply the pipe to preprocess the data using only PCA:
preprocessed_data = pipe["preprocessor"].transform(data)
# Apply the pipe's kmeans clusterer to predict labels:
predicted_labels = pipe["clusterer"]["kmeans"].labels_

# Evaluate the performance by calculating the silhouette coefficient:
silhouette_score(preprocessed_data, predicted_labels)
# np.float64(0.5723888163141256)


####################################################################################################
## Prepare 2-PC PCA / k-means plotting environment
####################################################################################################

# Create a pca dataframe using pandas, fill it with pipe's PCA PC1 and PCA PC2 return values:
pcadf = pd.DataFrame(pipe["preprocessor"].transform(data),columns=["PC1", "PC2"],)
# Add a column to store the kmeans clusterer's predictions
pcadf["predicted_cluster"] = pipe["clusterer"]["kmeans"].labels_


# Visualize with a scatterplot (save as PNG)
## Set plot attrs
plt.style.use("fivethirtyeight")
# custom_palette = sns.color_palette("tab10")
# cmap = sns.diverging_palette ( h_neg=10, h_pos=240, as_cmap=True)
# sns.set_style("white", {"axes.facecolor": "white"})
sns.set_style("white")
# Fill the scat variable with Seaborn scatterplot elements:
scat = sns.scatterplot(x="PC1",y="PC2",data=pcadf,hue="predicted_cluster",palette="tab10",)
# scat = sns.scatterplot(x="PC1",y="PC2",data=pcadf,hue="predicted_cluster",palette=custom_palette,cmap=cmap,)
# Add a TITLE element to scat (sns var)
scat.set_title("k-means cluster: ~185K alleles from 164 anemia-enriched loci")
# Configure LEGEND element
handles, labels  =  scat.get_legend_handles_labels()
scat.legend(handles, labels, loc='center right')
# Display plot
# plt.show()
# Save plot as a PNG file to disk
plt.savefig('kmeans_clustered_PC2_PC1_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

# Visualize with a scatterplot (save as SVG)
## Set plot attrs
plt.style.use("fivethirtyeight")
# custom_palette = sns.color_palette("tab10")
# cmap = sns.diverging_palette ( h_neg=10, h_pos=240, as_cmap=True)
# sns.set_style("white", {"axes.facecolor": "white"})
sns.set_style("white")
# Fill the scat variable with Seaborn scatterplot elements:
scat = sns.scatterplot(x="PC1",y="PC2",data=pcadf,hue="predicted_cluster",palette="tab10",)
# scat = sns.scatterplot(x="PC1",y="PC2",data=pcadf,hue="predicted_cluster",palette=custom_palette,cmap=cmap,)
# Add a TITLE element to scat (sns var)
scat.set_title("k-means cluster: ~1.1 million alleles from 2359 loci")
# Configure LEGEND element
handles, labels  =  scat.get_legend_handles_labels()
scat.legend(handles, labels, loc='center right')
# Display plot
# plt.show()
# Save plot as a PNG file to disk
plt.savefig('kmeans_clustered_PC2_PC1_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()