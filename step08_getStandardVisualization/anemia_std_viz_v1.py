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
alleles_df = pd.read_csv('./in/Anemia_199_genes_AlleleCount_4_ClinVar_classifiers.csv')
# Extract columns into list variables:
genes = alleles_df['gene_name'].tolist() 
rawBalleles = alleles_df['B'].tolist() 
rawLBalleles = alleles_df['LB'].tolist() 
rawLPalleles = alleles_df['LP'].tolist() 
rawPalleles = alleles_df['P'].tolist() 

# Pull out numeric only columns (remove the 'gene_name' column)
# non_numeric = ['gene_name']
# numeric_alleles_df = alleles_df.drop(non_numeric, axis=1)


#11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
#111111111111111111111111111111111111111             step 1: Table of raw allele counts           1111111111111111111111111111111
#11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111

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
alleles_df.to_csv('./out/anemia_ClinVar_B_LB_LP_P_raw_allele_counts.csv')


#22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
#222222222222222222             step 2: Table of categorical allele frequencies (normalized per locus)          22222222222222222
#22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222

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
dfAlleleFreqs.to_csv('./out/anemia_ClinVar_B_LB_LP_P_categorical_allele_freqs.csv')


#33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
#3333333333333333333333         step 3: bar plots of log2-transformed total allele counts (B+LB+LP+P)          333333333333333333
#33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
    
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
plt.savefig('./out/log2barplot.png')
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
plt.savefig('./out/log2barplot.svg')
# plt.show()
## Clear plot space
plt.clf()
plt.cla()
plt.close()



#44444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444
#444444444444444           step 4: histogram plots of log2-transformed total allele counts (B+LB+LP+P)          44444444444444444
#44444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444

    # 4. Histogram (10% bin sizes) of total allele counts (B+LB+LP+P)

## Histogram (PNG)
ax = sns.histplot(data=log2df_sorted_byLOG2, x="log2_alleleSUM", bins = 10, color='#802A2A', alpha=0.92)
plt.xlabel('ten percentile bins of log2-transformed allele sums', fontsize=10)
plt.ylabel('gene count', fontsize=10)
plt.savefig('./out/log2histogram.png')
# plt.show()
plt.clf()
plt.cla()
plt.close()

## Histogram (SVG)
ax = sns.histplot(data=log2df_sorted_byLOG2, x="log2_alleleSUM", bins = 10, color='#802A2A', alpha=0.92)
plt.xlabel('ten percentile bins of log2-transformed allele sums', fontsize=10)
plt.ylabel('gene count', fontsize=10)
plt.savefig('./out/log2histogram.svg')
# plt.show()
plt.clf()
plt.cla()
plt.close()

#55555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555
#5555555555555            step 5: PIE chart showing top20 genes for total allele counts (B+LB+LP+P)            555555555555555555
#55555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555555

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
plt.savefig('./out/pieTop20rawCounts.png')
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
plt.savefig('./out/pieTop20rawCounts.svg')
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


#66666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666
#6666666     step 6: Horizontal bar chart showing a sorted list of genes by top allele counts in a descending sort.     666666666
#66666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666

    # 6. Horizontal bar chart showing the sorted list of genes by top allele counts (descending sort).
        #   For each gene, raw number of each of 4 ClinVar categories: B, LB, LP, P.  Follow this color schema:
                #   CATEGORY: Benign                (dark blue)     #29386F
                #   CATEGORY: Likely benign         (light blue)    #DFEDFA
                #   CATEGORY: Likely pathogenic     (light red)     #FCE5EA
                #   CATEGORY: Pathogenic            (dark red)      #802A2A     


## Pseudocode for making this horizonatal bar chart of raw allele counts (4 cats)
    # 1) I want to plot, from left to right, P then LP then LB then B.
    # 2) To do this I want need to plot in this order:
            # B   as topAlleleSums
            # LB  as sum of LB+LP+P    OR as difference of SUM-B
            # LP  as sum of    LP+P    OR as difference of SUM-B-LB
            # P   as just P       P    OR as difference of SUM-B-LB-LP
    # 3) In doing so, I'm basically plotting bars on top of each other, and using layers to occlude data.
    # 4) The final result yields a horizontal stacked bar plot.
    

# Slice the top20 loci from the sorted log2df, and store in a pd variable
top20byAlleleCount = log2df_sorted_byLOG2[0:20]
topBvals = top20byAlleleCount['B'].tolist()
topLBvals = top20byAlleleCount['LB'].tolist()
topLPvals = top20byAlleleCount['LP'].tolist()
topPvals = top20byAlleleCount['P'].tolist()
topGeneNames = top20byAlleleCount.index        # index = gene names
topAlleleSums = top20byAlleleCount['alleleSUM'].tolist()

df = pd.DataFrame()
df['gene_name'] = topGeneNames
df['B'] = topBvals
df['LB'] = topLBvals
df['LP'] = topLPvals
df['P'] = topPvals
df['alleleSUM'] = topAlleleSums

i = 0
sum_LB_LP_P = []
sum_LP_P = []

for i in range(len(topGeneNames)):
    # Retrieve current genes raw allele counts per 1 of 4 ClinVar categories (B, LB, LP, and P):
    Bcount = int(df['B'][i])
    LBcount = int(df['LB'][i])
    LPcount = int(df['LP'][i])
    Pcount = int(df['P'][i])
    # Generate categorical sums needed for occlusion plotting methods:
    sum_LB_LP_P.append(LBcount + LPcount + Pcount)
    sum_LP_P.append(LPcount + Pcount)

# Write new columns to df pd variable (df represents the top20 genes by allele count sums across 4 cats):
df['sum_LB_LP_P'] = sum_LB_LP_P
df['sum____LP_P'] = sum_LP_P

## Plotting Sequence is CRITICAL to achieve accurate horizontal bar plot presentation since occlusion is key.
    # 2) To do this I want need to plot in this order:
            # B   as topAlleleSums                                              data = df['alleleSUM']      #29386F  
            # LB  as sum of LB+LP+P    OR as difference of SUM-B                data = df['sum_LB_LP_P']    #DFEDFA
            # LP  as sum of    LP+P    OR as difference of SUM-B-LB             data = df['sum____LP_P']    #FCE5EA
            # P   as just P       P    OR as difference of SUM-B-LB-LP          data = df['P']              #802A2A 

                #   CATEGORY: Benign                (dark blue)     #29386F
                #   CATEGORY: Likely benign         (light blue)    #DFEDFA
                #   CATEGORY: Likely pathogenic     (light red)     #FCE5EA
                #   CATEGORY: Pathogenic            (dark red)      #802A2A   


## I tested the following working tutorial on horizontal bar plots: # https://seaborn.pydata.org/examples/part_whole_bars.html

############################################ REFERENCE MATERIAL for step6 rawAlleleCount_top20genes (horiz barplot)

            # B   as topAlleleSums                                              data = df['alleleSUM']      #29386F  
            # LB  as sum of LB+LP+P    OR as difference of SUM-B                data = df['sum_LB_LP_P']    #DFEDFA
            # LP  as sum of    LP+P    OR as difference of SUM-B-LB             data = df['sum____LP_P']    #FCE5EA
            # P   as just P       P    OR as difference of SUM-B-LB-LP          data = df['P']              #802A2A 

#### Get a PNG of a rawAlleleCount_top20genes horizontal bar plot of 4 ClinVar categories
# Set gridspace
sns.set_theme(style="whitegrid")
# Initialize the matplotlib figure
f, ax = plt.subplots(figsize=(6, 8))
### Generate the 4 horizontal bar plots using OCCLUSION METHOD:
# BENIGN: Plot the total allele counts (which will later represent only B due to layer occlusion):
sns.set_color_codes("pastel")
sns.barplot(x="alleleSUM", y="gene_name", data=df,
            label="Benign", color="#29386F")
# LIKELY BENIGN: Plot the SUM of LB+LP+P (which will later represent only LB due to layer occlusion):
sns.set_color_codes("muted")
sns.barplot(x="sum_LB_LP_P", y="gene_name", data=df,
            label="Likely Benign", color="#DFEDFA")
# LIKELY PATHOGENIC: Plot the SUM of LP+P (which will later represent only Lp due to layer occlusion):
sns.set_color_codes("muted")
sns.barplot(x="sum____LP_P", y="gene_name", data=df,
            label="Likely Pathogenic", color="#FCE5EA")
# PATHOGENIC: Plot the P allele counts (which will later represent only P due to layer occlusion):
sns.set_color_codes("pastel")
sns.barplot(x="P", y="gene_name", data=df,
            label="Pathogenic", color="#802A2A")
# Set legend title
ax.legend(title="ClinVar Allele Category")
# Set axis labels
ax.set(ylabel='Gene Name', xlabel='Total Allele Counts')
# Save plot to disk
plt.savefig('./out/hBARtop20genes_byAlleleCounts.png')
# Show plot (debug state)
# plt.show()
# Clear plot space
plt.clf()
plt.cla()
plt.close()

#### Get an SVG of a rawAlleleCount_top20genes horizontal bar plot of 4 ClinVar categories
# Set gridspace
sns.set_theme(style="whitegrid")
# Initialize the matplotlib figure
f, ax = plt.subplots(figsize=(6, 8))
### Generate the 4 horizontal bar plots using OCCLUSION METHOD:
# BENIGN: Plot the total allele counts (which will later represent only B due to layer occlusion):
sns.set_color_codes("pastel")
sns.barplot(x="alleleSUM", y="gene_name", data=df,
            label="Benign", color="#29386F")
# LIKELY BENIGN: Plot the SUM of LB+LP+P (which will later represent only LB due to layer occlusion):
sns.set_color_codes("muted")
sns.barplot(x="sum_LB_LP_P", y="gene_name", data=df,
            label="Likely Benign", color="#DFEDFA")
# LIKELY PATHOGENIC: Plot the SUM of LP+P (which will later represent only Lp due to layer occlusion):
sns.set_color_codes("muted")
sns.barplot(x="sum____LP_P", y="gene_name", data=df,
            label="Likely Pathogenic", color="#FCE5EA")
# PATHOGENIC: Plot the P allele counts (which will later represent only P due to layer occlusion):
sns.set_color_codes("pastel")
sns.barplot(x="P", y="gene_name", data=df,
            label="Pathogenic", color="#802A2A")
# Set legend title
ax.legend(title="ClinVar Allele Category")
# Set axis labels
ax.set(ylabel='Gene Name', xlabel='Total Allele Counts')
# Save plot to disk
plt.savefig('./out/hBARtop20genes_byAlleleCounts.svg')
# Show plot (debug state)
# plt.show()
# Clear plot space
plt.clf()
plt.cla()
plt.close()

#77777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777
#777777777777777777777777777          step 7: PIE chart showing top20 genes for allele freqs (LP+P)       77777777777777777777777
#77777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777777
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
plt.savefig('./out/pieTop20normFreqs.png')
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
plt.savefig('./out/pieTop20normFreqs.svg')
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




#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
#8888888    step 8: Horizontal bar chart of top problematic allele categorical freqs (P over LP, over LB, over B)        88888888
#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888  

    # 6. Horizontal bar chart showing the sorted list of genes by top problematic allele frequencies (descending sort).
        #   For each gene, raw number of each of 4 ClinVar categories: B, LB, LP, P.  Follow this color schema:
                #   CATEGORY: Benign                (dark blue)     #29386F
                #   CATEGORY: Likely benign         (light blue)    #DFEDFA
                #   CATEGORY: Likely pathogenic     (light red)     #FCE5EA
                #   CATEGORY: Pathogenic            (dark red)      #802A2A     


## Pseudocode for making this horizonatal bar chart of raw allele counts (4 cats)
    # 1) I want to plot, from left to right, B.freq then LB.freq then LP.freq then P.freq.
    # 2) To do this I want need to plot in reverse order to exploit OCCLUSION:
            #  P.freq  as topAlleleSums
            # LP.freq  as sum of LP+LB+B    OR as difference of SUM-P
            # LB.freq  as sum of    LB+B    OR as difference of SUM-P-LP
            #  B.freq  as just P       B    OR as difference of SUM-P-LP-LB
    # 3) In doing so, I'm basically plotting bars on top of each other, and using layers to occlude data.
    # 4) The final result yields a horizontal stacked bar plot.
    

# Copy the original dfAlleleFreqs dataframe, after sorting by P% then LP% then LB % then B%
top20byProblematicAlleleFrequency = dfAlleleFreqs.sort_values(by=['P.freq', 'LP.freq', 'LB.freq', 'B.freq'], ascending=False)
# Slice to get the first 20 genes, and reinitialize top20 alleleFreq dataframe:
top20byProblematicAlleleFrequency = top20byProblematicAlleleFrequency[0:20]
# Extract listvars of allleFreq categories for later computation (necessary for occlusion plots):
topBfreqs = top20byProblematicAlleleFrequency['B.freq'].tolist()
topLBfreqs = top20byProblematicAlleleFrequency['LB.freq'].tolist()
topLPfreqs = top20byProblematicAlleleFrequency['LP.freq'].tolist()
topPfreqs = top20byProblematicAlleleFrequency['P.freq'].tolist()
topGeneNames = top20byProblematicAlleleFrequency.index                              # index = gene names
topAlleleFreqSums = top20byProblematicAlleleFrequency['alleleFreqSUM'].tolist()

# Pass to df (pandas dataframe variable) that I will use in plotting (data=df)
df = pd.DataFrame()
df['gene_name'] = topGeneNames
df['B.freq'] = topBfreqs
df['LB.freq'] = topLBfreqs
df['LP.freq'] = topLPfreqs
df['P.freq'] = topPfreqs
df['alleleSUM.freq'] = topAlleleFreqSums


i = 0
freqsum_LP_LB_B = []
freqsum_LB_B = []

for i in range(len(topGeneNames)):
    # Retrieve current genes allele frequencies per 1 of 4 ClinVar categories (B, LB, LP, and P):
    Bfreq = int(df['B.freq'][i])
    LBfreq = int(df['LB.freq'][i])
    LPfreq = int(df['LP.freq'][i])
    Pfreq = int(df['P.freq'][i])
    # Generate categorical sums needed for occlusion plotting methods:
    freqsum_LP_LB_B.append(LPfreq + LBfreq + Bfreq)
    freqsum_LB_B.append(LBfreq + Bfreq)

# Write new columns to df pd variable (df represents the top20 genes by allele frequencies across 4 cats):
df['freqsum_LP_LB_B'] = freqsum_LP_LB_B
df['freqsum____LB_B'] = freqsum_LB_B

## Plotting Sequence is CRITICAL to achieve accurate horizontal bar plot presentation since occlusion is key.
    # 2) To do this I want need to plot in reverse order to exploit OCCLUSION:
            #  P.freq  as topAlleleSums                                         data = df['alleleSUM.freq']     #802A2A
            # LP.freq  as sum of LP+LB+B    OR as difference of SUM-P           data = df['freqsum_LP_LB_B']    #FCE5EA
            # LB.freq  as sum of    LB+B    OR as difference of SUM-P-LP        data = df['freqsum____LB_B']    #DFEDFA        
            #  B.freq  as just P       B    OR as difference of SUM-P-LP-LB     data = df['B.freq']             #29386F       
    # 3) In doing so, I'm basically plotting bars on top of each other, and using layers to occlude data.
    # 4) The final result yields a horizontal stacked bar plot.

                #   CATEGORY: Benign                (dark blue)     #29386F
                #   CATEGORY: Likely benign         (light blue)    #DFEDFA
                #   CATEGORY: Likely pathogenic     (light red)     #FCE5EA
                #   CATEGORY: Pathogenic            (dark red)      #802A2A   


## I tested the following working tutorial on horizontal bar plots: # https://seaborn.pydata.org/examples/part_whole_bars.html

## REFERENCE MATERIAL FOR STEP 8 (HBAR PLOTS ON PROBLEMATIC ALLELE FREQS)
## Plotting Sequence is CRITICAL to achieve accurate horizontal bar plot presentation since occlusion is key.
    # 2) To do this I want need to plot in this precise order to exploit OCCLUSION:
            #  P.freq  as topAlleleSums                                         data = df['alleleSUM.freq']     #802A2A
            # LP.freq  as sum of LP+LB+B    OR as difference of SUM-P           data = df['freqsum_LP_LB_B']    #FCE5EA
            # LB.freq  as sum of    LB+B    OR as difference of SUM-P-LP        data = df['freqsum____LB_B']    #DFEDFA        
            #  B.freq  as just P       B    OR as difference of SUM-P-LP-LB     data = df['B.freq']             #29386F   

#### Get a PNG of a rawAlleleCount_top20genes horizontal bar plot of 4 ClinVar categories
# Set gridspace
sns.set_theme(style="whitegrid")
# Initialize the matplotlib figure
f, ax = plt.subplots(figsize=(6, 8))
### Generate the 4 horizontal bar plots using OCCLUSION METHOD:
# PATHOGENIC: Plot the total allele frequencies (which will later represent only P due to layer occlusion):
sns.set_color_codes("pastel")
sns.barplot(x="alleleSUM.freq", y="gene_name", data=df,
            label="Pathogenic", color="#802A2A")
# LIKELY PATHOGENIC: Plot the SUM of LB+LP+P (which will later represent only LB due to layer occlusion):
sns.set_color_codes("muted")
sns.barplot(x="freqsum_LP_LB_B", y="gene_name", data=df,
            label="Likely Pathogenic", color="#FCE5EA")
# LIKELY BENIGN: Plot the SUM of LP+P (which will later represent only Lp due to layer occlusion):
sns.set_color_codes("muted")
sns.barplot(x="freqsum____LB_B", y="gene_name", data=df,
            label="Likely Benign", color="#DFEDFA")
# BENIGN: Plot the P allele counts (which will later represent only P due to layer occlusion):
sns.set_color_codes("pastel")
sns.barplot(x="B.freq", y="gene_name", data=df,
            label="Benign", color="#29386F")
plt.subplots_adjust(left=0.186)  # Adjust the value as needed
# Set legend title
# ax.legend(title="ClinVar Allele\nFrequency Category")
plt.legend(shadow=True, fancybox=True, ncol=1, title="Allele Frequency Category")
# Set axis labels
ax.set(ylabel='Gene Name', xlabel='Allele Frequency')
# Save plot to disk
plt.savefig('./out/hBARtop20genes_byAlleleFreqs.png')
# Show plot (debug state)
# plt.show()
# Clear plot space
plt.clf()
plt.cla()
plt.close()




#### Get a SVG of a rawAlleleCount_top20genes horizontal bar plot of 4 ClinVar categories
# Set gridspace
sns.set_theme(style="whitegrid")
# Initialize the matplotlib figure
f, ax = plt.subplots(figsize=(6, 8))
### Generate the 4 horizontal bar plots using OCCLUSION METHOD:
# PATHOGENIC: Plot the total allele frequencies (which will later represent only P due to layer occlusion):
sns.set_color_codes("pastel")
sns.barplot(x="alleleSUM.freq", y="gene_name", data=df,
            label="Pathogenic", color="#802A2A")
# LIKELY PATHOGENIC: Plot the SUM of LB+LP+P (which will later represent only LB due to layer occlusion):
sns.set_color_codes("muted")
sns.barplot(x="freqsum_LP_LB_B", y="gene_name", data=df,
            label="Likely Pathogenic", color="#FCE5EA")
# LIKELY BENIGN: Plot the SUM of LP+P (which will later represent only Lp due to layer occlusion):
sns.set_color_codes("muted")
sns.barplot(x="freqsum____LB_B", y="gene_name", data=df,
            label="Likely Benign", color="#DFEDFA")
# BENIGN: Plot the P allele counts (which will later represent only P due to layer occlusion):
sns.set_color_codes("pastel")
sns.barplot(x="B.freq", y="gene_name", data=df,
            label="Benign", color="#29386F")
plt.subplots_adjust(left=0.186)  # Adjust the value as needed
# Set legend title
# ax.legend(title="ClinVar Allele\nFrequency Category")
plt.legend(shadow=True, fancybox=True, ncol=1, title="Allele Frequency Category")
# Set axis labels
ax.set(ylabel='Gene Name', xlabel='Allele Frequency')
# Save plot to disk
plt.savefig('./out/hBARtop20genes_byAlleleFreqs.svg')
# Show plot (debug state)
# plt.show()
# Clear plot space
plt.clf()
plt.cla()
plt.close()