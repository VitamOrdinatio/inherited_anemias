## The following python script performs multivariate statistics (MVS) when given a ClinVar raw allele counts table as a CSV file in this format:
    # n rows where n = total gene number
    # 5 columns with labels: gene_name   B   LB   LP    P
        # B = benign
        # LB = likely benign
        # LP = likely pathogenic
        # P = pathogenic
    # NOTA BENE: The CSV file must contain the original allele counts (that are not normalized frequencies).

# The output of this script follows this general MVS sequence.  PNG/SVG plots are written to working directory:
    # 1. Seaborn pairplot analyses of all 4 ClinVar categories graphed against one another
    # 2. t-SNE dimensional reduction
    # 3. Linear Regression / correlation coefficient analyses
    # 4. PCA
    # 5. PCA explained variance ratio
    # 6. PCA explained variance ratios as accumulated SUMS
    # 7. PCA component ruleset definitions
    # 8. k-means clustering
    # 9. 2-PC_k-means_clustering
        

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
cd /mnt/c/wslshare/github/inherited_anemias/step07_getMultivariateStats
python

############################################################ IMPORTS
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE

# LOAD DATA
alleles_df = pd.DataFrame()                             # Note: alleles_df will never be sorted
alleles_df = pd.read_csv('./in/Anemia_199_genes_AlleleCount_4_ClinVar_classifiers.csv')
# Pull out numeric only columns (remove the 'gene_name' column)
non_numeric = ['gene_name']
numeric_alleles_df = alleles_df.drop(non_numeric, axis=1)

###########################################################


alleles_df.shape
# (199, 5)                    # returns a tuple which means the df has 199 rows and 5 columns

alleles_df.head()              # view just the first five genes of 199 genes
  # gene_name    B   LB   LP    P
# 0     ABCA1  250  580   18   71
# 1     ABCG8   81  202   24   45
# 2    ACVRL1   94  191  168  373
# 3      ADA2   44  150   31  121
# 4      ALG8   71   86   20   32

alleles_df.describe()
# returns simple stats analysis:
                # B           LB          LP            P
# count  199.000000   199.000000  199.000000   199.000000
# mean    62.055276   289.552764   54.994975   158.894472
# std     97.449712   481.069102   94.265770   474.646486
# min      0.000000     1.000000    0.000000     1.000000
# 25%     18.000000    30.500000    5.000000    18.000000
# 50%     38.000000   130.000000   18.000000    48.000000
# 75%     69.000000   360.000000   62.500000   139.000000
# max    845.000000  3752.000000  638.000000  5005.000000

alleles_df.loc[0]   # gets first of 199 genes
# gene_name    ABCA1
# B              250
# LB             580
# LP              18
# P               71
# Name: 0, dtype: object

alleles_df.loc[1]   # gets second of 199 genes
# gene_name    ABCG8
# B               81
# LB             202
# LP              24
# P               45
# Name: 1, dtype: object

alleles_df.loc[1][1]
# np.int64(81)                  # 81 alleles of this categorical type
alleles_df.loc[1][2]
# np.int64(202)                  # 202 alleles of this categorical type
alleles_df.loc[1][3]
# np.int64(24)                 # 24 alleles of this categorical type
alleles_df.loc[1][4]
# np.int64(45)                 # 45 alleles of this categorical type


############################################################ NUMERIC DF INITALIZIATION
# Pull out numeric only columns (remove the 'gene_name' column)
# non_numeric = ['gene_name']
# numeric_alleles_df = alleles_df.drop(non_numeric, axis=1)

numeric_alleles_df.shape
# (199, 4)
numeric_alleles_df.head()
     # B   LB   LP    P
# 0  250  580   18   71
# 1   81  202   24   45
# 2   94  191  168  373
# 3   44  150   31  121
# 4   71   86   20   32

numeric_alleles_df.describe()
                # B           LB          LP            P
# count  199.000000   199.000000  199.000000   199.000000
# mean    62.055276   289.552764   54.994975   158.894472
# std     97.449712   481.069102   94.265770   474.646486
# min      0.000000     1.000000    0.000000     1.000000
# 25%     18.000000    30.500000    5.000000    18.000000
# 50%     38.000000   130.000000   18.000000    48.000000
# 75%     69.000000   360.000000   62.500000   139.000000
# max    845.000000  3752.000000  638.000000  5005.000000

############################################################ SNS PAIRPLOT
## Run seaborn's pairplot analysis
# sns.pairplot(numeric_alleles_df, diag_kind ='hist')
## Make sure on Win11/WSL2 to first run XLaunch, and ensure that 'disable access controls' is toggled ON in last menu dialog
# plt.show()
## View the plot from memory

# Run seaborn's pairplot analysis

# get an SVG version
sns.pairplot(numeric_alleles_df, diag_kind ='hist')
# Make sure on Win11/WSL2 to first run XLaunch, and ensure that 'disable access controls' is toggled ON in last menu dialog
# SAVE the pairplot to file system
plt.savefig('./out/raw_pairplot.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

# get a PNG version
sns.pairplot(numeric_alleles_df, diag_kind ='hist')
plt.savefig('./out/raw_pairplot.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############################################################ TSNE
# t-SNE
m = TSNE(learning_rate=50)
tsne_features = m.fit_transform(numeric_alleles_df)
len(tsne_features)
# 199
tsne_features[1:4,:]
# array([[ 0.34670743, -3.937564  ],
       # [-1.6433936 , -9.3039665 ],
       # [ 2.637077  , -2.3393624 ]], dtype=float32)
# Initialize tsne df
tsne_df = pd.DataFrame()
# Transfer tsne output into tsne df
tsne_df['x'] = tsne_features[:,0]
tsne_df['y'] = tsne_features[:,1]

# Setup an sns scatterplot to visualize tsne analysis

# get an SVG version
sns.scatterplot( x="x", y="y", data=tsne_df )
# sns.scatterplot( x="x", y="y", hue="y", palette=sns.color_palette("hls", 10), data=tsne_df, legend="full", alpha=0.3 )
plt.savefig('./out/TSNE_scatterplot.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

# get a PNG version
sns.scatterplot( x="x", y="y", data=tsne_df )
plt.savefig('./out/TSNE_scatterplot.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

tsne_df.var()
# x     11.904635
# y    110.879532
# dtype: float32

# Boxplots

# get an SVG version
sns.boxplot(tsne_df)
# plt.show()
plt.savefig('./out/TSNE_boxplot.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

# get a PNG version
sns.boxplot(tsne_df)
plt.savefig('./out/TSNE_boxplot.png')
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
           # B        LB        LP         P
# B   1.000000  0.792313  0.673274  0.834417
# LB  0.792313  1.000000  0.797543  0.762214
# LP  0.673274  0.797543  1.000000  0.756673
# P   0.834417  0.762214  0.756673  1.000000

## Heatmaps of linear regressions

# get an SVG version
cmap = sns.diverging_palette ( h_neg=10, h_pos=240, as_cmap=True)
sns.heatmap(numeric_alleles_df.corr(), center=0, cmap=cmap, linewidths=1, annot=True, fmt=".2f")
plt.savefig('./out/Correlation_heatmap.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

# get a PNG version
cmap = sns.diverging_palette ( h_neg=10, h_pos=240, as_cmap=True)
sns.heatmap(numeric_alleles_df.corr(), center=0, cmap=cmap, linewidths=1, annot=True, fmt=".2f")
plt.savefig('./out/Correlation_heatmap.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

# create a MASK to remove the upper triangle

# get an SVG version
corr = numeric_alleles_df.corr()
mask = np.triu( np.ones_like(corr, dtype = bool) )
cmap = sns.diverging_palette ( h_neg=10, h_pos=240, as_cmap=True)
sns.heatmap(numeric_alleles_df.corr(), mask=mask, center=0, cmap=cmap, linewidths=1, annot=True, fmt=".2f")
plt.savefig('./out/Correlation_heatmap_LT.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

# get a PNG version
corr = numeric_alleles_df.corr()
mask = np.triu( np.ones_like(corr, dtype = bool) )
cmap = sns.diverging_palette ( h_neg=10, h_pos=240, as_cmap=True)
sns.heatmap(numeric_alleles_df.corr(), mask=mask, center=0, cmap=cmap, linewidths=1, annot=True, fmt=".2f")
plt.savefig('./out/Correlation_heatmap_LT.png')
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
          # B        LB        LP         P
# 0  1.933497  0.605276 -0.393444 -0.185646
# 1  0.194895 -0.182455 -0.329633 -0.240562
# 2  0.328634 -0.205379  1.201815  0.452222
# 3 -0.185745 -0.290820 -0.255188 -0.080039
# 4  0.092020 -0.424193 -0.372174 -0.268019

## compare with original (pre-fit/transform), which represent raw allele counts:
numeric_alleles_df.head()
     # B   LB   LP    P
# 0  250  580   18   71
# 1   81  202   24   45
# 2   94  191  168  373
# 3   44  150   31  121
# 4   71   86   20   32

# Bring in PCA tools:
from sklearn.decomposition import PCA

# Declare a pca object of the PCA class
pca = PCA()

# Fit the standardized dataset to perform PCA:
pca.fit(numeric_alleles_df_std)

pcaDF = pd.DataFrame()
pcaDF['PC'] = ['PC1','PC2','PC3','PC4']

# Get the PCA explained variance ratio
print(pca.explained_variance_ratio_)
# This tells me individual PC contributions to the overall PCA.  Usually the first 2-3 components contain most of dataset variance
# [0.82726659 0.0863436  0.0555871  0.03080271]
#   PC1         PC2         PC3         PC4

i = 0
ExpVarianceRatios = []      # float list

for i in range(len(pca.explained_variance_ratio_)):
    ExpVarianceRatios.append(float(pca.explained_variance_ratio_[i]))

pcaDF['ExpVarRatios'] = ExpVarianceRatios
pcaDF = pcaDF.set_index('PC')
pcaDF.to_csv('./out/PCA_explained_variance_ratio.csv')

        #   Gene name
        #   CATEGORY: Benign                (dark blue)     #29386F
        #   CATEGORY: Likely benign         (light blue)    #DFEDFA
        #   CATEGORY: Likely pathogenic     (light red)     #FCE5EA
        #   CATEGORY: Pathogenic            (dark red)      #802A2A
        


## Barplot (PNG)
ax = sns.barplot(data=pcaDF, x='PC', y='ExpVarRatios', color='#29386F')
ax.set_xticklabels(ax.get_xticklabels(), ha="right")
plt.xlabel('principal component')
plt.ylabel('explained variance ratio')
plt.xticks(fontsize=8, rotation=45)  # Rotate labels by 45 degrees, set font to 8
plt.yticks(fontsize=8)
plt.subplots_adjust(bottom=0.2, left=0.114, right=0.945, top=0.802)  # Adjust the value as needed
plt.savefig('./out/PCA_expVARIANCEratio.png')
# plt.show()
## Clear plot space
plt.clf()
plt.cla()
plt.close()

## Barplot (SVG)
ax = sns.barplot(data=pcaDF, x='PC', y='ExpVarRatios', color='#29386F')
ax.set_xticklabels(ax.get_xticklabels(), ha="right")
plt.xlabel('principal component')
plt.ylabel('explained variance ratio')
plt.xticks(fontsize=8, rotation=45)  # Rotate labels by 45 degrees, set font to 8
plt.yticks(fontsize=8)
plt.subplots_adjust(bottom=0.2, left=0.114, right=0.945, top=0.802)  # Adjust the value as needed
plt.savefig('./out/PCA_expVARIANCEratio.svg')
# plt.show()
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############## work on PCA explained variance ratio CUMULATIVE SUMS:

# Get the accumulated PCA explaind variance ratio
print(pca.explained_variance_ratio_.cumsum())
# [0.8018952  0.89285166 0.94310647 0.97612715 0.99475924 1.        ]
## This tells me that 94.2% of the data is explained by the first 3 PCs of PCA

    
pcaDF = pd.DataFrame()
pcaDF['PC_range'] = ['PC1','PC1-PC2','PC1-PC3','PC1-PC4']

# Get the PCA explained variance ratio
print(pca.explained_variance_ratio_.cumsum())
# This tells me individual PC contributions to the overall PCA.  Usually the first 2-3 components contain most of dataset variance
# [0.82726659 0.0863436  0.0555871  0.03080271]
#   PC1         PC2         PC3         PC4

i = 0
CumSumExpVarianceRatios = []      # float list

for i in range(len(pca.explained_variance_ratio_.cumsum())):
    CumSumExpVarianceRatios.append(float(pca.explained_variance_ratio_.cumsum()[i]))

pcaDF['ExpVarRatioCumSum'] = CumSumExpVarianceRatios
pcaDF = pcaDF.set_index('PC_range')
pcaDF.to_csv('./out/PCA_explained_variance_ratio_CumSum.csv')

        #   Gene name
        #   CATEGORY: Benign                (dark blue)     #29386F
        #   CATEGORY: Likely benign         (light blue)    #DFEDFA
        #   CATEGORY: Likely pathogenic     (light red)     #FCE5EA
        #   CATEGORY: Pathogenic            (dark red)      #802A2A
        


## Barplot (PNG)
ax = sns.barplot(data=pcaDF, x='PC_range', y='ExpVarRatioCumSum', color='#29386F')
ax.set_xticklabels(ax.get_xticklabels(), ha="right")
plt.xlabel('principal component')
plt.ylabel('explained variance ratio')
plt.xticks(fontsize=8, rotation=45)  # Rotate labels by 45 degrees, set font to 8
plt.yticks(fontsize=8)
plt.subplots_adjust(bottom=0.2, left=0.114, right=0.945, top=0.802)  # Adjust the value as needed
plt.savefig('./out/PCA_expVARIANCEratioCUMsum.png')
# plt.show()
## Clear plot space
plt.clf()
plt.cla()
plt.close()

## Barplot (SVG)
ax = sns.barplot(data=pcaDF, x='PC_range', y='ExpVarRatioCumSum', color='#29386F')
ax.set_xticklabels(ax.get_xticklabels(), ha="right")
plt.xlabel('principal component')
plt.ylabel('explained variance ratio')
plt.xticks(fontsize=8, rotation=45)  # Rotate labels by 45 degrees, set font to 8
plt.yticks(fontsize=8)
plt.subplots_adjust(bottom=0.2, left=0.114, right=0.945, top=0.802)  # Adjust the value as needed
plt.savefig('./out/PCA_expVARIANCEratioCUMsum.svg')
# plt.show()
## Clear plot space
plt.clf()
plt.cla()
plt.close()



############## work on PCA component ruleset:

pcaDF = pd.DataFrame()
pcaDF['PC_ruleset'] = ['PC1','PC2','PC3','PC4']

# Get the PCA component ruleset:
print(pca.components_)
## This lists all principal components:
#       B           LB          LP            P             # VARIANCE EXPLANATION BY PRINCIPAL COMPONENTS:
# [[ 0.499075    0.50677172  0.48674359  0.50713632]        # equal parts
 # [-0.59806355  0.20431099  0.71412089 -0.30101223]        # B changes with P
 # [ 0.19593629  0.70541386 -0.2962713  -0.61337066]        # LP changes with P
 # [ 0.59569548 -0.45147617  0.40661449 -0.52533877]]       # LB changes with P


i = 0
B_ComponentRuleSet = []      # floats
LB_ComponentRuleSet = []      # floats
LP_ComponentRuleSet = []      # floats
P_ComponentRuleSet = []      # floats

for i in range(len(pca.components_)):
    B_ComponentRuleSet.append(float(pca.components_[i][0]))
    LB_ComponentRuleSet.append(float(pca.components_[i][1]))
    LP_ComponentRuleSet.append(float(pca.components_[i][2]))
    P_ComponentRuleSet.append(float(pca.components_[i][3]))
    


pcaDF['B.rules'] = B_ComponentRuleSet
pcaDF['LB.rules'] = LB_ComponentRuleSet
pcaDF['LP.rules'] = LP_ComponentRuleSet
pcaDF['P.rules'] = P_ComponentRuleSet
pcaDF = pcaDF.set_index('PC_ruleset')
pcaDF.to_csv('./out/PCA_component_rulset.csv')




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
# Populate the pandas dataframe with all four principal components:
pc_categories['PC1'] = pc[:,0]
pc_categories['PC2'] = pc[:,1]
pc_categories['PC3'] = pc[:,2]
pc_categories['PC4'] = pc[:,3]

############## GET ALL PAIRWISE PCA PC PLOTS           # CONSIDER A LOOP STRUCTURE

# Get PNG

############## PC1-GROUP
sns.scatterplot(data=pc_categories,y='PC1',x='PC1',alpha=0.4)
# plt.show()
plt.savefig('./out/PC1_PC1_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC1',x='PC2',alpha=0.4)
# plt.show()
plt.savefig('./out/PC1_PC2_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC1',x='PC3',alpha=0.4)
# plt.show()
plt.savefig('./out/PC1_PC3_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC1',x='PC4',alpha=0.4)
# plt.show()
plt.savefig('./out/PC1_PC4_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############## PC2-GROUP
sns.scatterplot(data=pc_categories,y='PC2',x='PC1',alpha=0.4)
# plt.show()
plt.savefig('./out/PC2_PC1_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC2',x='PC2',alpha=0.4)
# plt.show()
plt.savefig('./out/PC2_PC2_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC2',x='PC3',alpha=0.4)
# plt.show()
plt.savefig('./out/PC2_PC3_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC2',x='PC4',alpha=0.4)
# plt.show()
plt.savefig('./out/PC2_PC4_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############## PC3-GROUP
sns.scatterplot(data=pc_categories,y='PC3',x='PC1',alpha=0.4)
# plt.show()
plt.savefig('./out/PC3_PC1_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC3',x='PC2',alpha=0.4)
# plt.show()
plt.savefig('./out/PC3_PC2_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC3',x='PC3',alpha=0.4)
# plt.show()
plt.savefig('./out/PC3_PC3_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC3',x='PC4',alpha=0.4)
# plt.show()
plt.savefig('./out/PC3_PC4_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############## PC4-GROUP
sns.scatterplot(data=pc_categories,y='PC4',x='PC1',alpha=0.4)
# plt.show()
plt.savefig('./out/PC4_PC1_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC4',x='PC2',alpha=0.4)
# plt.show()
plt.savefig('./out/PC4_PC2_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC4',x='PC3',alpha=0.4)
# plt.show()
plt.savefig('./out/PC4_PC3_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC4',x='PC4',alpha=0.4)
# plt.show()
plt.savefig('./out/PC4_PC4_scatter.png')
## Clear plot space
plt.clf()
plt.cla()
plt.close()



# Get SVG

############## PC1-GROUP
sns.scatterplot(data=pc_categories,y='PC1',x='PC1',alpha=0.4)
# plt.show()
plt.savefig('./out/PC1_PC1_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC1',x='PC2',alpha=0.4)
# plt.show()
plt.savefig('./out/PC1_PC2_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC1',x='PC3',alpha=0.4)
# plt.show()
plt.savefig('./out/PC1_PC3_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC1',x='PC4',alpha=0.4)
# plt.show()
plt.savefig('./out/PC1_PC4_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############## PC2-GROUP
sns.scatterplot(data=pc_categories,y='PC2',x='PC1',alpha=0.4)
# plt.show()
plt.savefig('./out/PC2_PC1_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC2',x='PC2',alpha=0.4)
# plt.show()
plt.savefig('./out/PC2_PC2_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC2',x='PC3',alpha=0.4)
# plt.show()
plt.savefig('./out/PC2_PC3_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC2',x='PC4',alpha=0.4)
# plt.show()
plt.savefig('./out/PC2_PC4_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############## PC3-GROUP
sns.scatterplot(data=pc_categories,y='PC3',x='PC1',alpha=0.4)
# plt.show()
plt.savefig('./out/PC3_PC1_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC3',x='PC2',alpha=0.4)
# plt.show()
plt.savefig('./out/PC3_PC2_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC3',x='PC3',alpha=0.4)
# plt.show()
plt.savefig('./out/PC3_PC3_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC3',x='PC4',alpha=0.4)
# plt.show()
plt.savefig('./out/PC3_PC4_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()

############## PC4-GROUP
sns.scatterplot(data=pc_categories,y='PC4',x='PC1',alpha=0.4)
# plt.show()
plt.savefig('./out/PC4_PC1_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC4',x='PC2',alpha=0.4)
# plt.show()
plt.savefig('./out/PC4_PC2_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC4',x='PC3',alpha=0.4)
# plt.show()
plt.savefig('./out/PC4_PC3_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()
sns.scatterplot(data=pc_categories,y='PC4',x='PC4',alpha=0.4)
# plt.show()
plt.savefig('./out/PC4_PC4_scatter.svg')
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
# np.float64(0.5621971465820254)


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
scat.set_title("k-means cluster: ~112,534 alleles across 4 cats of 199 anemia-enriched loci")
# Configure LEGEND element
handles, labels  =  scat.get_legend_handles_labels()
scat.legend(handles, labels, loc='center right')
# Display plot
# plt.show()
# Save plot as a PNG file to disk
plt.savefig('./out/kmeans_clustered_PC2_PC1_scatter.png')
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
scat.set_title("k-means cluster: ~112,534 alleles across 4 cats of 199 anemia-enriched loci")
# Configure LEGEND element
handles, labels  =  scat.get_legend_handles_labels()
scat.legend(handles, labels, loc='center right')
# Display plot
# plt.show()
# Save plot as a SVG file to disk
plt.savefig('./out/kmeans_clustered_PC2_PC1_scatter.svg')
## Clear plot space
plt.clf()
plt.cla()
plt.close()