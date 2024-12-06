## The following python script performs multivariate statistics (MVS) when given a ClinVar raw allele counts table as a CSV file in this format:
    # n rows where n = total gene number
    # 7 columns with labels: gene_name   CC    B   LB   US   LP    P
        # CC = conflicting classifications
        # B = benign
        # LB = likely benign
        # US = uncertain signficance
        # LP = likely pathogenic
        # P = pathogenic
    # NOTA BENE: The CSV file must contain the original allele counts (that are not normalized frequencies).

# The output of this script follows this general MVS sequence.  PNG/SVG plots are written to working directory:
    # 1. Seaborn pairplot analyses of all 6 ClinVar categories graphed against one another
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

### After VENV vPCA construction, upon WSL bootup:
source gotovenv.sh
source vPCA/bin/activate
cd anemia/mitoCarta/
python

############################################################ IMPORTS
import pandas as pd
import numpy as np
# import seaborn as sns
# import matplotlib.pyplot as plt
# from sklearn.manifold import TSNE

# LOAD Anemia Loci
anemiaLOCI = pd.DataFrame()
anemiaLOCI = pd.read_csv('anemia_199genes_6ClinVar_classifiers.csv')

# Pull out list of anemia-enriched genes
anemiaGENES = anemiaLOCI['gene_name'].tolist()  

# Load MitoCarta Loci
mitocartaLOCI = pd.DataFrame()
mitocartaLOCI = pd.read_csv('Human_MitoCarta3.0_v1.csv')

# Pull out list of mitocarta gene names (as symbols)
mitocartaSymbols = mitocartaLOCI['Symbol'].tolist()
# Pull out list of mitocarta gene aliases (as synonyms)  
mitocartaSynonyms = mitocartaLOCI['Synonyms'].tolist()

# Create a single string variable in which you've appended all mitocarta symbols and synonyms
strMitoCartaGenes = ""
i = 0
# Collate all mitocarta gene info into a single string variable (that I can search against later)
for i in range(len(mitocartaSymbols)):
    strMitoCartaGenes = strMitoCartaGenes + "," + mitocartaSymbols[i] + "," + mitocartaSynonyms[i] 
    

len(strMitoCartaGenes)
# 29782  characters

# Convert strMitoCartaGenes to ALL-CAPS characters
strMitoCartaGenes = strMitoCartaGenes.upper()

# Set accumulator variable that tracks if gene is found in the MitoCarta list
geneIsInMitoCarta = []
j = 0
for j in range(len(anemiaGENES)):
    # Get current gene name
    gene = anemiaGENES[j]
    # Cross reference the MitoCarta list
    query = strMitoCartaGenes.find(gene)
    # Test query variable.  If -1, then gene was not found on MitoCarta list.
    if query == -1:
        geneIsInMitoCarta.append(0)
    elif query >= 0:
        geneIsInMitoCarta.append(1)


### Summary: There are 54 of 199 anemia-enriched genes that are found in the MitoCarta database!



###### Write Anemia Loci in MitoCarta status
# Write a new column to anemiaLOCI pd dataframe (has ClinVar data in it)
anemiaLOCI['inMitoCarta'] = geneIsInMitoCarta
# Set index to gene name:
anemiaLOCI = anemiaLOCI.set_index('gene_name')
# Write to disk
anemiaLOCI.to_csv('anemia_199genes_6ClinVar_MitoCartaStatus.csv')