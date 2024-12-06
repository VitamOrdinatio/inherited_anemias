## The goal of this script is to to generate SVG circles that are of relative proportions for the purposes of TDA node visualization.

## In its current implementation, the script simply performs a log2-transformation of the geneNumber size of each TDA node
## These log2(geneNumber) values are then manually shunted to BioRender to draw approximate circles with relative size differences.


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
cd /mnt/c/wslshare/github/inherited_anemias/step11_ClinVar_TDA
python

############################################################ IMPORTS
import pandas as pd
import numpy as np
import math
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import math

# Set the total node number (not index number).  Get this by manually opening the kmapper HTML output, and counting the nodes.
numNodes = 6
# Set each nodes' number of genes (guild membership).  Get this by manually examining the kmapper HTML output, and recording the guild membership for each node.  Guilds are TDA clusters, and guild members are genes.
numGenesPerNode = [107, 59, 36, 25, 14, 4]

# Calculate each nodes' log2-transformed geneNumber:
i = 0
log2numGenesPerNode = []

for i in range(len(numGenesPerNode)):
    # Fetch a node's gene number:
    geneNUM = numGenesPerNode[i]
    result = math.log2(geneNUM)
    log2numGenesPerNode.append(result)


# View log2 transformed geneNumber per node:
log2numGenesPerNode
# [6.741466986401147, 5.882643049361842, 5.169925001442312, 4.643856189774724, 3.807354922057604, 2.0]


newdf = pd.DataFrame()
newdf['nodeID'] = ['n1','n2','n3','n4','n5','n6']
newdf['numGenes'] = numGenesPerNode
newdf['log2(numGenes)'] = log2numGenesPerNode

newdf = newdf.set_index('nodeID')
newdf.to_csv('./out/anemia_199loci_ClinVar_4cats_TDA_nodeSizes.csv')