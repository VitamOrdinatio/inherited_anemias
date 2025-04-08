## Visualization of the interaction between genes (i.e., list of 199 anemia-enriched loci) and phenotypes (i.e., HPO enrichment terms regarding clinical manifestations)

## The goal of the following python script is to generate two word cloud figures:
    # 1) a word cloud of how frequent the 199 anemia-enriched loci are observed across 619 HPO terms
    # 2) a word cloud of how frequent the 619 HPO terms are encountered in the 199 loci dataset
    

##### LEFT OFF HERE.....GOING TO TEST THIS DATACAMP TUTORIAL:
    # https://www.datacamp.com/tutorial/wordcloud-python

    
### 1SHOT: BASH code to setup VENV pocket for vWordCloud:
# python3 -m venv vWordCloud
# source vWordCloud/bin/activate
# python3 -m pip install --upgrade pip
# python3 -m pip install numpy
# python3 -m pip install pandas
# python3 -m pip install matplotlib # will also install numpy
# python3 -m pip install pillow
# # standard wordcloud version is missing masks
# python3 -m pip install wordcloud        
# # get the latest wordcloud version baked into your VENV pocket.  This latest version has mask features.
# # pip install -e git+https://github.com/amueller/word_cloud.git
# deactivate



### Activate VENV pocket
source gotovenv.sh
source vWordCloud/bin/activate
cd /mnt/c/wslshare/github/inherited_anemias/step10_vizHPOterms/
python



#### Start Python codeblock:

# Start with loading all necessary libraries
import numpy as np
import pandas as pd
from os import path
from PIL import Image
from wordcloud import WordCloud, STOPWORDS, ImageColorGenerator

import matplotlib.pyplot as plt
# % matplotlib inline

import warnings
warnings.filterwarnings("ignore")

import random

# LOAD DATA

# HPO enrichment for select genes:
HPOgenesDF = pd.DataFrame()                             
HPOgenesDF = pd.read_csv('./in/anemia_199_loci_HPO_TDAprep.csv')

# Extract a list of genes from loaded dataset:
genes = HPOgenesDF['gene']
# Convert pandas series to a list of string elements:
genes = genes.tolist()
len(genes)      # 199 loci

# Pull out numeric only columns (remove the 'gene' column)
non_numeric = ['gene']
numeric_HPOgenesDF = HPOgenesDF.drop(non_numeric, axis=1)
# Extract a list of HPO terms from loaded dataset:
HPOterms = numeric_HPOgenesDF.columns
# Convert pandas index to a list of string elements
HPOterms = HPOterms.tolist()
len(HPOterms)   # 619 terms

# Count HPO term enrichment per single genetic locus:
row_sums = numeric_HPOgenesDF.sum(axis=1)
# Convert pandas series to a list of int elements:
row_sums = row_sums.tolist()
# Measure row_sums (gene = row)
len(row_sums)
# 199

# Count the number of genes are tagged per single HPO term:
col_sums = numeric_HPOgenesDF.sum(axis=0)
# Convert pandas series to a list of int elements:
col_sums = col_sums.tolist()
# Measure col_sums (column = HPO term)
len(col_sums)
# 619

# ## Useful variables:
# ## dataframes
# HPOgenesDF
# numeric_HPOgenesDF
# ## listVars
# genes       # list of string elements, each string = gene name
# row_sums    # list of INT elements, each element = SUM of how often gene name appears across all HPO terms
# HPOterms    # list of string elements, each string = HPO term name
# col_sums    # list of INT elements, each element = SUM of how often an HPO term appears across all genetic loci


# First, generate an input list for gene names word cloud synthesis

# Flush counters
i = 0
j = 0
# Set accumulator list variable:
geneNameEssay = []

for i in range(len(genes)):
    # get a gene name (str) from list of genes:
    gene = genes[i]
    # get the sum of this gene's name as seen across all HPO terms:
    observations = row_sums[i]
    j = 0
    for j in range(observations):
        geneNameEssay.append(gene)


len(geneNameEssay)      # 23,625 string elements.  Each gene = str element.
type(geneNameEssay)     # list datatype

# Randomly shuffle the geneNameEssay
random.shuffle(geneNameEssay)
random.shuffle(geneNameEssay)

# Convert the randomized geneName listVar to a string
strGeneNameEssay = ', '.join(geneNameEssay)


### Generate WordClouds on GENE NAME FREQUENCIES:


# ## PNG file block:

# # Pass to a dedicated text variable:
# text = strGeneNameEssay
# # Create and generate a word cloud image:
# # Changing optional word cloud arguments
# wordcloud = WordCloud(max_font_size=50, max_words=40, colormap='plasma', background_color="white").generate(text)
# plt.figure()
# plt.imshow(wordcloud, interpolation="bilinear")
# plt.axis("off")
# plt.savefig('./out/wordCloud_geneFreq_HPOenrichment.png')
# # plt.show()


## SVG file block:

# Pass to a dedicated text variable:
text = strGeneNameEssay
# Create and generate a word cloud image:
# Changing optional word cloud arguments
wordcloud = WordCloud(width=1600, height=800, max_words=40, colormap='plasma', background_color='#F7F2FE').generate(text)
plt.figure(figsize=(20,10))
plt.imshow(wordcloud, interpolation="bilinear")
plt.axis("off")
plt.savefig('./out/wordCloud_geneFreq_HPOenrichment.svg')
# plt.show()



############################################################

## Next, generate an input list for HPO term word cloud synthesis
# Flush counters
i = 0
j = 0
# Set accumulator list variable:
hpoNameEssay = []

for i in range(len(HPOterms)):
    # get an HPO term (str) from list of HPO terms
    hpoTerm = HPOterms[i]
    # get the sum of this HPO term's usage as seen across all genes:
    observations = col_sums[i]
    j = 0
    for j in range(observations):
        hpoNameEssay.append(hpoTerm)

len(hpoNameEssay)   # 23,626 gene names as string elements of a single list
type(hpoNameEssay)  # <class 'str'>

# Randomly shuffle the hpoNameEssay
random.shuffle(hpoNameEssay)
random.shuffle(hpoNameEssay)

# Convert the randomized hpoNameEssay listVar to a stringVar
strHPOnameEssay = ', '.join(hpoNameEssay)

len(strHPOnameEssay)        # 703370 characters

# String maintenance (remove underscore, and replace with space)
strHPOnameEssay = strHPOnameEssay.replace('_', ' ')
# May need to clean up common string terms later as well




############ Generate WordClouds on HPO TERM FREQUENCIES:



# ## PNG file block: 

# # Pass to a dedicated text variable:
# text = strHPOnameEssay.lower()
# # Changing optional word cloud arguments
# # Now, change some optional arguments of the word cloud like max_font_size, max_word, and background_color.
# # lower max_font_size, change the maximum number of word and lighten the background:
# # wordcloud = WordCloud(max_font_size=50, max_words=100, background_color="white").generate(text)

# # Create stopword list:
# stopwords = set(STOPWORDS)
# stopwords.update(['Abnormal', 'involving', 'system', 'Morphological', 'count', 'morphology', 'forming', 'Abnormality',
                  # 'concentration', 'central', 'tissue', 'physiology', 'activity', 'physiology', 'abnormal', 'abnormality',
                  # 'of', 'the'])
# # Generate a word cloud image
# wordcloud = WordCloud(max_font_size=50, max_words=40, stopwords=stopwords, background_color="white", colormap='plasma').generate(text)
# # Display the generated image:
# plt.imshow(wordcloud, interpolation='bilinear')
# plt.axis("off")
# plt.savefig('./out/wordCloud_hpoTermFreq_HPOenrichment.png')
# # plt.show()



## SVG file block: 

# Pass to a dedicated text variable:
text = strHPOnameEssay.lower()
# Changing optional word cloud arguments
# Now, change some optional arguments of the word cloud like max_font_size, max_word, and background_color.
# lower max_font_size, change the maximum number of word and lighten the background:
# wordcloud = WordCloud(max_font_size=50, max_words=100, background_color="white").generate(text)

# Create stopword list:
stopwords = set(STOPWORDS)
stopwords.update(['Abnormal', 'involving', 'system', 'Morphological', 'count', 'morphology', 'forming', 'Abnormality',
                  'concentration', 'central', 'tissue', 'physiology', 'activity', 'physiology', 'abnormal', 'abnormality',
                  'of', 'the'])
# Generate a word cloud image

# wordcloud = WordCloud(max_font_size=50, max_words=40, stopwords=stopwords, background_color="white", colormap='plasma').generate(text)

wordcloud = WordCloud(width=1600, height=800, max_words=40, stopwords=stopwords, colormap='plasma', background_color='#F7F2FE').generate(text)
plt.figure(figsize=(20,10))
# Display the generated image:
plt.imshow(wordcloud, interpolation='bilinear')
plt.axis("off")
plt.savefig('./out/wordCloud_hpoTermFreq_HPOenrichment.svg')
# plt.show()