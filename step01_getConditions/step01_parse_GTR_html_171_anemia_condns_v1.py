## Activate your Virtual Environment (vScrapeClinVar)
### navigate to .venv folder in wslshare
### run this command:
### source vScrapeClinvar/bin/activate
### enter python environment

## HTML syntax:
# td = columns (TD stands for table data)
# tr = rows (TR stands for table rows)

# dd = definition description
# dl = description list (contains many dt terms)
# dt = description terms (belong to a dl list)

# h2 = sub-heading tag

## Modified this geeksforgeeks.org tutorial: https://www.geeksforgeeks.org/beautifulsoup-scraping-link-from-html/

#Import libraries:
import requests
import pandas as pd
from bs4 import BeautifulSoup
import numpy as np
import re

####### Handling 1 of 3 NCBI GTR hematology queries:

## Initialize a list of URLs.  
## First 9 pages are results for a GTR query for "anemia" with OMIM / GeneReviews filters applied. 
list_of_URLs = ['https://www.ncbi.nlm.nih.gov/gtr/conditions/?term=anemia&filter=omim,gene_review', 'https://www.ncbi.nlm.nih.gov/gtr/conditions/?filter=omim%2Cgene_review&term=anemia&page=2', 'https://www.ncbi.nlm.nih.gov/gtr/conditions/?filter=omim%2Cgene_review&term=anemia&page=3', 'https://www.ncbi.nlm.nih.gov/gtr/conditions/?filter=omim%2Cgene_review&term=anemia&page=4', 'https://www.ncbi.nlm.nih.gov/gtr/conditions/?filter=omim%2Cgene_review&term=anemia&page=5', 'https://www.ncbi.nlm.nih.gov/gtr/conditions/?filter=omim%2Cgene_review&term=anemia&page=6', 'https://www.ncbi.nlm.nih.gov/gtr/conditions/?filter=omim%2Cgene_review&term=anemia&page=7', 'https://www.ncbi.nlm.nih.gov/gtr/conditions/?filter=omim%2Cgene_review&term=anemia&page=8', 'https://www.ncbi.nlm.nih.gov/gtr/conditions/?filter=omim%2Cgene_review&term=anemia&page=9']

## Willebrand Disease query with OMIM / GeneReviews filters applied returns 5 conditions:
# list_of_URLs = ['https://www.ncbi.nlm.nih.gov/gtr/conditions/?term=Willebrand&filter=omim,gene_review']

## "Hereditary factor" query with OMIM / GeneReviews filters applied returns 2 HEMOPHILIA conditions:
# list_of_URLs = ['https://www.ncbi.nlm.nih.gov/gtr/conditions/?term=Hereditary%20factor&filter=omim,gene_review']

###GTR scraping test run
# url = list_of_URLs[0]
# data = requests.get(url).text
# soup = BeautifulSoup(data, 'html.parser')

### Verifying tables and their classes
### At NCBI GTR, a query for 'anemia' with applied filters for both 'OMIM available' and 'GeneReviews available' yields 171 genetic conditions
####### This query result of 171 genetic conditions is arrayed via GTR across 9 html pages.  Each page can maximally carry 20 gen condn hits.
####### Each hit should have this HTML structure:
####### <div class="rprt" data-section="page-report">
#######     <div class="rprtnum">1.</div>
#######     <div class="rslt">
#######         <p class="title">
#######             <a href="/gtr/conditions/C3469521"
#######                ref="ncbi_uid=483333&link_uid=483333"
#######                data-ga-label="gtr-conditions-page">"Fanconi" <b>anemia</b>" complementation group A"
#######             </a>
#######         </p>
#######     </div>
####### </div>

### Verifying div tags, and their classes
# print('Classes of each div tag:')
# for div in soup.find_all('div'):
    # print(div.get('page_report'))

### facetInfo now has up to 20 hits
# facetInfo = soup.find_all("div", {'class':'rprt'})

### slice 1 of 20:
# hit = facetInfo[0]

# x = hit.find_all("a")
# len(x)
#### 7
# x[0]
#### <a href="/gtr/conditions/C3469521/" ref="ncbi_uid=483333&amp;link_uid=483333">Fanconi <b>anemia</b> complementation group A</a>
# x[1]
#### <a class="external-link" href="https://www.ncbi.nlm.nih.gov/books/NBK1116" target="blank" title="GeneReviews">GeneReviews</a>
# x[2]
#### <a href="/gtr/all/tests/?term=C3469521">Tests</a>
# x[3]
#### <a href="/gtr/all/labs/?term=C3469521">Labs</a>
# x[4]
#### <a href="https://www.ncbi.nlm.nih.gov/gtr/genes/?from_condition_id=483333">Genes</a>
# x[5]
#### <a href="/omim?LinkName=medgen_omim&amp;from_uid=483333">OMIM</a>
# x[6]
#### <a class="gene_reviews" href="/pubmed?LinkName=medgen_pubmed_genereviews&amp;from_uid=483333">GeneReviews</a>



### x[0].text has genetic condition name
# genCondn = x[0].text
### x[0].url has hyperlink

########################### 
# import necessary libraries 
# from bs4 import BeautifulSoup 
# import requests 
# import re 

# Initialize rolling lists for appending scrape data:
GeneticDisorderNAME = []        # Name of genetic disorder's condition/phenotype
GTRdisorderACCN = []            # NCBI GTR Accession number for a specific genetic disorder
GTRtestManifest = []            # NCBI GTR's hyperlink for available genetic tests for a specific disorder
GTRlabsManifest = []            # NCBI GTR's hyperlink for available genetic labs for a specific disorder
GTRgenesManifest = []           # NCBI GTR's hyperlink for known gene(s) underlying a specific disorder
redirectedOMIM_queryURL = []    # NCBI GTR's redirected hyperlink to OMIM for a specific genetic disorder
GeneReviewsURL = []             # NCBI GeneReview's page for a specific genetic disorder

##### Future loop starts here

# Cycle through ClinVar gene-search query result pages:
for url in list_of_URLs:
    try:
        # Request HTML page
        data = requests.get(url).text
        # Create soup object
        soup = BeautifulSoup(data, 'html.parser')

        #### Name of a precise genetic condition
        for link in soup.find_all('a', attrs={'href': re.compile("^/gtr/conditions")}): 
            # extract the partial url
            x = link.text
            # append to accumulating list variable
            GeneticDisorderNAME.append(x)
            # type(x)
            # print(x)   

        #### Accession Number on NCBI GTR for a precise genetic condition
        for link in soup.find_all('a', attrs={'href': re.compile("^/gtr/conditions")}): 
            # extract the partial url
            x = link.get('href')
            x = "https://www.ncbi.nlm.nih.gov/" + x
            # append to accumulating list variable
            GTRdisorderACCN.append(x)
            # type(x)
            # print(x)

        #### Available genetic TESTS for a specific disorder listed on NCBI GTR database
        for link in soup.find_all('a', attrs={'href': re.compile("^/gtr/all/tests/")}): 
            # extract the partial url
            x = link.get('href')
            x = "https://www.ncbi.nlm.nih.gov/" + x
            # append to accumulating list variable
            GTRtestManifest.append(x)
            # type(x)
            # print(x)   

        #### Available genetic testing LABS for a specific disorder listed on NCBI GTR database
        for link in soup.find_all('a', attrs={'href': re.compile("^/gtr/all/labs/")}): 
            # extract the partial url
            x = link.get('href')
            x = "https://www.ncbi.nlm.nih.gov/" + x
            # append to accumulating list variable
            GTRlabsManifest.append(x)
            # type(x)
            # print(x)   
            
        #### GENE(s) entry for a specific disorder listed on NCBI GTR database
        for link in soup.find_all('a', attrs={'href': re.compile("^https://www.ncbi.nlm.nih.gov/gtr/genes/")}): 
            # extract the partial url
            x = link.get('href')
            # append to accumulating list variable
            GTRgenesManifest.append(x)
            # type(x)
            # print(x)  

        #### OMIM redirect for a precise genetic condition
        for link in soup.find_all('a', attrs={'href': re.compile("^/omim")}): 
            # extract the partial url
            x = link.get('href')
            x = "https://www.ncbi.nlm.nih.gov/" + x    
            # append to accumulating list variable
            redirectedOMIM_queryURL.append(x)
            # type(x)
            # print(x) 

        #### GeneReviews entry a specific disorder listed on NCBI GTR database
        for link in soup.find_all('a', attrs={'class': "gene_reviews"}): 
            # extract the partial url
            x = link.get('href')
            x = "https://www.ncbi.nlm.nih.gov/" + x    
            # append to accumulating list variable
            GeneReviewsURL.append(x)
            # type(x)
            # print(x)


        # View first genetic disorder variables:
        # # Name of genetic disorder's condition/phenotype
        # GeneticDisorderNAME[0]
        # # NCBI GTR Accession number for a specific genetic disorder
        # GTRdisorderACCN[0]
        # # NCBI GTR's hyperlink for available genetic tests for a specific disorder
        # GTRtestManifest[0]
        # # NCBI GTR's hyperlink for available genetic labs for a specific disorder
        # GTRlabsManifest[0]
        # # NCBI GTR's hyperlink for known gene(s) underlying a specific disorder
        # GTRgenesManifest[0]
        # # NCBI GTR's redirected hyperlink to OMIM for a specific genetic disorder
        # redirectedOMIM_queryURL[0]
        # # NCBI GeneReview's page for a specific genetic disorder
        # GeneReviewsURL[0]

        # GeneticDisorderNAME[0]
        # GTRdisorderACCN[0]
        # GTRtestManifest[0]
        # GTRlabsManifest[0]
        # GTRgenesManifest[0]
        # redirectedOMIM_queryURL[0]
        # GeneReviewsURL[0]
        print(url, 'WAS SUCCESSFULLY SCRAPED!')        
    except:
        print(url, ' had an error')
        pass


dfGTR = pd.DataFrame()
dfGTR["Genetic_Disorder"] = GeneticDisorderNAME     # 171
dfGTR["GTR_Disorder_ACCN"] = GTRdisorderACCN        # 171
# dfGTR["GTR_Disorder_Tests"] = GTRtestManifest   # 167 len
# dfGTR["GTR_Disorder_Labs"] = GTRlabsManifest    # 167 len
# dfGTR["GTR_Disorder_Genes"] = GTRgenesManifest  # 163 len
dfGTR["OMIM_requery"] = redirectedOMIM_queryURL     # 171
dfGTR["Gene_Reviews"] = GeneReviewsURL              # 171

### Handling 1 of 3 output streams:
dfGTR.to_csv('NCBI_GTR_171genetic_cond_with_anemia.csv')
# dfGTR.to_csv('NCBI_GTR_5genetic_cond_with_willebrand.csv')
# dfGTR.to_csv('NCBI_GTR_2genetic_cond_with_hemophilia.csv')
