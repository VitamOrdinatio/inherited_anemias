## Handle collated GTR conditions list to extract underlying loci

## Handling VENV pocket:

#### Web-Scraping Environment
# python -m venv vScrapeClinVar
# source vScrapeClinVar/bin/activate
# #Make sure pip is up to date
# python -m pip install --upgrade pip
# #Install many useful packages for the vScrape environment
# python -m pip install --upgrade pip
# python -m pip install bs4
# python -m pip install pandas
# python -m pip install requests
# python -m pip install html5lib
# python -m pip install matplotlib
# python -m pip install numpy

#Import libraries:
import requests
import pandas as pd
from bs4 import BeautifulSoup
import numpy as np
import re

# target URL
# https://www.ncbi.nlm.nih.gov//gtr/conditions/C3469521/


# target NCBI GTR Conditions text to scrape:
# Associated genes section:
# HTML block:
    # <ul class=associate_genes gtr-reset-list>
        # <li>
            # <div>
                # <h3>
                    # <a href="/gtr/genes/"
                    #     ref="ncbi_uid="
                    #    data-ga-label="gtr-genes-page">FANCA
                    # </a>
                # </h3>
            # </div>
        # </li>
    # </ul>
   

#################################
# Initialize two list variables containing 178 genetic conditions:
### 171 from GTR query on "anemia" with OMIM / GeneReview filter facets applied
### 5   from GTR query on "willebrand" disease with OMIM / GeneReview filter facets applied
### 2   from GTR query on "hereditary factor" with OMIM / GeneReview filter facets applied (captures HEMOPHILIA)

# condLIST stores 178 names of NCBI GTR conditions that present with Anemia OR Willebrand disease OR Hemophilia
# linkLIST stores 178 GTR-condition hyperlinks of NCBI GTR conditions that present with Anemia OR Willebrand disease OR Hemophilia

condLIST = ['Fanconi anemia complementation group A', 'Hb SS disease', 'Fanconi anemia complementation group C', 'Fanconi anemia complementation group B', 'Fanconi anemia complementation group G', 'Fanconi anemia complementation group D2', 'Fanconi anemia complementation group I', 'Fanconi anemia complementation group E', 'Fanconi anemia complementation group F', 'Fanconi anemia complementation group L', 'beta Thalassemia', 'Diamond-Blackfan anemia 1', 'Fanconi anemia complementation group D1', 'Fanconi anemia complementation group J', 'Fanconi anemia complementation group N', 'Fanconi anemia complementation group O', 'Methylcobalamin deficiency type cblG', 'Megaloblastic anemia, thiamine-responsive, with diabetes mellitus and sensorineural deafness', 'Fanconi anemia complementation group P', 'Methylcobalamin deficiency type cblE', 'Diamond-Blackfan anemia 6', 'Diamond-Blackfan anemia 5', 'Hemolytic uremic syndrome, atypical, susceptibility to, 1', 'Majeed syndrome', 'Diamond-Blackfan anemia 10', 'Diamond-Blackfan anemia 7', 'Diamond-Blackfan anemia 9', 'Fanconi anemia, complementation group S', 'Diamond-Blackfan anemia 8', 'Diamond-Blackfan anemia 3', 'Thrombocytopenia, X-linked, with or without dyserythropoietic anemia', 'Fanconi anemia complementation group Q', 'Childhood onset GLUT1 deficiency syndrome 2', 'Diamond-Blackfan anemia 11', 'Insulin-dependent diabetes mellitus secretory diarrhea syndrome', 'Diamond-Blackfan anemia 4', 'Diamond-Blackfan anemia 12', 'Pearson syndrome', 'Congenital dyserythropoietic anemia, type I', 'Atypical hemolytic-uremic syndrome with MCP/CD46 anomaly', 'Fanconi anemia complementation group U', 'Familial juvenile hyperuricemic nephropathy type 2', 'Atypical hemolytic-uremic syndrome with C3 anomaly', 'Renal tubular acidosis, distal, 4, with hemolytic anemia', 'Atypical hemolytic-uremic syndrome with I factor anomaly', 'Atypical hemolytic-uremic syndrome with B factor anomaly', 'Diamond-Blackfan anemia 13', 'Atypical hemolytic-uremic syndrome with thrombomodulin anomaly', 'Congenital dyserythropoietic anemia type type 1B', 'Gaucher disease type I', 'alpha Thalassemia', 'Diamond-Blackfan anemia 14 with mandibulofacial dysostosis', 'Diamond-Blackfan anemia 15 with mandibulofacial dysostosis', 'Cobalamin C disease', 'Peutz-Jeghers syndrome', 'Fanconi anemia complementation group T', 'Metaphyseal chondrodysplasia, McKusick type', 'Multiple endocrine neoplasia, type 1', 'Gaucher disease type III', 'Wiskott-Aldrich syndrome', 'Abetalipoproteinaemia', 'Gaucher disease type II', 'Methylmalonic aciduria and homocystinuria type cblD', 'Nephronophthisis 1', 'Fanconi anemia complementation group R', 'Deficiency of ferroxidase', 'Diamond-Blackfan anemia 19', 'Diamond-Blackfan anemia 20', 'Neonatal intrahepatic cholestasis due to citrin deficiency', 'Diamond-Blackfan anemia 18', 'Methylmalonic aciduria and homocystinuria type cblF', 'Fanconi anemia complementation group V', 'Hyper-IgM syndrome type 1', 'Kearns-Sayre syndrome', 'Celiac disease, susceptibility to, 1', 'Fanconi anemia, complementation group W', 'Citrullinemia type II', 'Lysinuric protein intolerance', 'Jervell and Lange-Nielsen syndrome 1', 'Brain small vessel disease 1 with or without ocular anomalies', 'Cutaneous porphyria', 'Gaucher disease perinatal lethal', 'Autosomal dominant osteopetrosis 2', 'Nephronophthisis 3', 'Tangier disease', 'Nephronophthisis 4', 'Nephronophthisis 11', 'Infantile nephronophthisis', 'Thrombocytopenia 1', 'Renal tubular acidosis with progressive nerve deafness', 'Gaucher disease-ophthalmoplegia-cardiovascular calcification syndrome', 'X-linked severe congenital neutropenia', 'Jervell and Lange-Nielsen syndrome 2', 'Metaphyseal dysplasia without hypotrichosis', 'Anauxetic dysplasia 1', 'Retinal vasculopathy with cerebral leukoencephalopathy and systemic manifestations', 'Congenital defect of folate absorption', 'Generalized juvenile polyposis/juvenile polyposis coli', 'Methylmalonic acidemia with homocystinuria, type cblX', 'Autosomal dominant distal renal tubular acidosis', 'Autosomal dominant familial hematuria-retinal arteriolar tortuosity-contractures syndrome', 'Autosomal recessive osteopetrosis 4', 'Beta-thalassemia HBB/LCRB', 'Sitosterolemia 1', 'Prolidase deficiency', 'McLeod neuroacanthocytosis syndrome', 'Nephronophthisis 9', 'Sitosterolemia', 'Beta-thalassemia-X-linked thrombocytopenia syndrome', 'Diamond-Blackfan anemia 2', 'Nephronophthisis 15', 'Hereditary spherocytosis type 5', 'Nephronophthisis 7', 'DDX41-related hematologic malignancy predisposition syndrome', 'Jervell and Lange-Nielsen syndrome', 'Nephronophthisis 16', 'Nephronophthisis 13', 'Nephronophthisis 19', 'Wilms tumor 2', 'Nephronophthisis-like nephropathy 1', 'Nephronophthisis 18', 'MIRAGE syndrome', 'Wilms tumor 4', 'Pulmonary fibrosis and/or bone marrow failure syndrome, telomere-related, 7', 'Bile acid conjugation defect 1', 'Tyrosinemia type I', 'Propionic acidemia', 'Fabry disease', 'Shwachman-Diamond syndrome 1', 'Ehlers-Danlos syndrome, type 4', 'Infantile hypophosphatasia', 'Dyskeratosis congenita, X-linked', 'Niemann-Pick disease, type B', 'Chediak-Higashi syndrome', 'Telangiectasia, hereditary hemorrhagic, type 1', 'X-linked agammaglobulinemia', 'Wilson disease', 'DiGeorge syndrome', 'Mitochondrial complex IV deficiency, nuclear type 1', 'Recessive dystrophic epidermolysis bullosa', 'Radial aplasia-thrombocytopenia syndrome', 'Methylmalonic aciduria, cblA type', 'Telangiectasia, hereditary hemorrhagic, type 2', 'Neutropenia, severe congenital, 1, autosomal dominant', 'Methylmalonic aciduria, cblB type', 'Juvenile polyposis/hereditary hemorrhagic telangiectasia syndrome', 'Pulmonary fibrosis and/or bone marrow failure, Telomere-related, 1', 'Schimke immuno-osseous dysplasia', 'Dyskeratosis congenita, autosomal dominant 1', 'Diaphyseal dysplasia', 'Cerebroretinal microangiopathy with calcifications and cysts 1', 'Dyskeratosis congenita, autosomal dominant 3', 'Coenzyme Q10 deficiency, primary, 1', 'Hemochromatosis type 3', 'Autosomal recessive severe congenital neutropenia due to G6PC3 deficiency', 'Epidermolysis bullosa simplex 5B, with muscular dystrophy', 'Hemochromatosis type 2B', 'Pyridoxal phosphate-responsive seizures', 'Combined oxidative phosphorylation defect type 14', 'ALG8 congenital disorder of glycosylation', 'Revesz syndrome', 'Trimethylaminuria', 'Tubulointerstitial kidney disease, autosomal dominant, 2', 'COG1 congenital disorder of glycosylation', 'Thrombocytopenia 5', 'Aicardi-Goutieres syndrome 7', 'Cholesteryl ester storage disease', 'Autosomal dominant Kenny-Caffey syndrome', 'Ataxia-pancytopenia syndrome', 'Vasculitis due to ADA2 deficiency', 'Joubert syndrome with oculorenal defect', 'von Willebrand disease type 1', 'von Willebrand disease type 2', 'von Willebrand disease type 3', 'Glycogen storage disease due to glucose-6-phosphatase deficiency type IA', 'Glucose-6-phosphate transport defect', 'Hereditary factor IX deficiency disease', 'Hereditary factor VIII deficiency disease']


linkLIST = ['https://www.ncbi.nlm.nih.gov//gtr/conditions/C3469521/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0002895/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3468041/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1845292/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3469527/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3160738/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1836861/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3160739/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3469526/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3469528/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0005283/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2676137/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1838457/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1836860/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1835817/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3150653/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1855128/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0342287/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3469542/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1856057/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2931850/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2675859/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2749604/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1864997/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2750080/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2675512/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2750081/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4554406/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2675511/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1857719/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3550789/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3808988/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1842534/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3554042/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0342288/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2675860/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3809888/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0342784/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0271933/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2752040/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4310651/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2751310/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2752037/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C5436235/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2752039/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2752038/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4014641/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2752036/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3810185/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1961835/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0002312/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4225422/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4225411/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1848561/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0031269/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4084840/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0220748/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0025267/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0268251/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0043194/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0000744/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0268250/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1848552/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1855681/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4284093/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0878682/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C5193021/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C5193022/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1853942/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C5193020/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1848578/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4310652/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0398689/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0022541/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1859310/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4521564/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1863844/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0268647/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4551509/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4551998/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0162530/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1842704/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3179239/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1858392/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0039292/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1847013/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3150796/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1865872/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1839163/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0403554/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1856476/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1845987/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2676723/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1834821/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4551965/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1860518/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0342705/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1868081/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0796208/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/CN280572/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2673195/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1969106/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/CN322236/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2749759/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0268532/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0398568/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3151188/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0342907/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1839161/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1853666/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3541853/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2675192/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1969092/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4225174/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0022387/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3809320/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3280612/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4015542/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3887743/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3150419/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3890591/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4284088/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1832426/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C5830485/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C5543203/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0268490/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0268579/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0002986/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4692625/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0268338/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0268412/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1148551/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0268243/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0007965/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4551861/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0221026/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0019202/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0012236/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C5435656/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0079474/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0175703/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1855109/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1838163/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1859966/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1855102/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1832942/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3553617/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0877024/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4551974/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0011989/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4552029/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3151445/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3551954/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1858664/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2751630/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2931072/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1865616/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1864723/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4755312/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2931002/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1327916/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0342739/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1868139/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2931011/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4015537/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3888244/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0008384/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C4316787/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1327919/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C3887654/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1855675/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1264039/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1264040/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C1264041/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C2919796/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0268146/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0008533/', 'https://www.ncbi.nlm.nih.gov//gtr/conditions/C0019069/']



#############################
# target URL
# https://www.ncbi.nlm.nih.gov//gtr/conditions/C3469521/


# target NCBI GTR Conditions text to scrape:
# Associated genes section:
# HTML block:
    # <ul class=associate_genes gtr-reset-list>
        # <li>
            # <div>
                # <h3>
                    # <a href="/gtr/genes/"
                    #     ref="ncbi_uid="
                    #    data-ga-label="gtr-genes-page">FANCA
                    # </a>
                # </h3>
            # </div>
        # </li>
    # </ul>

################ MONOGENIC EXAMPLE
# url = 'https://www.ncbi.nlm.nih.gov/gtr/conditions/C3469521/'
# data = requests.get(url).text
# soup = BeautifulSoup(data, 'html.parser')

# elem = soup.find_all('a', attrs={'href': re.compile("^/gtr/genes/")}) 
# elem[0].text
## 'FANCA'

################ POLYGENIC EXAMPLE

#### beta Thalassemia is polygenic: HBB + HBB-LCR
# url for this polygenic condition:
# btURL = 'https://www.ncbi.nlm.nih.gov/gtr/conditions/C0005283/'
# btDATA = requests.get(btURL).text
# btSOUP = BeautifulSoup(btDATA, 'html.parser')

# btELEM = btSOUP.find_all('a', attrs={'href': re.compile("^/gtr/genes/")}) 

# btELEM[0].text
## 'HBB'
# btELEM[1].text
## 'HBB-LCR'

## Reminder: 
##          cycle through condLIST to get genetic condition names
##          cycle through linkLIST to scrape gene names
##          ##  cycle through bs4 search result

# len(condLIST)
## 171
# len(linkLIST)
## 171

## dummy link var:
# link = linkLIST[0]   # FANCA (monogenic: FANCA)

# link = linkLIST[10]  # beta Thalassemia (polygenic: HBB + HBB-LCR)


# Initialize accumulation variables
accumCONDS = []
accumGENES = []
i = 0

for link in linkLIST:
    try:
        #Request HTML page
        data = requests.get(link).text
        # Create soup object
        soup = BeautifulSoup(data, 'html.parser')
        elem = soup.find_all('a', attrs={'href': re.compile("^/gtr/genes/")}) 
        numGenes = len(elem)
        for gene in range(numGenes):
            accumCONDS.append(condLIST[i])
            accumGENES.append(elem[gene].text)
            print(elem[gene].text)
        # Increment counter
        print('link ' + str(i+1) + ' was scraped successfully')
        i = (i + 1)
    except:
        print('link ' + str(i+1) + ' had an error')
        i = (i + 1)
        pass

dfGTR_Genes = pd.DataFrame()
dfGTR_Genes["Genetic_Disorder"] = accumCONDS     # 178 unique genetic conditions
dfGTR_Genes["Gene_Name"] = accumGENES        
dfGTR_Genes.to_csv('NCBI_GTR_178_hematology_genetic_cond_GENE_etiology.csv')



