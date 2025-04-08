## Scrape NCBI ClinVar for allele metadata regarding a list of 199 anemia-enriched loci

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

#Import libraries:
import requests
import pandas as pd
from bs4 import BeautifulSoup
import numpy as np

#Define a function called search that looks for substrings within panda DataFrames:
def search(df: pd.DataFrame, substring: str, case: bool = False) -> pd.DataFrame:
    mask = np.column_stack([df[col].astype(str).str.contains(substring.lower(), case=case, na=False) for col in df])
    return df.loc[mask.any(axis=1)]
    
# test search function
##df = pd.DataFrame({'col1':['hello', 'world', 'Sun'], 'col2': ['today', 'sunny', 'foo'], 'col3': ['WORLD', 'NEWS', 'bar']})
##print(df)
##search(df, 'sun')

#Dummy list (test run)
#list_of_GeneNames=['AARS1','ABAT']

# Full list of 164 genes involved in anemia manifestations across 178 NCBI GTR hematology conditions
## 171 unique anemia conditions involving 159 genes
## 5   unique Willebrand conditions involving 3 genes (VWF, G6PC1, and SLC37A4)
## 2   unique hemophilia conditions involving 2 genes (F8, and F9) which are hereditary factors
list_of_GeneNames=['ABCA1', 'ABCG8', 'ACVRL1', 'ADA2', 'ALG8', 'ALPL', 'ANKS6', 'ATP6V1B1', 'ATP7B', 'BAAT', 'BRCA1', 'BRCA2', 'BRIP1', 'BTK', 'C3', 'CD40LG', 'CD46', 'CDAN1', 'CDIN1', 'CEP164', 'CEP83', 'CFB', 'CFH', 'CFHR1', 'CFHR3', 'CFI', 'CLCN7', 'COG1', 'COL3A1', 'COL4A1', 'COL7A1', 'COQ2', 'CP', 'CTC1', 'DCDC2', 'DDX41', 'DKC1', 'DNAJC21', 'ELANE', 'ENG', 'EPB42', 'ERCC4', 'ETV6', 'F8', 'F9', 'FAH', 'FAM111A', 'FANCA', 'FANCB', 'FANCC', 'FANCD2', 'FANCE', 'FANCF', 'FANCG', 'FANCI', 'FANCL', 'FARS2', 'FMO3', 'FOXP3', 'G6PC1', 'G6PC3', 'GATA1', 'GBA1', 'GLA', 'GLIS2', 'H19-ICR', 'HAMP', 'HBA1', 'HBA2', 'HBB', 'HBB-LCR', 'HCFC1', 'HLA-DQA1', 'HLA-DQB1', 'IFIH1', 'INVS', 'KCNE1', 'KCNQ1', 'LIPA', 'LMBRD1', 'LPIN2', 'LYST', 'MAD2L2', 'MEN1', 'MMAA', 'MMAB', 'MMACHC', 'MMADHC', 'MMP1', 'MTR', 'MTRR', 'MT-TN', 'MTTP', 'MT-TS1', 'MUC1', 'NAF1', 'NEK8', 'NPHP1', 'NPHP3', 'NPHP4', 'PALB2', 'PCCA', 'PCCB', 'PEPD', 'PLEC', 'PNPO', 'PRDX1', 'RAD51', 'RAD51C', 'RBM8A', 'REN', 'RFWD3', 'RMRP', 'RPL11', 'RPL15', 'RPL18', 'RPL26', 'RPL35', 'RPL35A', 'RPL5', 'RPS10', 'RPS15A', 'RPS17', 'RPS19', 'RPS24', 'RPS26', 'RPS28', 'RPS29', 'RPS7', 'SAMD9', 'SAMD9L', 'SBDS', 'SLC19A2', 'SLC25A13', 'SLC2A1', 'SLC37A4', 'SLC46A1', 'SLC4A1', 'SLC7A7', 'SLX4', 'SMAD4', 'SMARCAL1', 'SMPD1', 'SRP54', 'STK11', 'SURF1', 'TBX1', 'TERC', 'TERT', 'TFR2', 'TGFB1', 'THBD', 'TINF2', 'TMEM67', 'TREX1', 'TSR2', 'UBE2T', 'UROS', 'VWF', 'WAS', 'WDR19', 'XK', 'XPNPEP3', 'XRCC2']

# Cycle through ClinVar gene-search query result pages:
for gene in list_of_GeneNames:
    try:
        # Construct URL
        prefixURL = "https://www.ncbi.nlm.nih.gov/clinvar/?term="
        suffixURL = "%5bgene%5d"
        url = prefixURL + gene + suffixURL
        # Handshake with requests, extract page as text
        data = requests.get(url).text
        # Creating BeautifulSoup object
        soup = BeautifulSoup(data, 'html.parser')
        ### Verifying tables and their classes
        ### At NCBI ClinVar, I see 6 distinct tables
        # print('Classes of each table:')
        # for table in soup.find_all('table'):
            # print(table.get('class'))   
        ### Verifying descriptions lists (dl), and their classes
        ### At NCBI ClinVar, I see 3 distinct description lists:
        # print('Classes of each description list:')
        # for lists in soup.find_all('dl'):
            # print(lists.get('class'))
        ### Verifying sub-headings <h2> and their classes
        ### At NCBI Clinar, I see a single sub-heading <h2> tag:
        # print('Classes of each h2 sub-heading tag:')
        # for subheads in soup.find_all('h2'):
            # print(subheads.get('class'))
        #
        ### Verify existence of an HGVS table
        # print('Classes of each table:')
        # for table in soup.find_all('table'):
            # print(table.get('class'))
        # Filter soup.  Target is <li class="fil_val">
        facetInfo = soup.find_all("li", {'class':'fil_val'})
        # Pull out text values
        GeneName = gene
        # Germline alleles
        Germline = facetInfo[0].text
        Germline = Germline.strip()
        Germline = Germline.replace("Germline (","")
        Germline = Germline.replace(")","")
        # Somatic alleles
        Somatic = facetInfo[1].text
        Somatic = Somatic.strip()
        Somatic = Somatic.replace("Somatic (","")
        Somatic = Somatic.replace(")","")
        # Conflicting Classications
        Germline_ConflictingClassifications = facetInfo[2].text
        Germline_ConflictingClassifications = Germline_ConflictingClassifications.strip()
        Germline_ConflictingClassifications = Germline_ConflictingClassifications.replace("Conflicting classifications (","")
        Germline_ConflictingClassifications = Germline_ConflictingClassifications.replace(")","")        
        # Benign classifications
        Benign = facetInfo[3].text
        Benign = Benign.strip()
        Benign = Benign.replace("Benign (","")
        Benign = Benign.replace(")","")
        # Likely benign classifications
        LikelyBenign = facetInfo[4].text
        LikelyBenign = LikelyBenign.strip()
        LikelyBenign = LikelyBenign.replace("Likely benign (","")
        LikelyBenign = LikelyBenign.replace(")","")
        # Uncertain significance classifications
        UncertainSignificance = facetInfo[5].text
        UncertainSignificance = UncertainSignificance.strip()
        UncertainSignificance = UncertainSignificance.replace("Uncertain significance (","")
        UncertainSignificance = UncertainSignificance.replace(")","")
        # Likely pathogenic classifications
        LikelyPathogenic = facetInfo[6].text
        LikelyPathogenic = LikelyPathogenic.strip()
        LikelyPathogenic = LikelyPathogenic.replace("Likely pathogenic (","")
        LikelyPathogenic = LikelyPathogenic.replace(")","")
        # Pathogenic classifications
        Pathogenic = facetInfo[7].text
        Pathogenic = Pathogenic.strip()
        Pathogenic = Pathogenic.replace("Pathogenic (","")
        Pathogenic = Pathogenic.replace(")","")
        # P/LP vs. LB/B ratio
        P_LP_vs_LB_B = facetInfo[8].text
        P_LP_vs_LB_B = P_LP_vs_LB_B.strip()
        P_LP_vs_LB_B = P_LP_vs_LB_B.replace("P/LP vs LB/B (","")
        P_LP_vs_LB_B = P_LP_vs_LB_B.replace(")","")
        # P/LP vs. VUS ratio
        P_LP_vs_VUS = facetInfo[9].text
        P_LP_vs_VUS = P_LP_vs_VUS.strip()
        P_LP_vs_VUS = P_LP_vs_VUS.replace("P/LP vs VUS (","")
        P_LP_vs_VUS = P_LP_vs_VUS.replace(")","")
        # VUS vs LB/B ratio
        VUS_vs_LB_B = facetInfo[10].text
        VUS_vs_LB_B = VUS_vs_LB_B.strip()
        VUS_vs_LB_B = VUS_vs_LB_B.replace("VUS vs LB/B (","")
        VUS_vs_LB_B = VUS_vs_LB_B.replace(")","")
        # Frameshift alleles
        Frameshift = facetInfo[11].text
        Frameshift = Frameshift.strip()
        Frameshift = Frameshift.replace("Frameshift (","")
        Frameshift = Frameshift.replace(")","")
        # Missense alleles
        Missense = facetInfo[12].text
        Missense = Missense.strip()
        Missense = Missense.replace("Missense (","")
        Missense = Missense.replace(")","")
        # Nonsense alleles
        Nonsense = facetInfo[13].text
        Nonsense = Nonsense.strip()
        Nonsense = Nonsense.replace("Nonsense (","")
        Nonsense = Nonsense.replace(")","")
        # SpliceSite alleles
        SpliceSite = facetInfo[14].text
        SpliceSite = SpliceSite.strip()
        SpliceSite = SpliceSite.replace("Splice site (","")
        SpliceSite = SpliceSite.replace(")","")
        # ncRNA alleles
        ncRNA = facetInfo[15].text
        ncRNA = ncRNA.strip()
        ncRNA = ncRNA.replace("ncRNA (","")
        ncRNA = ncRNA.replace(")","")
        # NearGene alleles
        NearGene = facetInfo[16].text
        NearGene = NearGene.strip()
        NearGene = NearGene.replace("Near gene (","")
        NearGene = NearGene.replace(")","")
        # UTR alleles
        UTR = facetInfo[17].text
        UTR = UTR.strip()
        UTR = UTR.replace("UTR (","")
        UTR = UTR.replace(")","")
        # Deletion alleles
        Deletion = facetInfo[18].text
        Deletion = Deletion.strip()
        Deletion = Deletion.replace("Deletion (","")
        Deletion = Deletion.replace(")","")
        # Duplication alleles
        Duplication = facetInfo[19].text
        Duplication = Duplication.strip()
        Duplication = Duplication.replace("Duplication (","")
        Duplication = Duplication.replace(")","")
        # Indel alleles
        Indel = facetInfo[20].text
        Indel = Indel.strip()
        Indel = Indel.replace("Indel (","")
        Indel = Indel.replace(")","")
        # Insertion alleles
        Insertion = facetInfo[21].text
        Insertion = Insertion.strip()
        Insertion = Insertion.replace("Insertion (","")
        Insertion = Insertion.replace(")","")
        # SingleNucleotide alleles
        SingleNucleotide = facetInfo[22].text
        SingleNucleotide = SingleNucleotide.strip()
        SingleNucleotide = SingleNucleotide.replace("Single nucleotide (","")
        SingleNucleotide = SingleNucleotide.replace(")","")
        # ShortVariant alleles
        ShortVariant = facetInfo[23].text
        ShortVariant = ShortVariant.strip()
        ShortVariant = ShortVariant.replace("Short variant (< 50 bps) (","")
        ShortVariant = ShortVariant.replace(")","")
        # StructuralVariant alleles
        StructuralVariant = facetInfo[24].text
        StructuralVariant = StructuralVariant.strip()
        StructuralVariant = StructuralVariant.replace("Structural variant (>= 50 bps) (","")
        StructuralVariant = StructuralVariant.replace(")","")
        # Less1KB_singleGene alleles
        Less1KB_singleGene = facetInfo[25].text
        Less1KB_singleGene = Less1KB_singleGene.strip()
        Less1KB_singleGene = Less1KB_singleGene.replace("< 1kb, single gene (","")
        Less1KB_singleGene = Less1KB_singleGene.replace(")","")
        # More1KB_singleGene alleles
        More1KB_singleGene = facetInfo[26].text
        More1KB_singleGene = More1KB_singleGene.strip()
        More1KB_singleGene = More1KB_singleGene.replace("> 1kb, single gene (","")
        More1KB_singleGene = More1KB_singleGene.replace(")","")
        # More1KB_multipleGenes alleles
        More1KB_multipleGenes = facetInfo[27].text
        More1KB_multipleGenes = More1KB_multipleGenes.strip()
        More1KB_multipleGenes = More1KB_multipleGenes.replace("> 1kb, multiple genes (","")
        More1KB_multipleGenes = More1KB_multipleGenes.replace(")","")
        # PracticeGuideline alleles
        PracticeGuideline = facetInfo[28].text
        PracticeGuideline = PracticeGuideline.strip()
        PracticeGuideline = PracticeGuideline.replace("Practice guideline (","")
        PracticeGuideline = PracticeGuideline.replace(")","")
        # ExpertPanel alleles
        ExpertPanel = facetInfo[29].text
        ExpertPanel = ExpertPanel.strip()
        ExpertPanel = ExpertPanel.replace("Expert panel (","")
        ExpertPanel = ExpertPanel.replace(")","")
        # MultipleSubmitters alleles
        MultipleSubmitters = facetInfo[30].text
        MultipleSubmitters = MultipleSubmitters.strip()
        MultipleSubmitters = MultipleSubmitters.replace("Multiple submitters (","")
        MultipleSubmitters = MultipleSubmitters.replace(")","")
        # SingleSubmitter alleles
        SingleSubmitter = facetInfo[31].text
        SingleSubmitter = SingleSubmitter.strip()
        SingleSubmitter = SingleSubmitter.replace("Single submitter (","")
        SingleSubmitter = SingleSubmitter.replace(")","")
        # AtLeastOneStar alleles
        AtLeastOneStar = facetInfo[32].text
        AtLeastOneStar = AtLeastOneStar.strip()
        AtLeastOneStar = AtLeastOneStar.replace("At least one star (","")
        AtLeastOneStar = AtLeastOneStar.replace(")","")
        # RevStatus_ConflictingClassifications alleles
        RevStatus_ConflictingClassifications = facetInfo[33].text
        RevStatus_ConflictingClassifications = RevStatus_ConflictingClassifications.strip()
        RevStatus_ConflictingClassifications = RevStatus_ConflictingClassifications.replace("Conflicting classifications (","")
        RevStatus_ConflictingClassifications = RevStatus_ConflictingClassifications.replace(")","")

        # #Useful variables
        # Germline
        # Somatic
        # Germline_ConflictingClassifications
        # Benign
        # LikelyBenign
        # UncertainSignificance
        # LikelyPathogenic
        # Pathogenic
        # P_LP_vs_LB_B
        # P_LP_vs_VUS
        # VUS_vs_LB_B
        # Frameshift
        # Missense
        # Nonsense
        # SpliceSite
        # ncRNA
        # NearGene
        # UTR
        # Deletion
        # Duplication
        # Indel
        # Insertion
        # SingleNucleotide
        # ShortVariant
        # StructuralVariant
        # Less1KB_singleGene
        # More1KB_singleGene
        # More1KB_multipleGenes
        # PracticeGuideline
        # ExpertPanel
        # MultipleSubmitters
        # SingleSubmitter
        # AtLeastOneStar
        # RevStatus_ConflictingClassifications
        #
        #Create a list of dictionaries:
        list_of_dicts = [
             {"Gene name":gene,"Germline":Germline,"Somatic":Somatic,"Germline: Conflicting classiciations":Germline_ConflictingClassifications,"Benign":Benign,"Likely benign":LikelyBenign,"Uncertain significance":UncertainSignificance,"Likely pathogenic":LikelyPathogenic,"Pathogenic":Pathogenic,"P/LP vs. LB/B":P_LP_vs_LB_B,"P/LP vs. VUS":P_LP_vs_VUS,"VUS vs. LB/B":VUS_vs_LB_B,"Frameshift":Frameshift,"Missense":Missense,"Nonsense":Nonsense,"Splice site":SpliceSite,"ncRNA":ncRNA,"Near gene":NearGene,"UTR":UTR,"Deletion":Deletion,"Duplication":Duplication,"Indel":Indel,"Insertion":Insertion,"Single nucleotide":SingleNucleotide,"Short variant (< 50 bps)":ShortVariant,"Structural variant (>= 50bps)":StructuralVariant,"Less than 1kb, single gene":Less1KB_singleGene,"More than 1 kb, single gene":More1KB_singleGene,"More than 1 kb, multiple genes":More1KB_multipleGenes,"Practice guideline":PracticeGuideline,"Expert panel":ExpertPanel,"Multiple submitters":MultipleSubmitters,"Single submitter":SingleSubmitter,"At least one star":AtLeastOneStar,"Review status: Conflicting classifications":RevStatus_ConflictingClassifications}
        ]
        #Create a new pandas data-frame that stores a row entry of 10 useful variables:
        newdf = pd.DataFrame(list_of_dicts)
        #Write single ClinVar as CSV file
        newdf.to_csv('%s.csv' % gene)
        ##### This CODE BLOCK TRHOWS ERRORS:
        #Concat to growing data-frame
        ##dfBIG.concat([dfBIG, newdf], ignore_index=TRUE)
        print(gene, 'WAS SUCCESSFULLY SCRAPED!')        
    except:
        print(gene, ' had an error')
        pass