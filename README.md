highest

**STEP 1: Scrape NCBI GTR for genetic conditions**

We sought to comprehensively scrape all inherited hematological disorders with the highest curation statuses. To this end, we employed custom python scripts leveraging web-scraping modules like beautifulsoup4 and performed HTML scrape operations against the NCBI GTR database for each of the following queries 1) anemia, 2) willebrand and 3) hereditary factor.  In all instances, we scraped results with applied filters for OMIM and NCBI GeneReviews hits.

For anemia queries, we scraped 171 unique genetic conditions.
For willebrand queries, we scraped 5 unique genetic conditions.
For hereditary factor queries, we scraped 2 unique genetic conditions.

We then performed AWK operations to collate all 178 unique genetic conditions representing inherited hematological disorders of the highest curation statuses.

**STEP 2: Scrape NCBI GTR for genetic etiologies**

We next sought the underlying genetic etiologies (i.e., monogenic or polygenic) driving each of the inherited hematological disorders. After collation, we generated a list of 164 unique genes that are known to be involved in 170 unique hematological disorders.

**STEP 3: Assemble master locus list underlying inherited anemias**

Two scrape targets were not actually genes (i.e., they were regulatory control elements: H19-ICR and HBB-LCR) and were thus removed from our list of hematological genes, resulting in a list of 162 unique loci.

Two of these 162 scrape targets were mitochondrial genes encoded on the human mitochondrial chromosome (i.e., MT-TN and MT-TS1). A literature review revealed that such anemias are secondary manifestations that occassionally show up in mitochondrial disease patients. The nuclear-encoded POLG and POLG2 genes encode the mtDNA polymerase complex, and a loss-of-function in mtDNA polymerase activities underly POLG-related disorders, a subset of mitochondrial diseases. We manually collated all other mtDNA loci (n = 35 mitochondrial genes) as well as nuclear-encoded POLG1 and POLG2 (n = 2 nuclear genes) into the master locus list (n = 199 unique genes).

**STEP 4: Assemble master locus list underlying inherited anemias**

We scraped NCBI ClinVar for allele counts for each of 199 unique anemia-enriched genes. Alleles are classified in any one of six ClinVar categories: conflicting classifications (CC), benign (B), likely benign (LB), uncertain significance (US), likely pathogenic (LP), and pathogenic (P). Additional metadata for each scraped locus was extracted, including molecular consequences (i.e., frameshift, missense, nonsense, splice site, ncRNA, near gene, UTR, deletion, duplication, indel, insertion, single nucleotide) as well as large-scale chromosomal results (e.g., structural variant >= 50 bps).

**STEP 5: Get Human Phenotype Ontology (HPO) enrichment**

We employed g:profiler for gene set enrichment analysis (GSEA) using the HPO database. Our input comprised 199 unique loci and we extracted 619 HPO enrichment terms involving gene intersections. We generated a truth table of 199 gene rows by 619 HPO phenotypic columns. We next represented each gene in 619th-dimensional space, and performed topological data analysis (TDA) followed by CIRCOS visualization.

**STEP 6: Get MitoCarta gene status**

As the original scrape captured two mtDNA loci (i.e., MT-TN and MT-TS1), we wanted to examine how many other genes in our anemia-enriched gene list were also found in the MitoCarta database.  MitoCarta tracks 1,136 human genes that are nuclear-encoded but exhibit strong evidence for mitochondrial functions (i.e., localization, biochemistry, etc).  A total of 54 of the 199 genes are found in MitoCarta.  Of these 54, we manually added 37 (i.e., POLG1, POLG2, and 35 mtDNA genes). Thus, 17 of the original 162 scraped genes from NCBI GTR are found in the MitoCarta database.

**STEP 7: Perform multivariate statistical analyses on ClinVar allele categorical distributions**

Across all 199 anemia-enriched genetic loci, we scraped 192,296 unique alleles from each of six NCBI ClinVar categories (i.e., CC, B, LB, US, LP, and P). We removed alleles in the CC and US categories, to arrive at a dataset of 199 loci made up of 112,534 raw allele counts across four NCBI ClinVar categories (i.e., B, LB, LP, and P).  We next visualized our 112K allelic dataset by performing multivariate statistical analyses, including pairplots, t-SNE, linear regressions, correlation coefficients, principal component analysis (PCA), PCA explained variance ratios, unfiltered k-means clustering, and two-component (2PC) k-means clustering.

**STEP 8: Perform standard visualization analyses on ClinVar allele categorical distributions**

Across all 199 anemia-enriched genetic loci, we scraped 192,296 unique alleles from each of six NCBI ClinVar categories (i.e., CC, B, LB, US, LP, and P). We removed alleles in the CC and US categories, to arrive at a dataset of 199 loci made up of 112,534 raw allele counts across four NCBI ClinVar categories (i.e., B, LB, LP, and P).  We next visualized our 112K allelic dataset by generating the following output:

   1) Two CSV files that track either raw allele counts or normalized allelic categorical frequencies.
 
   2) A barplot of log2-transformed total allele counts across all 199 loci (PNG + SVG formats)
 
   3) A histogram plot of log2-transformed total allele counts across ten-percentile bins for all 199 loci (PNG + SVG formats)
 
   4) A pie chart of a slice of the top 20 loci after sorting the list of 199 genes by total allele counts (derived from 4 ClinVar categories)
 
   5) A horizontal barplot showing the relative number of alleles across 4 categories for each of the top 20 loci by total allele counts
 
   6) A pie chart of a slice of the top 20 loci after sorting the list of 199 genes by normalized allele frequencies (derived from 4 ClinVar categories)
 
   7) A horizontal barplot showing the relative number of alleles across 4 categories for each of the top 20 loci by normalized allele frequencies

**STEP 9: Prepare HPO gene ontology enrichment terms for topological data analysis (TDA)**

We manually loaded a list of 199 anemia-enriched loci into the online gene ontology (GO) / gene set enrichment analysis algorithm (GSEA) known as g:profiler. We next manually downloaded the various GO datasets for our anemia loci. The Human Phenotype Ontology (HPO) database contains ~18,000 unique phenotype terms, and after g:profiler analysis of our 199 anemia loci, we arrived at 619 statistically-significant, unique HPO terms.  To perform TDA on this dataset, we wrote custom python scripts to transform the g:profiler HPO intersections (which contains lists of gene names) for each of 619 HPO terms into a truth table. The resulting output is thus a dataframe of 199 rows (genes) by 619 columns (HPO terms) where each cell tracks 0 for no enrichment and 1 for enrichment.

**STEP 10: Graphically visualize HPO GSEA terms using word clouds**

To crudely visualize a truth table of 199 genes by 619 phenotypes, we employed pythonic wordclouds. We generated wordclouds for the top 40 of 199 genes and also the top 40 of 619 phenotypes found across the entire HPO GSEA set.

**STEP 11: Perform ClinVar TDA analysis**

We next generated allele frequencies on a per-locus basis for each of the 199 anemia-enriched loci. This permits mapping each gene in 4D-space, where each 1 of 4 axes services 1 of 4 ClinVar allele categories (B, LB, LP, and P). We then execute python's keplermapper to generate an interactive HTML that showcases how these 199 anemia loci interact in 4D space based on their ClinVar allelic categorical frequency distributions. A CSV file contains precise gene rosters for each TDA cluster. The output directory contains a log2-transformation of TDA node sizes which is useful for generating custom graphics in BioRender to approximate TDA node relationships.

**STEP 12: Visualize ClinVar TDA output using CIRCOS plots**

CIRCOS plots are useful to showcase the overall connectivity map across all TDA clusters. We used pyCirclize to generate CIRCOS plot digestions of our ClinVar TDA analysis of anemia loci.

**STEP 13: Perform HPO TDA analysis**

This step contains the keplermapper output for the HPO GSEA truth table constructed in step 09. Keplermapper settings are located in the "Anemia199_HPO_TDA_runParameters.txt" file. Precise gene rosters for each TDA cluster is found in the "Anemia199HPOBands.csv" file. An interactive HTML file that visualizes three TDA components (aka tracks) with varying nodes and gene members is provided in the "Anemia199HPOMapper.html" file. The output folder contains log2-transformations of node sizes which is useful for BioRender illustrations to summarize the HPO TDA output. A geneHOOKS CSV file was also produced to help visualize the biology domain expertise.

**STEP 14: Visualize HPO TDA output using CIRCOS plots**

CIRCOS plots are useful to showcase the overall connectivity map across all TDA clusters. We used pyCirclize to generate CIRCOS plot digestions of our HPO TDA analysis of anemia loci.


**STEP 15: Get ClinVar TDA node definitions**

We next extracted all genes found on the entire gene roster for any given TDA cluster (aka node). For each TDA node, we then calculated the mean ClinVar allelic frequency for each ClinVar allelic category. For the ClinVar TDA dataset, we thus generated 6 different barplots that approximates an estimation of how each ClinVar TDA node differs. SVG and PNG versions of the seaborn barplots are found in the output folder.
		
