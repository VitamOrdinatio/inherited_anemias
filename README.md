HELLO

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

Across all 199 anemia-enriched genetic loci, we scraped 192,296 total alleles from each of six NCBI ClinVar categories (i.e., CC, B, LB, US, LP, and P). We performed multivariate statistical analyses, including pairplots, t-SNE, linear regressions, correlation coefficients, principal component analysis (PCA), PCA explaind variance ratios, unfiltered k-means clustering, and two-component (2PC) k-means clustering.


