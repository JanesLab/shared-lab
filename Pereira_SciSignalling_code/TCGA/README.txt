TCGA-BRCA.htseq_fpkm.tsv is fpkm data for all TCGA breast cancer cases downloaded directly from Xena (https://xena.ucsc.edu/)

TCGA-BRCA.htseq_fpkm_geneid.csv is the same as above matrix but with associated Gene IDs (matched from ENSG IDs)

FPKM_TNBC_path.csv has FPKM values for all 122 TNBC cases determined by pathology report (ER, PR negative and HER2 not amplified using info in BRCA_clinicalMatrix, Sham sent me this clinical matrix)

Code in TCGA_TNBCcases.R goes through steps to produce all these files