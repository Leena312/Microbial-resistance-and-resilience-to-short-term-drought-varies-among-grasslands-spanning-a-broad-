# Microbial-resistance-and-resilience-to-short-term-drought-varies-among-grasslands-spanning-a-broad-
Title: Microbial resistance and resilience to short-term drought varies among grasslands spanning a broad climatic gradient

Authors: Leena Vilonen, Jing Yuan, Melinda Smith, and Pankaj Trivedi

Journal: Currently unpublished


In this repository is the data and code used for the paper titled above. Here are the file names and what they contain:

Philip files:

SGSphylip.phylip -> phylip file for creating 16S tree files for the SGS site. See methods in paper for details.

HYSphylip.phylip -> phylip file for creating 16S tree files for the HYS site. See methods in paper for details.

CHYphylip.phylip -> phylip file for creating 16S tree files for the CHY site. See methods in paper for details.

KNZphylip.phylip -> phylip file for creating 16S tree files for the KNZ site. See methods in paper for details.

R files:

16S-chapter4.R -> All the 16S sequencing data analysis is contained in this R file for effects of drought (resistance to drought)

16S-chapter4-recovery.R -> All the 16S sequencing data analysis is contained in this R file for effects after the drought (resilience to drought)

ITS-chapter4.R -> All the ITS sequencing data analysis is contained in this R file for effects of drought (resistance to drought)

ITS-chapter4-recovery.R -> All the ITS sequencing data analysis is contained in this R file for effects after the drought (resilience to drought)

GHP Nitrogen Paper 1.r -> All the data analysis for qPCR, nitrogen, and enzymes

OTU tables and metadata:

16SMappingFile.txt -> the mapping file for the 16S data. Formatted for mctoolsr.

zotutab_16S.txt.zip -> the zipped OTU table for the 16S data. Formatted for mctoolsr.

ITSMappingFile.txt -> the mapping file the ITS data. Formatted for mctoolsr.

zotutab_ITS.txt -> OTU table for the ITS data. Formatted for mctoolsr.

Data Spreadsheets:

Treatment Labels.xlsx -> this has all the treatment labels for the samples. Important note, two projects worth of data are in this. This data sheet is important for culling the data down for this project as seen in the R files.

qPCR GHP.xlsx -> This has the qPCR values for each sample.

GHPEnzymeActivity.xlsx -> This has all the enzyme values for each sample. See the methods for the specific enzymes.

GHP Inorganic N for R. xlsx -> This has all the inorganic N values for each sample. See the methods for more details.
