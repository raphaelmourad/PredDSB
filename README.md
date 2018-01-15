# Predicting double-strand DNA breaks using epigenome marks or DNA at kilobase resolution


# Data and R scripts to predict double-strand DNA breaks (DSBs) from epigenomic and chromatin data, from DNA motif occurrence data, or from DNA motif occurrence and DNA shape data. 
# Note that DSB data provided are from Lensing et al. (Nature Methods, 2016). Here data only comprise 20.000 DSBs which are a subset of the 84.946 DSBs detected by DSBCapture in order to reduce computational burden for predictions. But, in the article, all DSB sites were used. That explains why prediction performances when using this code are lower than in the article, especially for DNA motif occurrence-based predictions, and DNA motif occurrence and DNA shape-based predictions. 

# There are three main folders: 
# - The folder "data" contains DSB data from  Normal Human Epidermal Keratinocytes (NHEK) cells (subfolder "DSB", file "breakome_DSBcap_hg19_20kseq.bed"). Non-DSB sites are also provided in the file "breakome_DSBcap_hg19_20kseq_neg.bed". The folder also contains the epigenomic and chromatin data (subfolder "Epigenome"), precomputed DNA motif occurrence data (subfolder "Motif") and precomputed DNA shape data (subfolder "Shape"). The folder "data" also contains DSB data (subfolder "DSB_U2OS") and epigenomic and chromatin data (subfolder "Epigenome_U2OS") for Human Bone Osteosarcoma Epithelial (U2OS) cells.
# - The folder "script" contains five R scripts: "comp_shape.R" to compute DNA shape, "create_controlDSB.R" to randomly draw non-DSB sites that are similar to DSB sites, "predictDSB.R" to predict DSBs in NHEK, "predictDSBU2OS.R" to predict DSBs in U2OS and "miscFunctions.R" that includes fonctions for the main script "predictDSB.R". 
# - The folder "results" contains five subfolders: epigenome prediction results (subfolder "predEpigenome"), DNA motif prediction results (subfolder "predMotif"), DNA motif+shape prediction results (subfolder "predMotif+Shape"), NHEK epigenome model-training for U2OS predictions (subfolder "predEpigenomeForU2OS") and U2OS epigenome predictions using NHEK-trained model (subfolder "predEpigenomeU2OS"). 


# Installation of R packages
install.packages(c("pROC","glmnet","ranger","Matrix","gkmSVM"))
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Hsapiens.UCSC.hg19")
biocLite("BSgenome.Hsapiens.UCSC.hg19.masked")
biocLite("JASPAR2016")
biocLite("TFBSTools")
biocLite("DNAshapeR")
biocLite("rtracklayer")
biocLite("Biostrings")

# To predict DSBs in NHEK cells, run the R script predictDSB.R. 

# To predict DSBs in U2OS cells using a NHEK-trained model, run the R script predictDSB.R. 

Reference:
- Stefanie V. Lensing, Giovanni Marsico, Robert Hansel-Hertsch, Enid Y. Lam, David Tannahill, and
Shankar Balasubramanian. DSBCapture: in situ capture and sequencing of DNA breaks. Nature
Methods, 13(10):855–857, August 2016.

Contact:
raphael.mourad@ibcg.biotoul.fr
