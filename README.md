# Studying 3D genome evolution using genomic sequence

**Overview**

We propose a novel approach to study the 3D genome evolution in vertebrates using the genomic sequence only, e.g. without the need for Hi-C data. The approach is simple and relies on comparing the distances between convergent and divergent CTCF motifs (ratio R).

**Systems Requirements**

The scripts were written in R language. 

To run the scripts, you need several R packages. To install the packages:
`install.packages(c("circlize","bootstrap","phytools"))` \
`source("https://bioconductor.org/biocLite.R")` \
`biocLite("GenomicRanges")` 

**Usage**

Before computing ratio R, one must use FIMO to scan for CTCF motifs. FIMO outputs a "tsv" file that is used by the "compDistFun.R" function to compute R. You can use FIMO online (http://meme-suite.org/tools/fimo), select the genome assembly and upload the CTCF MEME file "data/CTCF_meme/MA0139.1.meme".

You can also use precomputed CTCF motifs from FIMO for selected assemblies: hg19, hg38, bosTau8, ce11, dm6, mm10, rn6 and xenTro7. These motifs are available in the folder "data/CTCF_motif".

In this package, there are three main folders: 
- The folder "data" contains: the CTCF motif PWM in MEME format for FIMO (subfolder "CTCF_meme"), CTCF motifs called by FIMO (subfolder "CTCF_motif"), CTCF motif DeepBind scores (subfolder "CTCF_deepbind"), CTCF motif conservation scores (subfolder "CTCF_cons"), CTCF ChIP-seq peaks from GM12878 ENCODE cells (subfolder "CTCF_peak"), 3D domain borders (subfolder "3DDomainBorder"), chromosome regions (subfolder "region"), phylogenetic tree (subfolder "tree").
- The folder "script" contains seven R scripts: "compDistFun.R" is the function to compute R, "compute_R_species.R" to compute R in different species, "compute_R_Rp_Rc_human.R" to compute different R's in human, "compute_R_regions_human.R" to compute R for different chromosome regions, "compute_R_peak_domainBorders_human.R" to compute R at CTCF peaks and 3D domain borders, "test_differenceR_2species.R" to test the difference of R values between two species and "ancestral_R_reconstruction.R" to use ancestral R reconstruction. 
- The folder "results" contains three subfolders: precomputed R for different species (subfolder "matPvalPM"), DNA motif prediction results (subfolder "predMotif"), precomputed distances between consecutive motifs depending on orienation (subfolder "matPM") and ancestral R reconstruction results (subfolder "phylo"). 

**References**

**Contact**:
raphael.mourad@ibcg.biotoul.fr
