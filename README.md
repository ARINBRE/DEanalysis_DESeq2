# DE analysis of integer count RNA-seq data using package DESeq2

## Objective
This educational module demonstrates the steps needed to perform differential expression (DE) analysis using Bioconductor package DESeq2. The analysis applies to integer count matrices generated from RNA-seq data.

## Introduction
Both R code and the same R code within a Jupyter Notebook supported by a conda environment are available. The modeule uses the Pickrell dataset* of sequenced cDNA libraries generated from 69 lymphoblastoid cell lines that were derived from Yoruban Nigerian individuals as part of the HapMap project. The dataset has samples from unrelated male and female subjects. An expression set object that has the assay data and group labels (Gender) is available from Bioconductor package tweeDEseqCountData. The assay data is a numeric matrix of integer counts. These counts were generated by aligning RNA-seq reads to the human genome and counting alignments per genes (used Ensembl gene annotations). The resulting matrix of integer counts represent gene expression profiles for the 69 samples. Bioconductor package DESeq2 accepts integer counts and perform differential expression analysis.

*Pickrell JK, Marioni JC, Pai AA, Degner JF, Engelhardt BE, Nkadori E, Veyrieras JB, Stephens M, Gilad Y, Pritchard JK. Understanding mechanisms underlying human gene expression variation with RNA sequencing. Nature 2010, 464:768-772.doi: 10.1038/nature08872. PMID: 20220758; PMCID: PMC3089435.

## Setup and Installation
The module was conducted using R version 4.4.2 in a Jupyter Notebook with R kernel. The associated Anaconda environment for the Jupyter Notebook is provided in file 'environment_R_4.4.2.yaml'. If you run the R code, ensure the proper version is installed on your machine along with the required libraries. To install Bioconductor packages, packge **BiocManager** is required. You can install it using the following line:

```R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
```

Then use package BiocManager to install the Bioconductor packages needed in this module:

```R
if (!requireNamespace("tweeDEseqCountData", quietly = TRUE)) BiocManager::install("tweeDEseqCountData")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
```

If you prefer to work with Jupyter Notebook, import the yaml environment file into your Anaconda, install R and the needed packages as shown above. Alternatively, you can create your own conda environment. For example:
```
conda create -n R_4.4.2
```
Then activate the environment and install the desired version of R:
```
conda activate R_4.4.2
conda install -c conda-forge r-base=4.4.2 r-essentials -y 
```
To install the R kernel, you can start R in a terminal and use install.packages:
```
conda activate R_4.4.2
R
```
```R
> install.packages("IRkernel")
```

If installing any package directly from R using install.packages("package_name") proves to be problematic due to dependencies that require a specific package version, you may install such packages using conda after activating your environment. For example:
```
conda activate r_env
conda install -c r r-matrixStats=1.5.0
```

 All package versions used in this module are shown at the end of the Jupyter notebook.
 
