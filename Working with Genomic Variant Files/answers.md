---
title: "Untitled"
author: "Bonface Onyango"
date: "2023-07-21"
output: ''
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

**Q1: How many positions are found in this region in the VCF file?**

A: The file contains 69 positions

**Q2: How many samples are included in the VCF file?**

A: the samples included are 3489

**Q4: How many positions are there with `AC=1`? Note that you cannot simply count lines since the output of `bcftools filter` includes the VCF header lines. You will need to use `bcftools query` to get this number.**

A: The positions with allele count of one are 1075 positions

**Q5: What is the ratio of transitions to transversions (ts/tv) in this file?**

A: The ratio of transversion to transition is 3.47

# Install necessary libraries

```{r}
install.packages("BiocManager")  # Install BiocManager package (if not installed)
BiocManager::install("maftools")  # Install the maftools package from Bioconductor

```

```{r}
BiocManager::install("maftools")  # Install the maftools package from Bioconductor

```

```{r}
library(maftools)

# view available TCGA cohorts, check that LGG is there
tcgaAvailable()

# load the available LGG cohort
lgg <- tcgaLoad(study = "LGG")

# view a summary of this file
lgg


```

```{r}

```

```         
```

## Mutation Annotation Format (MAF) files

```{r}
# Read the MAF file into maftools

# Calculate the number of variants per sample
variants_per_sample <- countMutations(lgg, "Tumor_Sample_Barcode")

```

```{r}
# Load the required libraries
library(maftools)

# Load the LGG cohort
lgg <- tcgaLoad(study = "LGG")

# Oncoplot of the top five mutated genes
oncoplot(lgg, top = 5)

# Boxplot of the transition-to-transversion ratio
boxplot_TransitionTransversion(lgg)

# Plot comparing mutation load in LGG cohort to other TCGA cohorts (log scale)
tcgaCompare(maflite = lgg, cohort = "LGG", logscale = TRUE)
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
