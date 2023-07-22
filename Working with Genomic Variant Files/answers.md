---
title: "answers"
author: "Bonface Onyango"
date: "2023-07-21"
output: html_document
---

**Q1: How many positions are found in this region in the VCF file?**

A: The file contains 69 positions

**Q2: How many samples are included in the VCF file?**

A: the samples included are 3489

**Q4: How many positions are there with `AC=1`? Note that you cannot simply count lines since the output of `bcftools filter` includes the VCF header lines. You will need to use `bcftools query` to get this number.**

A: The positions with allele count of one are 1075 positions

**Q5: What is the ratio of transitions to transversions (ts/tv) in this file?**

A: The ratio of transversion to transition is 3.47

## Mutation Annotation Format (MAF) files

#### Install necessary libraries

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)


```

```{r}
install.packages("BiocManager")  # Install BiocManager package (if not installed)
BiocManager::install("maftools")  #
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

```{Install the maftools package from Bioconductor}

```

```{r}
# sample summary
getSampleSummary(lgg)

```

**Q6: What is the median number of variants per sample in this data set?**

A: From the summary statistics, the median number of variants per sample is 28

```{r}
plotmafSummary(maf = lgg, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

```

# Generating plots with maftools

An oncoplot of the top five mutated genes

```{r}


#oncoplot for top five mutated genes.
oncoplot(maf = lgg, top = 5)

```

A boxplot of the transistion-to-transversion ratio

```{r}
laml.titv = titv(maf = lgg, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = laml.titv)
```

A plot comparing the mutation load in this LGG cohort to other TCGA cohorts. Use log scale.

```{r}
laml.mutload = tcgaCompare(maf = lgg, cohortName = 'Example-LAML', logscale = TRUE, capture_size = 50)

```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{r cars}
summary(cars)
```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}
plot(pressure)
```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
