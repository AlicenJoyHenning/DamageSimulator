![DamageSimulator](docs/damagesimulator.svg)

<br>

## Context

Given the limited availability of ground-truth-damaged data, this package was created to simulate damaged cells in single-cell RNA sequencing. It operates by removing a user-defined proportion of cells from an input dataset and applying a damage perturbation function before reintroducing the cells to the original data. Here, perturbation is defined by community-accepted characteristics of damage:
* High mitochondrial content
* Relatively lower cytoplasmic content

<br>

## Quickstart

Install using the `devtools` package in `R`, 
```R
devtools::install_github("AlicenJoyHenning/damagesimulator", build_vignettes = TRUE)
```

Use your own count matrix or use those provided interanlly through `scRNAseq`, 
```R
library(scRNAseq)

# Retrieve dataset & extract counts 
sce <- fetchDataset(name = "he-organs-2020", version = "2023-12-21", path = "blood")
counts <- sce@assays@data$counts
rownames(counts) <- rownames(sce)
colnames(counts) <- colnames(sce)
counts <- as(counts, "sparseMatrix") 
count_matrix <- counts
```

And run the `simulate_damage` function specifying your desired proportion of damage and which damage population you are targeting, 
```R
results <- simulate_damage(count_matrix = counts,
                           damage_proportion = 0.2,  # 20 % of input cells will be damaged by the function 
                           beta_proportion = 1)      # 100 % of the damaged cells will exhibit an advanced profile of damage 
```


<br>

## Concluding remarks

We acknowledge that this understanding of damage is limited and hope to gain more experimentally verified studies characterizing damage in scRNA-seq. But, as the first damage simulator of its kind, we propose ```damagesimulator``` as a valuable starting point for further damage investigations.

<br>

---
