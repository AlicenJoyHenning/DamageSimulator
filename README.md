<br>

<p align="center">
  <img src="https://github.com/AlicenJoyHenning/damagesimulator/blob/main/inst/extdata/logo.png" alt="limiric_logo" height="140" width="640">
</p>

<br>

Given the lack of ground-truth damaged data, this package was created to simulate damaged cells in single-cell RNA sequencing. It operates by removing a user-defined proportion of cells from an input dataset and applying a damage perturbation function before reintroducing the cells to the original data. Here, perturbation is defined by community-accepted characteristics of damage:
* High mitochondrial content
* Relatively lower cytoplasmic content

<br>

We acknowledge that this understanding of damage is limited and hope to gain more experimentally verified studies characterizing damage in scRNA-seq. But, as the first damage simulator of its kind, we propose ```damagesimulator``` as a valuable starting point for further damage investigations.


<br>

---
