<br>

<p align="center">
  <img src="https://github.com/AlicenJoyHenning/damagesimulator/blob/main/inst/extdata/logo.png" alt="limiric_logo" height="140" width="640">
</p>

<br>

Given the lack of ground-truth damaged data, this package was created to simulate damaged cells in single-cell RNA sequencing in an effort to evaluate damaged cell detection strategies.
It operates by removing a user-defined proportion of cells, applying a damage perturbation function, and reintroducing them into the original dataset. Here, perturbation is defined by community-accepted characteristics of damage:
* High mitochondrial content
* Relatively lower cytoplasmic content

<br>

However, we acknowledge that this understanding of damage is limited, as will the simulated cells it generates. We hope to gain experimentally verified studies characterizing damage in scRNA-seq, but until then, as the first damage simulator of its kind, we propose ```damagesimulator``` as a starting point for further damage investigations.


<br>

---
