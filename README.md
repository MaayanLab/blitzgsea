<img title="a title" alt="blitzGSEA" src="https://github.com/MaayanLab/blitzgsea/blob/main/icon/bgsea_small.png" width=200>

# blitzGSEA Introduction

This Python package provides a computationally performant <b>G</b>ene <b>S</b>et <b>E</b>nrichment <b>A</b>nalysis (GSEA) implementation of the pre-rank algorithm [1]. GSEApy was used as the reference for the running sum and enrichment score calculation [2]. The algorithm estimates the enrichment score (ES) distribution of the null model by fitting data to gamma distibutions instead of calculating permutations for each gene set. blitzGSEA calculates p-values with much higher accuracy than other reference implementations available in Python.

Gene set libraries can directly be loaded from Enrichr (<a href="https://maayanlab.cloud/Enrichr" target="_blank">https://maayanlab.cloud/Enrichr</a>). For this use the `blitzgsea.enrichr.get_library()` function. All libraries can also be listed with `blitzgsea.enrichr.print_libraries()`.

blitzGSEA provides plotting functions to generate publication ready figures similar to the original GSEA-P software. `blitzgsea.plot.running_sum()` plots an enrichment plot for a single gene set and `blitzgsea.plot.top_table()` plots the top `n` gene sets in a compact table. 

# Installation

blitzGSEA is currently only available as a Python package in this GitHub repository. You can install the blitzGSEA Python package and its dependencies through pip by using the following command:

```
$ pip install git+https://github.com/MaayanLab/blitzgsea.git
```

# Run enrichment analysis using blitzGSEA

blitzGSEA depends on two input files. 1) a gene signature and 2) a gene set library. The gene set library is a dictionary with the name of the gene set as key and a list of gene ids as values. Gene set libraries can be loaded directly from Enrichr. The signature should be a pandas dataframe with two columns [0,1]. The first column should contain the gene ids (matching the gene ids in the gene set library).

| index | 0	| 1 |
|:-----|:-------------:|------:|
| 1	| ADO	| -7.833439 |
| 2	| CHUK	| -7.800920 |
| 3	| GOLGA4	| -7.78722 |
| ... | ... | ... |

The gene set library is a dictionary with the gene set names as key and lists of gene ids as values.

```python
{
'ERBB2 SIGNALING PATHWAY (GO:0038128)': ['CDC37',
                                          'PTPN12',
                                          'PIK3CA',
                                          'SOS1',
                                          'CPNE3',
                                          'EGF',
                                          ...
                                         ],
'HETEROTYPIC CELL-CELL ADHESION (GO:0034113)': ['DSC2',
                                                 'ITGAV',
                                                 'ITGAD',
                                                 'LILRB2',
                                                 ...
                                                ],
...
}
```

### Python example

This short example will download two files (signature and gene set library). The gene set library consists of KEGG pathways and the signature is an example signature of differential gene expression of muscle samples from young and old donors. Differential gene expression was computed with Limma Voom.

```python
import blitzgsea as blitz
import urllib.request
import pandas as pd

# download example gene expression signature
url = "https://github.com/MaayanLab/blitzgsea/raw/main/testing/ageing_muscle_gtex.tsv"
urllib.request.urlretrieve(url, "ageing_muscle_gtex.tsv")

# read signature as pandas dataframe
signature = pd.read_csv("ageing_muscle_gtex.tsv")

# list available gene set libraries in Enrichr
blitz.enrichr.print_libraries()

# use enrichr submodule to retrieve gene set library
library = blitz.enrichr.get_library("KEGG_2021_Human")

# run enrichment analysis
result = blitz.gsea(signature, library)
```

### Optional Parameters

The main function of `blitzgsea.gsea()` supports several optional parameters. The default parameters should work well for most use cases. 

| parameter name | type | default	| description |
|:-----|:---------|:-------------|:------|
| `permutations`	| int | 2000	| Number of randomized permutations to estimate ES distributions. |
| `min_size` | int | 5 | Minimum number of genes in geneset. |
| `max_size` | int | Inf | Maximal number of genes in gene set. |
| `anchors`	| int | 20 | Number of gene set size distributions calculated. Remaining are interpolated. |
| `processes`	| int | 4	| Number of parallel threads. Not much gain after 4 threads. |
| `symmetric` | bool | False | Use same distribution parameters for negative and positive ES. If `False` estimate them separately. |
| `signature_cache` | bool | True | Use precomputed anchor parameters. |
| `plotting`| bool | False | Plot estimated anchor parametes. |
| `verbose` | bool | False | Toggle additonal output. |
| `seed` | int | 0 | Random seed. Same seed will result in identical result. If seed equal `-1` generate random seed. |

### Plotting enrichment results

blitzGSEA supports several plotting functions. `blitzgsea.plot.running_sum()` and `blitzgsea.plot.top_table()` can be used after enrichment has been performed. `blitzgsea.plot.running_sum()` shows the running sum of an individual gene set. It has a `compact` mode in which the image will be more readable if small. `blitzgsea.plot.top_table()` shows the top `n` enriched gene sets and displays the results in a table, with normalized enrichment score (NES) and the distribution of hits relative to the gene ranking of the signature.

### Example
```python

import blitzgsea as blitz
import urllib.request
import pandas as pd

# download example gene expression signature
url = "https://github.com/MaayanLab/blitzgsea/raw/main/testing/ageing_muscle_gtex.tsv"
urllib.request.urlretrieve(url, "ageing_muscle_gtex.tsv")

# read signature as pandas dataframe
signature = pd.read_csv("ageing_muscle_gtex.tsv")

# list available gene set libraries in Enrichr
blitz.enrichr.print_libraries()

# use enrichr submodule to retrieve gene set library
library = blitz.enrichr.get_library("KEGG_2021_Human")

# run enrichment analysis
result = blitz.gsea(signature, library)

# plot the enrichment results and save to pdf
fig = blitz.plot.running_sum(signature, "CELL ADHESION MOLECULES", library, result=result, compact=False)
fig.savefig("running_sum.png", bbox_inches='tight')

fig_compact = blitz.plot.running_sum(signature, "CELL ADHESION MOLECULES", library, result=result, compact=True)
fig_compact.savefig("running_sum_compact.png", bbox_inches='tight')

fig_table = blitz.plot.top_table(signature, library, result, n=15)
fig_table.savefig("top_table.png", bbox_inches='tight')

```

The resulting plots will look like the examples below:

#### running_sum.pdf

<div style="bachground-color: white">
<img title="a title" alt="blitzGSEA sunning_sum" src="https://github.com/MaayanLab/blitzgsea/blob/main/icon/running_sum.png" width=300>
</div>

#### running_sum_compact.pdf
<img title="a title" alt="blitzGSEA sunning_sum" src="https://github.com/MaayanLab/blitzgsea/blob/main/icon/running_sum_compact.png" width=300>

#### top_table.pdf
<img title="a title" alt="blitzGSEA sunning_sum" src="https://github.com/MaayanLab/blitzgsea/blob/main/icon/top_table.png" width=300>

### Sample shuffling
This is the sample shuffling algorithm from GSEApy. It performs a t-test to build signatures for phenotype shuffled groups. The input is a gene expression dataframe, which should be normalized for library size. `groups` is a list containing 0 or 1 describing the corresponding group for the samples in `exprs`. The index of `exprs` are the gene ids matching the gene set library. 

```
blitz.shuffle.gsea(exprs, library, groups, permutations=50, seed=1)
```

| parameter name | type | default	| description |
|:-----|:---------|:-------------|:------|
| `exprs`	| pd.DataFrame | NA	| Normalized gene expression matrix. |
| `library` | dictionary | NA | Gene set library. |
| `groups` | list | NA | Phenotype group labels of samples. Labels are 0 or 1. |
| `permutations` | int | 1000 | Number of permutations. |
| `seed`	| int | 1 | Random state seed. |

# Dependencies
Python 3.6+

# References

[1] Subramanian, Aravind, Heidi Kuehn, Joshua Gould, Pablo Tamayo, and Jill P. Mesirov. "GSEA-P: a desktop application for Gene Set Enrichment Analysis." Bioinformatics 23, no. 23 (2007): 3251-3253.

[2] Fang, Z., A. Wolf, Y. Liao, A. McKay, F. Fröhlich, J. Kimmel, L. Xiaohui, and sorrge. "zqfang/gseapy: gseapy-v0. 10.3." (2021).
