<img title="a title" alt="blitzGSEA" src="https://github.com/MaayanLab/blitzgsea/blob/main/icon/bgsea_small.png" width=200>

# blitzGSEA Introduction

This Python package provides a computationally performant <b>G</b>ene <b>S</b>et <b>E</b>nrichment <b>A</b>nalysis (GSEA) implementation of the pre-rank algorithm [1]. GSEApy was used as the reference for the running sum and enrichment score calculation [2]. The algorithm estimates the enrichment score (ES) distribution of the null model by fitting data to gamma distibutions instead of calculating permutations for each gene set. blitzGSEA calculates p-values with much higher accuracy than other reference implementations available in Python.

Gene set libraries can directly be loaded from Enrichr (<a href="https://maayanlab.cloud/Enrichr" target="_blank">https://maayanlab.cloud/Enrichr</a>). For this use the `blitzgsea.enrichr.get_library()` function. All libraries can also be listed with `blitzgsea.enrichr.print_libraries()`.

blitzGSEA provides plotting functions to generate publication ready figures similar to the original GSEA-P software.

# Installation

blitzGSEA is currently only available as a Python package in this GitHub repository. You can install the blitzGSEA Python package and its dependencies through pip by using the following command:

```
$ pip install git+https://github.com/MaayanLab/blitzgsea.git
```

# Run enrichment analysis using blitzGSEA

blitzGSEA depends on two input files. 1) a gene signature and 2) a gene set library. The gene set library can be any signature with genes and weights associated with them. The signature should be a pandas dataframe with two columns [0,1]. The first column should contain the gene ids (matching the gene ids in the gene set library).

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

# download example GMT file from Enrichr
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
| `anchors`	| int | 20 | Number of gene set size distributions calculated. Remaining are interpolated. |
| `processes`	| int | 4	| Number of parallel threads. Not much gain after 4 threads. |
| `symmetric` | bool | False | Use same distribution parameters for negative and positive ES. If `False` estimate them separately. |
| `plotting`| bool | False | Plot estimated anchor parametes. |
| `verbose` | bool | False | Toggle additonal output. |
| `seed` | int | 0 | Random seed. Same seed will result in identical result. If seed equal `-1` generate random seed. |



# References

[1] Subramanian, Aravind, Heidi Kuehn, Joshua Gould, Pablo Tamayo, and Jill P. Mesirov. "GSEA-P: a desktop application for Gene Set Enrichment Analysis." Bioinformatics 23, no. 23 (2007): 3251-3253.

[2] Fang, Z., A. Wolf, Y. Liao, A. McKay, F. Fröhlich, J. Kimmel, L. Xiaohui, and sorrge. "zqfang/gseapy: gseapy-v0. 10.3." (2021).
