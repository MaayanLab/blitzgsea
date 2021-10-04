<img title="a title" alt="blitzGSEA" src="https://github.com/MaayanLab/blitzgsea/blob/main/icon/bgsea_small.png" width=200>

# blitzGSEA Introduction

This Python package provides a computationally performant <b>G</b>ene <b>S</b>et <b>E</b>nrichment <b>A</b>nalysis (GSEA) implementation of the pre-rank algorithm [1]. The algorithm estimates the enrichment score (ES) distribution of the null model by fitting data to gamma distibutions instead of calculating permutations for each gene set. blitzGSEA calculates p-values with much higher accuracy than other reference implementations available in Python.

blitzGSEA provides plotting functions to generate publication ready figures similar to the oroginal GSEA-P software.

# Installation

blitzGSEA is currently only available as a Python package in this GitHub repository. You can install the blitzGSEA Python package and its depensencies through pip by using the following command:

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

Python

```python
import blitzGSEA as bgsea

result = bgsea.gsea(signature, library)
```
