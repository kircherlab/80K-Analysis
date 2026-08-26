## Massively parallel characterization and predictive modelling of neuronal regulatory variation

Kilian Salomon1,\*,#, Chengyu Deng2,3\*, Pyaree Mohan Dash 1, Theofilos Chalkiadakis1, Qinrui Li2,3, Ziwei Chen4, Nicholas F. Page2,3, Mustafa Helal5, Sebastian Röner6,8, Anshul Kundaje4,7, Claudia Langenberg6,8, Maik Pietzner8,9, Jay Shendure10,11,12,13, Max Schubach1, Nadav Ahituv2,3,# and Martin Kircher1,5,#

\* These authors contributed equally
\# correspondence should be addressed to kilian.salomon@bih-charite.de, nadav.ahituv@ucsf.edu, martin.kircher@uni-luebeck.de

1.	Computational Genome Biology, Exploratory Diagnostic Sciences, Berlin Institute of Health at Charité - Universitätsmedizin Berlin, Berlin, Germany
2.	Department of Bioengineering and Therapeutic Sciences, University of California San Francisco, San Francisco, CA, USA
3.	Institute for Human Genetics, University of California San Francisco, San Francisco, CA, USA
4.	Department of Computer Science, School of Engineering, Stanford University
5.	University of Luebeck, Institute of Human Genetics, University Hospital Schleswig-Holstein, Campus Luebeck, Germany
6.	Computational Medicine, Medical and Health Data Sciences, Berlin Institute of Health at Charité - Universitätsmedizin Berlin, Berlin, Germany
7.	Department of Genetics, School of Medicine, Stanford University
8.	Precision Healthcare University Research Institute, Queen Mary University of London, London, UK.
9.	Health Data Modelling, Medical and Health Data Sciences, Berlin Institute of Health at Charité - Universitätsmedizin Berlin, Berlin, Germany
10.	Department of Genome Sciences, University of Washington, Seattle, WA, USA
11.	Seattle Hub for Synthetic Biology, Seattle, WA, USA
12.	Brotman Baty Institute for Precision Medicine, Seattle, WA, USA
13.	Howard Hughes Medical Institute, Seattle, WA, USA
 
### Abstract
Disease-associated variants reside frequently in noncoding cis-regulatory elements (CREs), yet their functional consequences remain poorly understood. We performed a large-scale lentiMPRA in human excitatory neurons, quantifying the impact of >46,000 naturally occurring variants across >27,000 candidate CREs near 524 disease-associated genes. These data improved regulatory variant effect predictions beyond state-of-the-art models. Significant allelic effects occurred at comparable rates across common, rare, and singleton variants, demonstrating that, within MPRA-measurable effects, population frequency carries limited information about per-variant regulatory impact. Variant effect detectability and magnitude were governed primarily by baseline activity of the enclosing regulatory element and local sequence context. Regulatory effects were distributed across numerous transcription factors rather than concentrated in master regulators, consistent with a combinatorial enhancer architecture. We establish a large-scale functional variant catalog and provide a complementary benchmark and resource for developing and evaluating models of noncoding regulatory variation.

### Repository Contents

This repository contains the analysis and visualization code used to generate the figures in the manuscript. Code is organized by analysis type rather than by figure number:

- **`cCRE_based_analysis/`**: Design overview of the MPRA (gene/cCRE selection, control sets) and analysis of cCRE regulatory activity, including brain/tissue cCRE enrichment (Fig. 1, Fig. 2).
- **`variant_based_analysis/`**: Allelic variant effect analysis across the allele-frequency spectrum, transcription factor binding site (TFBS) enrichment, and variant locus visualization (Fig. 3, Fig. 4).
- **`modeling/`**: Benchmarking of sequence-based model predictions (e.g. Enformer, AlphaGenome, CADD, ChromBPNet) against measured MPRA variant effects (Fig. 5). See `modeling/README.md` for details on the model prediction sources.
- **`00_helpful_functions/`**: Shared utility functions imported across the analysis notebooks.

Each analysis subfolder contains the Jupyter notebook(s) used to produce the corresponding figures, together with a local `data/` directory holding the processed input tables the notebooks read from and a `global80K_config.yaml` file resolving named file paths used within the notebooks.
MPRA datasets generated in this study have been submitted to the IGVF portal (http://data.igvf.org) under the following accession number: IGVFDS1419ZPHD. The ATAC-Seq performed in NGN2-derived excitatory neurons used in this study is available on request. The visualization of open-chromatin peaks is available at the Gene Expression Omnibus under accession number GSE113480.

### Environment setup

Three separate conda environments cover the code in this repository (see `environment/`):

| Environment | Covers | Create with |
|---|---|---|
| `mpra80k_python` | All Python notebooks/scripts (`cCRE_based_analysis/`, `variant_based_analysis/`, `modeling/`, `00_helpful_functions/`) | `conda env create -f environment/python_environment.yml` |
| `mpra80k_element_variant_processing` | `perform_bcalm_elements_cli.R` / `perform_bcalm_variants_cli.R` (BCalm/mpralm element- and variant-level activity calls) | `conda env create -f environment/element_variant_processing_environment.yml` |
| `mpra80k_variant_plotting` | `variant_based_analysis/variant_loci_visualization/` Gviz locus plots | `conda env create -f environment/variant_plotting_environment.yml` |

The BCalm environment additionally requires cloning [kircherlab/BCalm](https://github.com/kircherlab/BCalm) (v0.9.0 used here) and passing its path as the `<bcalm_path>` CLI argument which is loaded at runtime via `devtools::load_all()`, not installed as a package. External data dependencies too large to bundle in this repo (>25MB) are listed in `bigger_than_25mb.tsv`; files referenced from an absolute local path are listed in `not_accessible.tsv`.
