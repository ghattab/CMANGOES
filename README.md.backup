# CMANGOES

**Carbon-based Multi-level AtomicNeiGhborhOod EncodingS.**


## Manuscript

This package is created for the following paper:

***"CMANGOES: CMANGOES: Carbon-based Multi-level AtomicNeiGhborhOod EncodingS"*** by Georges Hattab, Nils Neumann, Aleksandar Anžel, Dominik Heider

**Paper badge placeholder, link to the PDF placeholder**

Please cite the paper as:
``` Bibtex citation placeholder
```

---
Abstract:

> Exploring new ways to discover organic molecules is critical for developing therapies. With recent advances in bioinformatics, virtual screening of databases is possible. However, biochemical data must be encoded using computer algorithms to make them machine-readable, taking into account distance and similarity measures to support tasks such as similarity search. Motivated by the ubiquity of the carbon element and the structured patterns that emerge in biochemical molecules, we propose a parametric fingerprinting algorithm that integrates neighborhood information at multiple levels. Here we present CMANGOES, carbon-based multi-level atomic neighborhood encodings, a package that integrates and implements a parametric fingerprint. It implements a walk along the carbon chain of an organic molecule to compute different representations of the feature encoding of a molecule. The resulting encodings are task-specific, reproducible, and readily formatted for various domain tasks including machine learning tasks. Encodings created using this package were evaluated using a 10-fold stratified cross validation for binary classification with 8 peptide data sets and 6 different encodings. CMANGOES is built on open-source software and is implemented as a Python package (cmangoes).
The source code is available at [https://github.com/ghattab/CMANGOES](https://github.com/ghattab/cmangoes).

**Paper image placeholder**

## Dependancy

The code is written in Python 3.8.12 and tested on Linux with the following libraries installed:

|Library|Version|
|---|---|
|biopython|1.78|
|ipython|7.27.0|
|matplotlib|3.4.3|
|networkx|2.6.3|
|numpy|1.21.2|
|openbabel|3.1.1|
|pandas|1.3.4|
|pysmiles|1.0.1|



## Data
The data used in the **Example 1** comes from the following paper:

> **Integration of time-series meta-omics data reveals how microbial ecosystems respond to disturbance**, Herold, M., Martínez Arbas, S., Narayanasamy, S. et al. Nat Commun 11, 5281(2020).
https://doi.org/10.1038/s41467-020-19006-2.

It is stored at [Data/cached/example_1/](./Data/cached/example_1) in either a raw format or as a [pickle](https://docs.python.org/3/library/pickle.html) object.

The data used in the **Example 2** comes from the following paper:

> **Short- and Long-Term Transcriptomic Responses of Escherichia coli to Biocides: a Systems Analysis**, Merchel Piovesan Pereira, B., Wang, X., & Tagkopoulos, I. (2020). Applied and environmental microbiology, 86(14), e00708-20.
https://doi.org/10.1128/AEM.00708-20.

It is stored at [Data/cached/example_2/](./Data/cached/example_2) in a raw format.


## Code
|Script|Description|
|---|---|
|[Code/](./Code/)|contains all scripts necessary to run the tool.
|[Code/main.py](./Code/main.py)|contains the code that implements the whole pipeline.
|[Code/Dimensionality_Reduction.Rmd](./Code/Dimensionality_Reduction.Rmd)|contains the code used for ML tasks in the paper, specifically the dimensionality reduction step of ML pipeline.
|[Code/Machine_Learning.Rmd](./Code/Machine_Learning.Rmd)|contains the code used for ML tasks in the paper.
|[Code/Preprocessing.Rmd](./Code/Preprocessing.Rmd)|contains the code used for ML tasks in the paper, specifically the preprocessing step of ML pipeline.

## Getting started


## Installation & Running
### Stable
The easiest way to install the tool is to use our latest Docker image:

```
docker pull aanzel/movis:latest
docker run --publish 8501:8501 --detach --name movis aanzel/movis:latest
```


### Unstable
*Caution! Use at your own risk!*

You could also clone this repository, build a docker container yourself, and run it locally. This is not recommended as we might introduce unstable features that might not end in the next release of MOVIS. Below is a sequence of instructions (for Linux-based systems) to run the **unstable** version of MOVIS:

```
git clone https://github.com/AAnzel/MOVIS.git
cd MOVIS
docker build -t movis-local:unstable .
docker run --publish 8501:8501 --detach --name movis movis-local:unstable
```

## License

Licensed under the MIT License ([LICENSE](./LICENSE)).

### Contribution

Any contribution intentionally submitted for inclusion in the work by you, shall be licensed under the MIT Licence.