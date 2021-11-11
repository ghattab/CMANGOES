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
The example data comes from the following papers:

> **Computer-aided prediction of antigen presenting cell modulators for designing peptide-based vaccine adjuvants**, Nagpal, G., Chaudhary, K., Agrawal, P. et al., J Transl Med (2018), 16 (181).
https://doi.org/10.1186/s12967-018-1560-1.

Dataset was downloaded from https://webs.iiitd.edu.in/raghava/vaxinpad/sequences.phphttps://webs.iiitd.edu.in/raghava/vaxinpad/sequences.php. It is stored at [Data/ace_vaxinpad/](./Data/ace_vaxinpad/).

> **State of the art prediction of HIV-1 protease cleavage sites**, Rögnvaldsson T., You L., Garwicz D., Bioinformatics (2015) 31 (8): 1204-1210. https://doi.org/10.1093/bioinformatics/btu810.

Dataset was downloaded from https://archive.ics.uci.edu/ml/datasets/HIV-1+protease+cleavage. It is stored at [Data/hiv_protease/](./Data/hiv_protease/). 




## Code
|Script|Description|
|---|---|
|[Code/](./Code/)|contains all scripts necessary to run the tool.
|[Code/cmangoes.py](./Code/cmangoes.py)|contains the code that implements the whole pipeline.
|[Code/Dimensionality_Reduction.Rmd](./Code/Dimensionality_Reduction.Rmd)|contains the code used for ML tasks in the paper, specifically the dimensionality reduction step of ML pipeline.
|[Code/Machine_Learning.Rmd](./Code/Machine_Learning.Rmd)|contains the code used for ML tasks in the paper.
|[Code/Preprocessing.Rmd](./Code/Preprocessing.Rmd)|contains the code used for ML tasks in the paper, specifically the preprocessing step of ML pipeline.

## Getting started
```
usage: cmangoes [-h] [--level {1,2,12}] [--image {0,1}] [--show_graph SHOW_GRAPH]
                [--output_path OUTPUT_PATH]
                input_file {b,d} {c,s}

cmangoes: Carbon-based Multi-level Atomic Neighborhood Encodings

positional arguments:
  input_file            A required path-like argument
  {b,d}                 A required character argument that specifies an encoding to be used. b is for
                        binary, d is for discretized
  {c,s}                 A required character argument that specifies a padding to be used. c is for
                        centered, s is for shifted

optional arguments:
  -h, --help            show this help message and exit
  --level {1,2,12}      An optional integer argument that specifies the upper boundary of levels that
                        should be considered. Default: 12 (levels 1 and 2). Example: level 1 returns
                        only first-level neighbors, and level 2 return only second-level neighbors.
                        Level 12 returns level 1 and level 2 neighbors
  --image {0,1}         An optional integer argument that specifies whether images should be created or
                        not. Default: 0 (without images)
  --show_graph SHOW_GRAPH
                        An optional integer argument that specifies whether a graph representation
                        should be created or not. Default: 0 (without representation). The user should
                        provide the number between 1 and the number of sequences in the parsed input
                        file. Example: if number 5 is parsed for this option, a graph representation of
                        the 5th sequence of the input file shall be created and placed in the
                        corresponding images folder
  --output_path OUTPUT_PATH
                        An optional path-like argument. For parsed paths, the directory must exist
                        beforehand. Default: ./CMANGOES_Results
```

## Installation & Running
Place yourself in the [Code](./Code) directory, then run the following command `python cmangoes.py --help` to see how to run the tool. The output of this option is presented in the previous subsection.

You can also use the standalone CMANGOES executable on GNU Linux systems. To do that, place yourself in the [Executable](./Executable) directory and then run the tool with the following command `./cmangoes input_file {b, d} {c, s}` or just use `./cmangoes --help` to see how to run the tool.

## License

Licensed under the MIT License ([LICENSE](./LICENSE)).

### Contribution

Any contribution intentionally submitted for inclusion in the work by you, shall be licensed under the MIT Licence.