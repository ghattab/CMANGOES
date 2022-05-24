# CMANGOES

**Carbon-based Multi-level AtomicNeiGhborhOod EncodingS.**


## Manuscript

This package is created for the following paper:

***"A parametric approach to molecular encodings of carbon-based multilevel atomic neighborhoods"*** by Georges Hattab, Nils Neumann, Aleksandar Anžel, Dominik Heider

**Paper badge placeholder, link to the PDF placeholder**

Please cite the paper as:
``` Bibtex citation placeholder
```

---
Abstract:

> TBD

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
The example data we used to evaluate CMANGOES comes from the [peptidereactor](https://doi.org/10.1093/nargab/lqab039) [repository](https://github.com/spaenigs/peptidereactor/tree/master/data). We selected 30 data sets from that repository and placed them at [Data/Original_datasets/](Data/Original_datasets/). Each data set has a separate README file that contains the additional information of that data set.

After encoding those data sets with CMANGOES, we used the [peptidereactor](https://doi.org/10.1093/nargab/lqab039) to evaluate our method against all other encoding methods available in the peptidereactor. Evaluation results of CMANGOES are placed at [Data/Visualization_data/data/](Data/Visualization_data/data/). Evaluation results of other encoding methods provided by the peptidereactor are located at [Data/Visualization_data/peptidereactor_vis_data/](Data/Visualization_data/peptidereactor_vis_data/). We combined both evaluation results to create visualizations present in our manuscript.



## Code
|Script|Description|
|---|---|
|[Code/](./Code/)|contains all scripts necessary to run the tool.
|[Code/cmangoes.py](./Code/cmangoes.py)|contains the code that implements the whole pipeline.
|[Code/batch_encoding.py](./Code/batch_encoding.py)|contains the code that uses CMANGOES to batch-encode selected data sets from the [peptidereactor](https://doi.org/10.1093/nargab/lqab039) [repository](https://github.com/spaenigs/peptidereactor/tree/master/data).
|[Code/visualize.ipynb](./Code/visualize.ipynb)|contains the code that creates visualizations of the final results.
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
