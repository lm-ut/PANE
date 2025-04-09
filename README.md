# PANE - Principal component Ancestry proportions using NNLS Estimation

*PANE* ('Pa-Neh' - bread in Italian)  leverages PCA and NNLS (Non-Negative Least Squares) to assess the ancestral composition of admixed individuals with high accuracy and reliability.   


In this page you will find instructions for a basic usage, for a more in-depth tutorial click [here](https://lm-ut.github.io/PANE/articles/Tutorial.html)

## Installation

```
install.packages("devtools") # if you have not installed "devtools" package
devtools::install_github("lm-ut/PANE")
library("PANE")
```

*PANE* requires the following R packages: 

```
install.packages("dplyr")
install.packages("nnls")
install.packages("ggplot2")
```

## Introduction

*PANE* estimations are based on a PCA where the target and source groups are available. Along with those, we suggest to use additional groups to better define the PC space. Once that the PC space is defined, a set of NNLS is then applied on the PC coordinates, effectively summarizing the genetic ancestry. 

For its most basic usage PANE needs:  

- A PCA matrix (a dataframe with N PCs)
- A list of target and source groups (or samples)

### Basic Usage Example
  
```
$ pca = read_eigen(pca_input = 'data/TOY.pca.evec')
$ PANE_result = pane(pca_input = pca, sources = c('GST1','GST2'), admixed=c('70GST1.30GST2'))
```


#### Reading the PCA Matrix


*PANE* package has two functions to read the PCA matrix, ```read_eigen()``` and ```read_flash()```.  
* ```read_eigen()``` will read a PCA that has been created with smartpca from the [EIGENSOFT](https://github.com/DReichLab/EIG) software.  
* ```read_flash()``` instead, will read a PCA that has been created with [flashpca](https://github.com/gabraham/flashpca) software.
   
The goal of both functions is to set the PCA file as follows:  
  
| POP | IND   | PCN |
| --- | ----  | ------- |
| SRC1 | SRC1_1 | 0.01 |
| SRC1 | SRC1_2 | 0.02 |
| ADM1 | ADM1_1 | 0.08 |
  
If neither ```read_eigen()``` nor ```read_flash()``` is for you, you might want to simply use ```read.table()```, and set the file so that it has the aforementioned look.  


#### Running PANE

  
The function ```pane()``` requires also a list of the target and reference groups: you can provide the list in two ways.

1) You can provide an 'AS_file': a file with the list of the Admixed groups (A) and the Source groups (S). If you want to use PANE sample-wise rather than group-wise, simply adjust the PCA file so that in the 'POP' column is identical to the 'IND' column, and set the AS_file with the samples list, rather than the group list.   
The AS_file is a two-columns file with the population list on the first column, and the 'A/S' information on the second column. The 'A/S' information stands for Admixed (A) or Source (S). For each population/group we will indicate whether PANE should consider it as a Source (S) or as an admixed target (A), the file looks like this:  

| POP | A/S |
| --- | --- |
| SRC1 | S |
| SRC2 | S |
| ADM1 | A |

To read the AS_file, a simple ```read_table(file, header=T)``` will be sufficient.   

```
$ AS_file = read.table('data/Example_AS', header=TRUE)
```

With the PCA and AS_file loaded, we are finally ready to run PANE as follows:

```
$ PANE_result = pane(pca_input = pca, as_file = AS_file)
```
  
2) You can avoid relying on the AS_file if you wish, using a vector of the target and source groups directly in pane() function, as follows:
  
```
$ PANE_result = pane(pca_input = pca, sources = c('GST1','GST2'), admixed=c('70GST1.30GST2'))
```

#### Writing PANE output
  
Finally, if you want to save *PANE* results on a table-like format, you can use ```write_pane()```:
  
```
$ PANE_result <- pane(pca_input = pca, sources = c('GST1','GST2'), admixed=c('70GST1.30GST2'))
$ write_pane(PANE_result, output_name = 'my_dir/my_pane_results.txt')
```

## Cite PANE

If you use PANE, please cite [our paper](https://doi.org/10.1186/s13059-025-03491-z)

## Contact

For questions and bug reports please contact [LM](mailto:ludovica.molinaro@kuleuven.be).
