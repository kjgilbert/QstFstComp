# QstFstComp

[![Build Status](https://travis-ci.org/kjgilbert/QstFstComp.png?branch=master)](https://travis-ci.org/kjgilbert/QstFstComp)

An R function to compare the *Q<sub>ST</sub>* of a single phenotypic trait to the mean *F<sub>ST</sub>* of series of marker loci. `QstFstComp` calculates the distribution of *Q<sub>ST</sub>* – *F<sub>ST</sub>* under a model assuming neutrality of both the phenotypic trait and the genetic markers from which *F<sub>ST</sub>* is estimated.

[Gilbert and Whitlock (*In Press*)](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12303/abstract) describes the use of these procedures and their derivation. If you use this method, please cite:

	Gilbert KJ and MC Whitlock (*In Press*) *Q<sub>ST</sub>* *F<sub>ST</sub>* comparisons with unbalanced half-sib designs. *Molecular Ecology Resources*.

The method (and this code) is based on [Whitlock and Guillaume (2009)](http://www.genetics.org/content/183/3/1055) which may also be cited.


The package can be installed using [devtools](https://github.com/hadley/devtools), which itself can be installed from CRAN with

```
install.packages("devtools")
```

Once devtools is installed, run

```
library(devtools)
install_github("kjgilbert/QstFstComp")
library(QstFstComp)
```
and the package will be installed and open. Functions for calculating *F<sub>ST</sub>* according to [Weir and Cockerham (1984)](http://www.jstor.org/discover/10.2307/2408641?uid=2&uid=4&sid=21104217684983) depend on the package [`hierfstat`](http://cran.r-project.org/web/packages/hierfstat/index.html) by Jerome Goudet.
<!-- or can alternatively be installed with devtools::install github("kjgilbert/QstFstComp") where there should be an underscore between "install" and "github"-->

There is only one function to run the analysis from the package: `QstFstComp`.  See the help page `?QstFstComp` for more information once installed and loaded.


## Documentation

The procedure requires two input data files described in the following sections and can then be run with the steps outlined below. 

#### Step 1: *F<sub>ST</sub>* INPUT FILE
Data from genetic markers are expected to be in a .csv format, with the first row giving column names.

*F<sub>ST</sub>* is calculated in two different ways, and the user must choose which method to use. The default is the [Weir and Cockerham (1984)](http://www.jstor.org/discover/10.2307/2408641?uid=2&uid=4&sid=21104217684983) method (`AFLP = FALSE`), as first implemented in [Whitlock and Guillaume (2009)](http://www.genetics.org/content/183/3/1055). If using this method, the genetic data (fst.dat) must be a data frame where each row is an individual and the first column indicates population of origin and the following columns represent genotypes at loci (as in hierfstat package). 

See the the example data file for multiallelic loci: 
```
data(gtrunchier) # example dataset from hierfstat 
gtrunchier[,-1]  # contains one column extraneous to this analysis
``` 
or for biallelic loci:
```
data(biallelic)
```

If the genetic marker data are from AFLP markers, the user must specify `AFLP = TRUE` as an input parameter, and the input format for data must be a .csv file with populations in columns and loci in rows. The entries of this file should be ![q_hat](https://github.com/kjgilbert/QstFstComp/raw/master/q_hat.png) values as calculated in [Lynch and Milligan (1994)](http://www.indiana.edu/~lynchlab/PDF/Lynch63.pdf). The first row should be a header with population identifiers, and the first column should be locus identifiers followed by columns with the ![q_hat](https://github.com/kjgilbert/QstFstComp/raw/master/q_hat.png) values for each locus. To the right of these columns, add another set of columns with the ![q_hat](https://github.com/kjgilbert/QstFstComp/raw/master/q_hat.png) variances for each locus, in the same population order as the ![q_hat](https://github.com/kjgilbert/QstFstComp/raw/master/q_hat.png) values. A table of such ![q_hat](https://github.com/kjgilbert/QstFstComp/raw/master/q_hat.png) values and ![q_hat](https://github.com/kjgilbert/QstFstComp/raw/master/q_hat.png) variances can be obtained from the program [AFLP-SURV](http://www.ulb.ac.be/sciences/lagev/aflp-surv.html), by Xavier Vekemans. 

See `data(aflp)` as an example.

#### Step 2: *Q<sub>ST</sub>* INPUT FILE

Trait data for the estimation of *Q<sub>ST</sub>* (qst.dat) should be in a .csv file with column names in the first row. Different formatting is required depending on whether the user wishes to analyze data sampled from a half-sib dam model (dams nested within populations with all offspring coming from separate sires) or a half-sib sire model (dams nested within sires, sires nested within populations, with offspring coming from each dam that shares a sire across other dams).

For the half-sib dam model (`breeding.design = "half.sib.dam"`), there must be 3 columns of data where each row after the header contains the trait data for one individual and its identifiers for the population and dam it originated from:
- column 1 identifies the population of origin. 	Each population should have a unique name or number.
- column 2 identifies the dam (mother) of the offspring. Each dam should have a unique name or number (i.e., a dam from one population should never have the same name as a dam from another population).
- column 3 is the numerical value of the trait in question for the individual identified.

See the example data file:
```
data(hsdam)
```

For the half-sib sire model (`breeding.design = "half.sib.sire"`), there must be 4 columns of data where each row after the header contains the data for one individual: its identifiers for the population, sire, and dam it originated from, and the value of its trait:
- column 1 identifies the population of origin. Each population should have a unique name or number.
- column 2 identifies the sire (father) of the offspring. Each sire should have a unique name.
- column 3 identifies the dam (mother) of the offspring. Each dam should have a unique name.
- column 4 is the numerical value of the trait in question for the individual identified.

See the example data file:
```
data(hssire)
```


#### Step 3: Running the analysis

(1) Load the file containing information about the genetic marker frequencies to an object called  “MarkerInfo”

```
MarkerInfo <- read.csv(“/Your_FST_FileName.csv”) 
```

Replace the text in quotes with the directory path for your *F<sub>ST</sub>* file.

(2) Load information about the breeding experiment and the trait to a data frame called “TraitInfo"

```
TraitInfo <- read.csv(“/Your_QST_FileName.csv”) 
```

Replace the text in quotes with the directory path for your *Q<sub>ST</sub>* file.

(3) Run the function “QstFstComp” using these two data objects: 

```
QstFstComp(fst.dat = MarkerInfo, qst.dat = TraitInfo, numpops = XXX, breeding.design = "half.sib.dam", nsim = 10000)
```

Replace `XXX` with the number of populations included in your study. (This should be the same number for both the *F<sub>ST</sub>* and *Q<sub>ST</sub>* data sets.) Replace `"half.sib.dam”` with `“half.sib.sire”` if you are using the half-sib sire model. If using the half-sib dam model, but with different relatedness of offspring, this can be accommodated by changing the default parameter, `dam.offspring.relatedness=0.25` from the default of 1/4 to the desired value.

Add a parameter `output = “full”` to have the function produce a longer list of output values. See the help page `?QstFstComp` for more details on output values and example analyses. The distribution of *Q<sub>ST</sub>*-*F<sub>ST</sub>* can be plotted from the vector of values which is automatically output to a text file in the current working directory.

