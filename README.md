# ADAM: Adaptive Differential Abundance Method

‘ADAM’ is an R package designed to implement an Adaptive approach for
Differential Abundance Analysis. It provides a streamlined workflow for
analyzing data, selecting the best-performing method, and generating
robust results.

This package is still being actively developed and optimalised for user
experience. \## Installation

To install the `ADAM` package directly from GitHub, follow these steps:

### 1. Ensure you have the `devtools` package installed:

    #install.packages("devtools")

### 2. Install `ADAM` from GitHub and load ADAM package:

    devtools::install_github("lucp9827/ADAM")

    ## Skipping install of 'ADAM' from a github remote, the SHA1 (249d3625) has not changed since last install.
    ##   Use `force = TRUE` to force installation

    library(ADAM)

    ## Registered S3 method overwritten by 'rmutil':
    ##   method         from
    ##   print.response httr

### Dependencies

The `ADAM` package depends on several external packages. These are
automatically installed when you install `ADAM`, but ensure you have the
following available in your R environment: - `dplyr` - `signtrans` -
`phyloseq` - `ALDEx2` - `ANCOMBC` - `DESeq2` - `MicrobiomeStat` -
`RioNorm2` - `corncob` - `microbiome`

#### You can install any missing dependencies manually:

    #install.packages(c("dplyr", "signtrans", "phyloseq"))

#### For Bioconductor packages (e.g., `ALDEx2`):

    #if (!requireNamespace("BiocManager", quietly = TRUE)) {
    #  install.packages("BiocManager") BiocManager::install(c("ALDEx2", "ANCOMBC", "DESeq2"))}

## Using the ADAM\_main Function

The main function of the package, `ADAM_main`, performs the full ADAM
workflow. Here’s a step-by-step guide to use it.

We use example data for scenario 3.3, which has 100 simulated dataset
using the SPSimseq framework. These datasets have 250 taxa, 75 samples
per group with 10% DA taxa which has been selected using a threshold of
LFC=1.5 and in the setting with low sparisty. You can replace this with
your dataset.

    load("simsB_3.3_p.RData")
    data <- sim.data.bulk.p_B[[1]]

Parameters: - `norm_method`: Normalization method (default: `"TSS"`). -
`perc_strata`: Proportion of taxa of the original data in strata
(default: `0.2`). - `OTU_replace`: Logical; Sampling with replacement
(default: `TRUE`). - `threshold_value`: Logical; whether to use adaptive
thresholding (default: `TRUE`). - `alpha`: Significance level for FDR
control (default: `0.05`). - `B`: Number of bootstrap/permutation
datasets (default: `2`). - `path`: File path to save results (default:
current working directory). - `dataset`: Numeric identifier for the
dataset (default: `1`).

Note: Create a Folder in your working directory titled ‘Variations’,
which is used to store the variations of the Test datasets that are
created. Intermediate results (Test dataset, performance of DA methods
on Variations of Test datasets, raw p-values of each test) are saved in
your working directory.

### Running the Function

    library(ADAM)

    # Defining Parameters
    norm_method <- "TSS" 
    perc_strata <- 0.2 
    OTU_replace <- TRUE 
    threshold_value <- TRUE 
    alpha <- 0.05 
    B <- 1 
    path <- getwd() 
    dataset <- 1

    final_results <- ADAM_main( data = data, norm_method = norm_method, perc_strata = perc_strata, OTU_replace = OTU_replace, threshold_value = threshold_value, alpha = alpha, B = B, path = path, dataset = dataset )

    ## Data-driven thresholds were applied based on log fold changes:
    ## 
    ## I1 threshold (high impact): Median log fold change + 3 * MAD =
    ## 1.14680366616798
    ## 
    ## I01 threshold (baseline): Median log fold change =
    ## 0.366795230243319
    ## 
    ## I02 threshold (moderate impact): Median log fold change + 1.5 * MAD =
    ## 0.756799448205649

    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.
    ## Warning in .local(object): Coercing from data.frame class to character matrix 
    ## prior to building taxonomyTable. 
    ## This could introduce artifacts. 
    ## Check your taxonomyTable, or coerce to matrix manually.

    ## aldex.clr: generating Monte-Carlo instances and clr values

    ## conditions vector supplied

    ## operating in serial mode

    ## computing center with all features

    ## aldex.ttest: doing t-test

    ## aldex.effect: calculating effect sizes

    ## Registered S3 methods overwritten by 'registry':
    ##   method               from 
    ##   print.registry_field proxy
    ##   print.registry_entry proxy

    ## 'ancombc' has been fully evolved to 'ancombc2'. 
    ## Explore the enhanced capabilities of our refined method!

    ## Loading required package: foreach

    ## Loading required package: rngtools

    ## 0  features are filtered!
    ## The filtered data has  150  samples and  252  features will be tested!
    ## Pseudo-count approach is used.
    ## Fit linear models ...
    ## Completed.

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

    ## [1] 1
    ## [1] 1
    ## [1] 2
    ##      [,1] [,2] [,3] [,4] [,5] [,6]
    ## [1,]    1   64   96   80   69   67
    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
    ## [1,]   64    1   96   73   69   67   22
    ## [2,]   64    1   96   73   69   67   80
    ## [1] 3
    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7]
    ## [1,]   64    1   96   73   69   67   80
    ##      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
    ## [1,]   80    1   96   73   69   67   64   22
    ## [1] 4

    ## converting counts to integer mode

    ## using pre-existing size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 36 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

    ## [1] 1
    ## [1] 1
    ## [1] 2

The `ADAM_main` function returns a list:

-   `final_res`: Performance metrics (Type I error rate, sensitivity,
    FDR).

-   `method`: The statistical method chosen as best DA method

-   `test`: Raw results of the method applies (statistic, raw pvalues,
    adjusted pvalues (BH-procedure))

### Troubleshooting ———————————————————

-   **Missing Dependencies**: Ensure all required packages are installed
    (see the Dependencies section).

-   **Data Format**: Input data should conform to the expected format
    for phyloseq objects, with ‘group’ used in the sample data
    information to determine group allocation

#### Contributing

Contributions are welcome! Please submit issues or pull requests to
improve the package or its documentation.

This package is still being actively developed and optimalised for user
experience.