# sellandscape
This repository contains a custom script and a small dataset used to decompose dN/dS values into lineage-specific, gene-specific components and residuals

## Dependencies
The script is written using Perl, and was tested on *perl v 5.32.1*. The only side package that is used in decomposition analyses is "Math" to compute CDF. Generally it should be pre-installed with Perl itself, but if it is missing, you can use conda to install it:
```
conda install bioconda::perl-math-cdf
```
## How to use it
The input should be a tab-separated file, where first column corresponds to the lineage ID (ATGC), second column to the gene ID (COG) and third column to dN/dS values. The example run is following:
```
perl decomposition_code_local.pl input_file -clean > output_file
```

Note - the *-clean* flag is optional - it removes outliers based on their residual files and recompute constraints. As a toy example, file *test_data_decomposition.txt* contains data for 20 ATGCs. The following command will decompose their dN/dS values:

```
perl decomposition_code_local.pl test_data_decomposition.txt -clean > test_data_output.out
```

As the optimization is fully determenistic, the script should produce exactly the same results as in *test_data_output.out* regardless of your PC configurations. Hence all results that are deposited https://zenodo.org/records/17575555 can be reproduced by this script.

For any technical questions, please contact: roman.kogay@nih.gov or wolf@ncbi.nlm.nih.gov
