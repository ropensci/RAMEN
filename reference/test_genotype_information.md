# Genotype metadata example data set

Genotype chromosome 21 position metadata obtained from a private
genotyping data set.

## Usage

``` r
test_genotype_information
```

## Format

### `test_genotype_information`

A data frame with 8539 rows and 3 columns:

- CHROM:

  Chromosome

- POS:

  Probe genomic position (h19)

- ID:

  SNP ID

## Examples

``` r
summary(test_genotype_information)
#>        CHROM           POS                   ID      
#>  Length   :8539   Min.   :10873592   Length   :8539  
#>  N.unique :   1   1st Qu.:20974287   N.unique :8539  
#>  N.blank  :   0   Median :25404407   N.blank  :   0  
#>  Min.nchar:   5   Mean   :25550382   Min.nchar:  15  
#>  Max.nchar:   5   3rd Qu.:30278092   Max.nchar:  56  
#>                   Max.   :35619929                   
```
