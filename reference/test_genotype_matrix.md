# Genotype matrix example

Genotype matrix example using a gene-dosage model, which encodes the
SNPs ordinally depending on the genotype allele charge, such as 2 (AA),
1 (AB) and 0 (BB). Valus were drawn from a binomial distribution with
size 2 and probability 0.5.

## Usage

``` r
test_genotype_matrix
```

## Format

### `test_genotype_matrix`

A data frame with 8,539 rows and 30 columns:

- *rownames*:

  SNP ID

- ID1:30:

  Individual's 1 to 30 genotypes; column names correspond to individual
  IDs

## Examples

``` r
summary(test_genotype_matrix)
#>       ID1              ID2             ID3             ID4      
#>  Min.   :0.0000   Min.   :0.000   Min.   :0.000   Min.   :0.00  
#>  1st Qu.:1.0000   1st Qu.:0.000   1st Qu.:0.000   1st Qu.:0.00  
#>  Median :1.0000   Median :1.000   Median :1.000   Median :1.00  
#>  Mean   :0.9966   Mean   :1.001   Mean   :1.001   Mean   :0.99  
#>  3rd Qu.:1.0000   3rd Qu.:2.000   3rd Qu.:2.000   3rd Qu.:1.00  
#>  Max.   :2.0000   Max.   :2.000   Max.   :2.000   Max.   :2.00  
#>       ID5              ID6              ID7             ID8       
#>  Min.   :0.0000   Min.   :0.0000   Min.   :0.000   Min.   :0.000  
#>  1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:1.000   1st Qu.:1.000  
#>  Median :1.0000   Median :1.0000   Median :1.000   Median :1.000  
#>  Mean   :0.9861   Mean   :0.9977   Mean   :1.006   Mean   :1.004  
#>  3rd Qu.:1.0000   3rd Qu.:2.0000   3rd Qu.:2.000   3rd Qu.:1.000  
#>  Max.   :2.0000   Max.   :2.0000   Max.   :2.000   Max.   :2.000  
#>       ID9              ID10            ID11            ID12      
#>  Min.   :0.0000   Min.   :0.000   Min.   :0.000   Min.   :0.000  
#>  1st Qu.:0.0000   1st Qu.:1.000   1st Qu.:1.000   1st Qu.:0.000  
#>  Median :1.0000   Median :1.000   Median :1.000   Median :1.000  
#>  Mean   :0.9916   Mean   :1.014   Mean   :1.003   Mean   :0.998  
#>  3rd Qu.:1.0000   3rd Qu.:2.000   3rd Qu.:2.000   3rd Qu.:1.000  
#>  Max.   :2.0000   Max.   :2.000   Max.   :2.000   Max.   :2.000  
#>       ID13             ID14            ID15            ID16      
#>  Min.   :0.0000   Min.   :0.000   Min.   :0.000   Min.   :0.000  
#>  1st Qu.:0.0000   1st Qu.:1.000   1st Qu.:1.000   1st Qu.:1.000  
#>  Median :1.0000   Median :1.000   Median :1.000   Median :1.000  
#>  Mean   :0.9987   Mean   :1.006   Mean   :1.001   Mean   :1.004  
#>  3rd Qu.:2.0000   3rd Qu.:2.000   3rd Qu.:1.000   3rd Qu.:2.000  
#>  Max.   :2.0000   Max.   :2.000   Max.   :2.000   Max.   :2.000  
#>       ID17            ID18            ID19             ID20      
#>  Min.   :0.000   Min.   :0.000   Min.   :0.0000   Min.   :0.000  
#>  1st Qu.:1.000   1st Qu.:1.000   1st Qu.:0.0000   1st Qu.:0.000  
#>  Median :1.000   Median :1.000   Median :1.0000   Median :1.000  
#>  Mean   :1.004   Mean   :1.008   Mean   :0.9967   Mean   :1.002  
#>  3rd Qu.:2.000   3rd Qu.:2.000   3rd Qu.:2.0000   3rd Qu.:2.000  
#>  Max.   :2.000   Max.   :2.000   Max.   :2.0000   Max.   :2.000  
#>       ID21            ID22             ID23             ID24       
#>  Min.   :0.000   Min.   :0.0000   Min.   :0.0000   Min.   :0.0000  
#>  1st Qu.:1.000   1st Qu.:0.0000   1st Qu.:0.0000   1st Qu.:0.0000  
#>  Median :1.000   Median :1.0000   Median :1.0000   Median :1.0000  
#>  Mean   :1.002   Mean   :0.9938   Mean   :0.9855   Mean   :0.9961  
#>  3rd Qu.:1.000   3rd Qu.:1.0000   3rd Qu.:1.0000   3rd Qu.:1.0000  
#>  Max.   :2.000   Max.   :2.0000   Max.   :2.0000   Max.   :2.0000  
#>       ID25            ID26            ID27            ID28      
#>  Min.   :0.000   Min.   :0.000   Min.   :0.000   Min.   :0.000  
#>  1st Qu.:0.000   1st Qu.:1.000   1st Qu.:1.000   1st Qu.:1.000  
#>  Median :1.000   Median :1.000   Median :1.000   Median :1.000  
#>  Mean   :1.007   Mean   :1.007   Mean   :1.007   Mean   :1.013  
#>  3rd Qu.:2.000   3rd Qu.:1.000   3rd Qu.:2.000   3rd Qu.:2.000  
#>  Max.   :2.000   Max.   :2.000   Max.   :2.000   Max.   :2.000  
#>       ID29            ID30      
#>  Min.   :0.000   Min.   :0.000  
#>  1st Qu.:0.000   1st Qu.:1.000  
#>  Median :1.000   Median :1.000  
#>  Mean   :1.002   Mean   :1.012  
#>  3rd Qu.:2.000   3rd Qu.:2.000  
#>  Max.   :2.000   Max.   :2.000  
```
