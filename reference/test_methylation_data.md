# Methylation data matrix example

Simulated M values of the 3000 probes selected in *test_array_manifest*
for 30 individuals. Values were converted from beta values, which were
drawn from a bimodal Beta distribution.

## Usage

``` r
test_methylation_data
```

## Format

### `test_methylation_data`

A data frame with 3,000 rows and 30 columns:

- *rownames*:

  Probe IDs (column TargetID in the EPIC array)

- ID1:30:

  DNAme profile of individuals 1 to 30; column names correspond to
  individual IDs

## Examples

``` r
summary(test_methylation_data)
#>       ID1              ID2              ID3              ID4        
#>  Min.   :-1.895   Min.   :-1.708   Min.   :-1.626   Min.   :-3.679  
#>  1st Qu.: 1.673   1st Qu.: 1.691   1st Qu.: 1.624   1st Qu.: 1.602  
#>  Median : 2.766   Median : 2.724   Median : 2.773   Median : 2.716  
#>  Mean   : 3.030   Mean   : 3.011   Mean   : 3.026   Mean   : 2.957  
#>  3rd Qu.: 4.098   3rd Qu.: 4.125   3rd Qu.: 4.109   3rd Qu.: 4.042  
#>  Max.   :13.968   Max.   :15.591   Max.   :16.356   Max.   :13.406  
#>       ID5              ID6              ID7              ID8        
#>  Min.   :-2.359   Min.   :-1.954   Min.   :-2.066   Min.   :-2.088  
#>  1st Qu.: 1.577   1st Qu.: 1.587   1st Qu.: 1.615   1st Qu.: 1.616  
#>  Median : 2.680   Median : 2.759   Median : 2.694   Median : 2.713  
#>  Mean   : 2.973   Mean   : 3.001   Mean   : 2.970   Mean   : 2.994  
#>  3rd Qu.: 4.019   3rd Qu.: 4.124   3rd Qu.: 4.071   3rd Qu.: 4.136  
#>  Max.   :15.837   Max.   :19.309   Max.   :13.942   Max.   :14.721  
#>       ID9              ID10             ID11             ID12       
#>  Min.   :-1.377   Min.   :-2.583   Min.   :-2.194   Min.   :-1.935  
#>  1st Qu.: 1.631   1st Qu.: 1.745   1st Qu.: 1.667   1st Qu.: 1.693  
#>  Median : 2.774   Median : 2.855   Median : 2.784   Median : 2.806  
#>  Mean   : 3.006   Mean   : 3.079   Mean   : 3.024   Mean   : 3.065  
#>  3rd Qu.: 4.055   3rd Qu.: 4.131   3rd Qu.: 4.036   3rd Qu.: 4.174  
#>  Max.   :16.020   Max.   :13.898   Max.   :13.906   Max.   :14.932  
#>       ID13             ID14             ID15             ID16       
#>  Min.   :-1.788   Min.   :-2.870   Min.   :-2.341   Min.   :-2.691  
#>  1st Qu.: 1.679   1st Qu.: 1.616   1st Qu.: 1.646   1st Qu.: 1.673  
#>  Median : 2.819   Median : 2.727   Median : 2.820   Median : 2.769  
#>  Mean   : 3.063   Mean   : 2.968   Mean   : 3.056   Mean   : 3.037  
#>  3rd Qu.: 4.044   3rd Qu.: 4.045   3rd Qu.: 4.107   3rd Qu.: 4.099  
#>  Max.   :14.722   Max.   :16.234   Max.   :15.453   Max.   :16.481  
#>       ID17             ID18             ID19             ID20       
#>  Min.   :-2.215   Min.   :-2.289   Min.   :-2.208   Min.   :-2.209  
#>  1st Qu.: 1.648   1st Qu.: 1.640   1st Qu.: 1.651   1st Qu.: 1.591  
#>  Median : 2.723   Median : 2.776   Median : 2.780   Median : 2.771  
#>  Mean   : 2.983   Mean   : 3.013   Mean   : 3.054   Mean   : 3.004  
#>  3rd Qu.: 4.029   3rd Qu.: 4.051   3rd Qu.: 4.099   3rd Qu.: 4.074  
#>  Max.   :11.995   Max.   :15.660   Max.   :16.406   Max.   :12.814  
#>       ID21              ID22              ID23             ID24        
#>  Min.   :-8.0482   Min.   :-8.6276   Min.   :-9.369   Min.   :-8.5430  
#>  1st Qu.:-3.3710   1st Qu.:-3.3573   1st Qu.:-3.371   1st Qu.:-3.3773  
#>  Median :-2.5229   Median :-2.4975   Median :-2.543   Median :-2.5046  
#>  Mean   :-2.6163   Mean   :-2.5975   Mean   :-2.659   Mean   :-2.6449  
#>  3rd Qu.:-1.7383   3rd Qu.:-1.7104   3rd Qu.:-1.808   3rd Qu.:-1.7757  
#>  Max.   : 0.8294   Max.   : 0.5153   Max.   : 1.023   Max.   : 0.6307  
#>       ID25             ID26              ID27              ID28       
#>  Min.   :-8.019   Min.   :-9.5781   Min.   :-10.931   Min.   :-9.204  
#>  1st Qu.:-3.356   1st Qu.:-3.3971   1st Qu.: -3.390   1st Qu.:-3.289  
#>  Median :-2.529   Median :-2.5561   Median : -2.491   Median :-2.475  
#>  Mean   :-2.621   Mean   :-2.6450   Mean   : -2.627   Mean   :-2.581  
#>  3rd Qu.:-1.761   3rd Qu.:-1.7877   3rd Qu.: -1.759   3rd Qu.:-1.731  
#>  Max.   : 1.152   Max.   : 0.9988   Max.   :  1.066   Max.   : 1.296  
#>       ID29              ID30        
#>  Min.   :-8.5215   Min.   :-8.1865  
#>  1st Qu.:-3.4093   1st Qu.:-3.3412  
#>  Median :-2.5437   Median :-2.5153  
#>  Mean   :-2.6430   Mean   :-2.6377  
#>  3rd Qu.:-1.7914   3rd Qu.:-1.7777  
#>  Max.   : 0.6812   Max.   : 0.7559  
```
