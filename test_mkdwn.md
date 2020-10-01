test markdown
================
luke durell
10/1/2020

Hello this is from Rmarkdown.

``` r
#- here's some code that does something
r <- 5
r + 5
```

    ## [1] 10

``` r
for (i in 4:(r+4)) {
  message(paste0("This is the ", i+4, "th iteration!"))
}
```

    ## This is the 8th iteration!

    ## This is the 9th iteration!

    ## This is the 10th iteration!

    ## This is the 11th iteration!

    ## This is the 12th iteration!

    ## This is the 13th iteration!
