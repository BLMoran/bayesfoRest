# Build display labels for per-study random effects

Build display labels for per-study random effects

## Usage

``` r
make_study_labels(data)
```

## Arguments

- data:

  Analysis dataset containing the (unique) `Author` key and, ideally,
  `Author_original` and `year`.

## Value

A named character vector: names = Author keys (as they appear in
`r_Author[<key>,Intercept]`); values = "Author (year)".
