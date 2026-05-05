# Read TraCE21ka datasets

Read original TraCE21ka data, returning a stars object. This function
differs from the ones in the dsclim package, which return the data in
the format of the climate4r framework.

## Usage

``` r
read_trace(folder, var, sf = NULL)
```

## Arguments

- folder:

  character string indicating the path to the folder with the TraCE21ka
  dataset.

- var:

  character string indicating the name of the variable that is to be
  readed.

- sf:

  a sf vector object for which the TraCE21ka data should be limited
  (cropped).

## Value

The function returns an stars object, usually in the form of raster
data, except if the argument sf is provided, in which case the object is
a vector stars object.

## Examples

``` r
if (FALSE) { # \dontrun{
read_trace("data/TraCE21ka/", "TSMX")
} # }
```
