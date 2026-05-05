# Load cropped dsclim files

Load a cropped version of dsclim (downscaled TraCE21ka) .nc files.

## Usage

``` r
.read_dsclim(file, sf = NULL, proxy)
```

## Arguments

- file:

  A string with full name (path and file name) of the .nc file to be
  loaded from the hard disk drive.

- sf:

  A sf object from the sf package to be used as template or mask to crop
  the downscaled trace .nc file.

- proxy:

  A logical value indicating whether the data should be loaded as a
  stars_proxy object (TRUE) or as a stars object (FALSE).

## Value

A stars object with the downscaled trace data.
