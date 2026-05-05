# Load dsclim files

Load dsclim (downscaled TraCE21ka) files by years.

## Usage

``` r
read_dsclim(
  folder,
  var,
  y_start,
  y_end,
  rcp = NULL,
  gcm = NULL,
  calendar_dates = FALSE,
  sf = NULL,
  proxy = TRUE
)
```

## Arguments

- folder:

  A character string with the path to the folder where the dsclim
  (downscaled TraCE21ka) files, ordered in folders by variables.

- var:

  A character string with the name of the variables to be loaded.

- y_start:

  A number with the first year to be loaded. Years are provided as cal
  BP (1950 calendar year = 1 cal BP)

- y_end:

  A number with the last year to be loaded.

- rcp:

  When loading future data (1991-2100), the Representative Concentration
  Pathway that want to be loaded. Valid values are: "rcp2.6", "rcp4.5",
  "rcp6.0", or "rcp8.5".

- gcm:

  When loading future data (1991-2100), the General Circulation Models
  that want to be loaded. Valid values are: "CESM1-CAM5",
  "CSIRO-Mk3-6-0", or "IPSL-CM5A-MR".

- calendar_dates:

  A logical value indicating whether dates should be corrected to
  calendar dates in the output.

- sf:

  A sf object to be used for cropping the object, can be points or
  polygons.

- proxy:

  A logical value indicating whether the data should be loaded as a
  stars_proxy object (TRUE) or as a stars object (FALSE).

## Value

A stars object with the dsclim (downscaled TraCE21ka) data, cropped
using the sf object if provided.

## Details

Years should be in the format calibrated Before Present (e.g. 0 for
loading year 1950 in the gregorian calendar).

## Examples

``` r
if (FALSE) { # \dontrun{
read_dsclim("data/dsclim/",
            var = "tasmax",
            -22000,
            -19950,
            proxy = FALSE)
read_dsclim("data/dsclim/", var = "tasmax",
            41,
            150,
            rcp = "rcp2.6",
            gcm = "CESM1-CAM5",
            calendar_dates = TRUE,
            proxy = FALSE)
} # }
```
