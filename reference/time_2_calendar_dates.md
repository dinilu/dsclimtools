# Modify time to calendar years

Modify a stars object to change the time dimension to calendar years
between specified starting and ending years.

## Usage

``` r
time_2_calendar_dates(data, y_start, y_end, by = "1 month")
```

## Arguments

- data:

  A stars object with the data to be modified.

- y_start:

  A number with the starting year of the data (in calibrated Before
  Present format).

- y_end:

  A number with the ending year of the data (in calibrated Before
  Present format).

- by:

  A number or string specifying the interval of the dates to be used in
  the new stars object.

## Value

A stars object as in data argument but with changed time dimension.

## Examples

``` r
if (FALSE) { # \dontrun{
data <- read_dsclim("data/dsclim/",
                    var = "tasmax",
                    41,
                    150,
                    rcp = "rcp2.6",
                    gcm = "CESM1-CAM5",
                    calendar_dates = FALSE,
                    proxy = FALSE)
time_2_calendar_dates(data,
                      41,
                      150,
                      by = "1 month")
} # }
```
