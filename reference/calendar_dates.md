# Create calendar dates

Create calendar dates by specified time intervals.

## Usage

``` r
calendar_dates(y_start, y_end, by = "1 month")
```

## Arguments

- y_start:

  A number with the first year for the required dates.

- y_end:

  A number with the last year for the required dates.

- by:

  A number or string (as in lubridate package) to be used as interval to
  get the calendar dates. Default to "1 month" which is returning
  monthly dates (first day of the month).

## Value

A vector with dates from the y_start to the y_end years at the time
intervals specified by "by" argument.

## Examples

``` r
if (FALSE) { # \dontrun{
calendar_dates(41, 150, by = "1 month")
calendar_dates(41, 150, by = "1 year")
calendar_dates(-22000, 150, by = "1 year")
} # }
```
