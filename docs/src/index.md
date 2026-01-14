# NDBC.jl

## Overview

NDBC.jl provides a thin, unit-aware interface to NOAA NDBC historical buoy data.
It downloads the wave spectral files, parses them into time- and frequency-indexed
AxisArrays, and exposes helpers for common metadata queries.

## Quick start

```julia
using NDBC

buoy = 46050

available = NDBC.available()
omnidirectional = NDBC.available_omnidirectional(buoy)

spec = NDBC.request(buoy, 2021, false)
swden = NDBC.request_omnidirectional(buoy, 2021, false)
```

## Detailed usage

### Discover available data

```julia
using NDBC

all_files = NDBC.available()
# columns: buoy, year, b_file

buoy_files = NDBC.available(46050)
```

Some years have an alternate "b" file with corrected data. Set `b_file=true`
when requesting those files.

### Request wave spectra

```julia
using NDBC

data = NDBC.request(46050, 2021, false)
# AxisArray with axes: time, frequency, parameter

using Dates
single_time = data[time=DateTime("2021-01-01T00:40:00")]
```

The returned AxisArray uses Unitful units:

- Spectral density: m^2/Hz
- Mean directions: degrees
- Directional Fourier coefficients: dimensionless

### Read local files

```julia
using NDBC

local = NDBC.read("data/41001i2015.txt.gz")
local_swdir2 = NDBC.read("data/41001i2015.txt.gz", "swdir2")
```

### Station metadata

```julia
using NDBC

meta = NDBC.metadata(46050)
# keys include: "Latitude", "Longitude", "Water depth", "Watch circle radius"
```

## Reference

See the [API reference](@ref API) for all exported functions.

## Theory

NDBC wave files contain spectral density and directional information by frequency.
The `swden` files are omnidirectional spectra S(f). The `swdir` and `swdir2` files
store mean direction and secondary direction estimates by frequency. The `swr1`
and `swr2` files store normalized Fourier coefficients that characterize the
spread of the directional distribution. NDBC.jl preserves this structure and
returns AxisArrays indexed by time and frequency so you can slice and combine
spectra directly in Julia.

!!! note
    Requests download data from NOAA over HTTPS; network access is required.
