# BuoyData.jl

## Overview

BuoyData.jl provides a thin, unit-aware interface to historical buoy data.
The NDBC submodule provides access to NOAA National Data Buoy Center historical data.
It queries available data, ingests metadata, downloads the wave spectral files, 
and parses them into time-indexed AxisArrays of WaveSpectra.jl datatypes.

## Quick start

```julia
import BuoyData

buoy = 46050

available = BuoyData.NDBC.available()
omnidirectional = BuoyData.NDBC.available(buoy, :omnidirectional_spectrum)

spec = BuoyData.NDBC.request(buoy, 2021, false, :thredds)
swden = BuoyData.NDBC.request(buoy, 2021, false, :omnidirectional_spectrum)
```

## Detailed usage

### Discover available data

```julia
import BuoyData: NDBC

all_files = NDBC.available()
# columns: buoy, year, b_file

buoy_files = NDBC.available(46050)
```

Some years have an alternate "b" file with corrected data. Set `b_file=true`
when requesting those files.

### Request wave spectra

```julia
import BuoyData: NDBC

data = NDBC.request(46050, 2021, false)
# AxisArray with axes: time. 
# WaveSpectra.Spectra object with axes frequency and direction

using Dates
single_time = data[time=DateTime("2021-01-01T00:40:00")]
```

The returned AxisArray contains Unitful WaveSpectra datatypes.

### Read local files

```julia
import BuoyData: NDBC

local_swdir2 = NDBC.read("examples/data/41001i2015.txt.gz")
local_swdir2 = NDBC.read("examples/data/41001i2015.txt.gz", "swdir2")
```

Requested data (except metadata and available data) is automatically cached locally in
"~/.cache/JuliaOceanWaves/BuoyData/" to prevent excessive and expensive API calls. 
Data is automatically saved and read when calling the request function, though files could
be read maunally as well using ``NDBC.read``:

```julia
import BuoyData: NDBC

cached_file = "~/.cache/JuliaOceanWaves/BuoyData/NDBC/historical/46050w2021.txt.gz"
46050_swden = NDBC.read(cached_file)
```

### Station metadata

```julia
import BuoyData: NDBC

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
spread of the directional distribution. 
All spectral and directional parameters are stored in a single NetCDF file when 
downloading data from NOAA's THREDDS server.
NDBC.jl preserves this structure and
returns AxisArrays of WaveSpectra Spectra or OmnidirectionalSpectra structures
indexed by time.

As recorded in the [NDBC documentation](https://www.ndbc.noaa.gov/faq/measdes.shtml),
the directional spectrum, $S(f,\theta)$, is defined by the omnidirectional spectrum, 
$S_o(f)$, and the directional spreading function, $D(f,\theta)$:

$$S(f,\theta) = S_o(f) * D(f,\theta)$$

where

$$D(f,\theta) = \frac{1}{\pi}\left(0.5 + R_1 cos(\theta - \alpha_1) + R_2 cos(\theta - \alpha_2)\right)$$

and $R_1, R_2$ are the normalized (nondimensional) polar coordinates of the Fourier 
coefficients ("swr1" and "swr2") respectively. $\alpha_1, \alpha_2$ are mean and principal 
wave directions ("swdir1" and "swdir2") respectively.

!!! note
    Requests download data from NOAA over HTTPS; network access is required.
    Data is automatically cached to decrease computational and network expenses. 
