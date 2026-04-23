# Quick Start

`BuoyData.NDBC` is the main entry point for NOAA National Data Buoy Center data.
It can discover available records, request historical wave spectra, and read cached
or local files.

```julia
import BuoyData: NDBC

all_available = NDBC.available()
buoy_available = NDBC.available(46050)
spectra = NDBC.request(46050, 2021, false)
meta = NDBC.metadata(46050)
```

The returned spectral data uses `AxisArrays`, `Unitful`, and `WaveSpectra` types
so frequency, direction, and spectral density units remain attached to the data.
