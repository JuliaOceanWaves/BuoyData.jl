# NDBC.jl

[![CI](https://github.com/JuliaOceanWaves/NDBC.jl/actions/workflows/CI.yaml/badge.svg)](https://github.com/JuliaOceanWaves/NDBC.jl/actions/workflows/CI.yaml)
[![Coverage](https://codecov.io/gh/JuliaOceanWaves/NDBC.jl/graph/badge.svg)](https://codecov.io/gh/JuliaOceanWaves/NDBC.jl)

Tools for downloading and working with NOAA NDBC buoy data in Julia.

## Installation

```julia
import Pkg
Pkg.add(url="https://github.com/JuliaOceanWaves/NDBC.jl")
```

## Quick start

```julia
using NDBC

buoy = 46050

available = NDBC.available()
omnidirectional = NDBC.available_omnidirectional(buoy)

spec = NDBC.request(buoy, 2021, false)
swden = NDBC.request_omnidirectional(buoy, 2021, false)
```

## Examples

See the scripts under `examples/` and `test/` for end-to-end usage, including plotting.

## Tests

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

Some tests download data and may require network access.
