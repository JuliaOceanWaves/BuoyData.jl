# BuoyData.jl

[![CI](https://github.com/JuliaOceanWaves/BuoyData.jl/actions/workflows/CI.yaml/badge.svg)](https://github.com/JuliaOceanWaves/BuoyData.jl/actions/workflows/CI.yaml)
[![Coverage](https://codecov.io/gh/JuliaOceanWaves/BuoyData.jl/graph/badge.svg)](https://codecov.io/gh/JuliaOceanWaves/BuoyData.jl)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaoceanwaves.github.io/BuoyData.jl/dev/)

Tools for downloading and working with buoy data in Julia.
Docs: https://juliaoceanwaves.github.io/BuoyData.jl/dev/

## Installation

```julia
import Pkg
Pkg.add(url="https://github.com/JuliaOceanWaves/BuoyData.jl")
```

## Quick start

```julia
import BuoyData: NDBC

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
