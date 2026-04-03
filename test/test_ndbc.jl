#!/usr/bin/env julia
# Regular Julia script version of examples/ndbc.jl (Pluto-free), adapted for tests.
# Run from repo root: julia --project=test test/ndbc_script.jl

import Dates
import Unitful
import HDF5
using Test
import BuoyData: NDBC
using AxisArrays

write_gold = false
gold_path = normpath(joinpath(@__DIR__, "gold", "omnidata_head.hdf5"))
n_gold = 3
tol = 1e-6

# --- Settings ---
buoy = 46050  # PACWAVE
year = 2021
bfile = false

data_dir = normpath(joinpath(@__DIR__, "..", "examples", "data"))

# Test what is available
avail = NDBC.available()
buoys = unique(avail[!, :buoy])

available_swden = NDBC.available(:omnidirectional_spectrum)
available_pacwave = NDBC.available(buoy, :omnidirectional_spectrum; retries = 1) # test passing keyword arguments to HTTP.request
data_swden = NDBC.request(buoy, year, bfile, :omnidirectional_spectrum)

# See metadata for a given buoy
metadata = NDBC.metadata(buoy)
metadata = NDBC.metadata(46053; retries = 1) # test passing keyword arguments to HTTP.request
lats = Unitful.ustrip(metadata["Latitude"])
lons = Unitful.ustrip(metadata["Longitude"])

# Request a variety of data
data_historical = NDBC.request(buoy, year, bfile)
data_thredds = NDBC.request(
    buoy, year, bfile, :omnidirectional_spectrum, :thredds; retries = 1)
d1 = NDBC.read(joinpath(data_dir, "41001i2015.txt.gz"))
d2 = NDBC.read(joinpath(data_dir, "41001i2015.txt"))
d3 = NDBC.read(joinpath(data_dir, "41001i2015.txt.gz"), "swdir2")
d4 = NDBC.read_netcdf(
    joinpath(
        homedir(), ".cache", "JuliaOceanWaves", "BuoyData", "NDBC", "thredds", "46050w2021.nc"), "swden") # tests that data is cached and can be read directly

years = copy(available_pacwave.year)
years = years[(end - 1):end] # trim years for the sake of the test speed

omnidata = NDBC.request(buoy, years[1], bfile, :omnidirectional_spectrum)
data = omnidata.data
time = omnidata.axes[1][:]
for i in 2:length(years)
    global omnidata = NDBC.request(buoy, years[i], bfile, :omnidirectional_spectrum)
    global data = vcat(data, omnidata.data)
    global time = vcat(time, omnidata.axes[1][:])
end
omnidata = AxisArray(data, time)

n_data = max(1, min(n_gold, size(omnidata, 1)))
gold_new = omnidata[1]
if n_data > 1
    for i in 2:n_data
        global gold_new = hcat(gold_new, omnidata[i])
    end
end
gold_new = Unitful.ustrip(gold_new)

if write_gold
    mkpath(dirname(gold_path))
    HDF5.h5open(gold_path, "w") do f
        f["omnidata_head"] = gold_new
    end
    println("Wrote gold file to: ", gold_path)
end

@test isfile(gold_path)
if isfile(gold_path)
    gold_old = HDF5.h5open(gold_path, "r") do f
        read(f["omnidata_head"])
    end

    @test size(gold_old) == size(gold_new)

    for idx in CartesianIndices(gold_old)
        @test isapprox(gold_old[idx], gold_new[idx]; atol = tol, rtol = tol)
    end
end
