#!/usr/bin/env julia
# Regular Julia script version of examples/ndbc.jl (Pluto-free), adapted for tests.
# Run from repo root: julia --project=test test/ndbc_script.jl

import Dates
import Unitful
import HDF5
using Test

include(joinpath(@__DIR__, "..", "src", "NDBC.jl"))

write_gold = false
enable_plotting = false
gold_path = normpath(joinpath(@__DIR__, "gold", "omnidata_head.hdf5"))
n_gold = 3
tol = 1e-6

# --- Settings ---
buoy = 46050  # PACWAVE
year = 2021
bfile = false

data_dir = normpath(joinpath(@__DIR__, "..", "examples", "data"))

swden = NDBC.available_omnidirectional()
pw_swden = NDBC.available_omnidirectional(buoy)
data_swden = NDBC.request_omnidirectional(buoy, year, bfile)

avail = NDBC.available()
avail_pw = NDBC.available(buoy)
data = NDBC.request(buoy, year, bfile)

d1 = NDBC.read(joinpath(data_dir, "41001i2015.txt.gz"))
d2 = NDBC.read(joinpath(data_dir, "41001i2015.txt"))
d3 = NDBC.read(joinpath(data_dir, "41001i2015.txt.gz"), "swdir2")

metadata = NDBC.metadata(buoy)
buoys = unique(NDBC.available()[!, :buoy])

lats = fill(NaN, length(buoys))
lons = fill(NaN, length(buoys))
for i in eachindex(buoys)
    global data = NDBC.metadata(buoys[i])
    lats[i] = Unitful.ustrip(data["Latitude"])
    lons[i] = Unitful.ustrip(data["Longitude"])
end
mask = .!isnan.(lats)
ulats = lats[mask]
ulons = lons[mask]

if enable_plotting
    import WGLMakie
    import GeoMakie
    import Plots

    fig = WGLMakie.Figure(resolution = (1300, 1300))
    ga = GeoMakie.GeoAxis(
        fig[1, 1];
        dest = "+proj=moll"
    )
    WGLMakie.image!(ga, -180 .. 180, -90 .. 90, rotr90(GeoMakie.earth());
        interpolate = false, inspectable = false)
    WGLMakie.scatter!(ga, ulons, ulats, color = "tomato2", inspectable = true,
        inspector_label = (f, i, p) -> buoys[i])
    WGLMakie.DataInspector(fig)
    display(fig)

    S = data_swden[time = Dates.DateTime("2021-01-01T00:40:00")]
    display(Plots.plot(S.axes[1].val, S))
else
    println("Skipping plotting (enable_plotting = false).")
end

years = copy(avail_pw.year)
deleteat!(years, findall(x -> (x == 2014 || x == 2015), years)) # data issues in 2014/2015
years = years[(end - 1):end] # trim years for the sake of the test speed

alldata = NDBC.request(buoy, years[1], bfile)
for i in 2:length(years)
    global alldata = vcat(alldata, NDBC.request(buoy, years[i], bfile))
end

omnidata = NDBC.request_omnidirectional(buoy, years[1], bfile)
for i in 2:length(years)
    global omnidata = vcat(omnidata, NDBC.request_omnidirectional(buoy, years[i], bfile))
end

if enable_plotting
    display(Plots.plot(Unitful.ustrip.(omnidata.data)))
end

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
