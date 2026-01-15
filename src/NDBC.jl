
module NDBC

using Dates

using HTTP
using TranscodingStreams, CodecZlib
using DelimitedFiles
using DataFrames
using Unitful
using Unitful: Hz, m
using DimensionfulAngles: °ᵃ as °
using AxisArrays
using Spectra

function _available(parameter::AbstractString)
    # scrape website
    url = "https://www.ndbc.noaa.gov/data/historical/" * parameter * "/"
    raw = filter(x -> occursin(".txt.gz", x), split(String(HTTP.get(url).body)))
    # parse
    filenames = map(x -> String(split(x, "\"")[2]), raw)
    buoys = map(x -> x[1:5], filenames)
    years = map(x -> String(split(x, ".")[1][7:end]), filenames)
    # create DataFrame
    data = DataFrame("buoy" => buoys, "year" => years, "b_file" => false)
    # remove entries with bad file names, currently only "42002w2008_old.txt.gz"
    regular_file(y::String) = length(y) == 4
    b_file(y::String) = length(y) == 5 && y[1] == 'b'
    filter!(:year => y -> regular_file(y) || b_file(y), data)
    # b-files
    b_files = filter(:year => b_file, data)
    for row in eachrow(b_files)
        ibuoy = row.buoy
        iyear = row.year[2:end]
        # TODO: there's probably a more efficient way of doing this:
        data.b_file = @. ifelse(data.buoy == ibuoy && data.year == iyear, true, data.b_file)
    end
    filter!(:year => y -> regular_file(y), data)
    data[!, :year] = parse.(Int, data[!, :year])
    sort!(data)
end


function _available(parameter::AbstractString, buoy::Union{AbstractString,Int})
    data = _available(parameter)
    _filterbuoy(data, buoy)
end


function _request(parameter::AbstractString, buoy::Union{AbstractString,Int}, year::Int, b_file::Bool=false)
    # get data
    sep_dict = Dict("swden" => "w", "swdir" => "d", "swdir2" => "i", "swr1" => "j", "swr2" => "k")
    sep = b_file ? sep_dict[parameter] * "b" : sep_dict[parameter]
    filename = string(buoy) * sep * string(year) * ".txt.gz"
    url = "https://www.ndbc.noaa.gov/data/historical/" * parameter * "/" * filename
    raw = transcode(GzipDecompressor, HTTP.get(url).body)
    _read(raw, parameter)
end


function _filterbuoy(data::DataFrame, buoy::Union{AbstractString,Int})
    filter!(row -> row.buoy == string(buoy), data)
    select!(data, Not(:buoy))
end


function _read(file::Union{AbstractString,Vector{UInt8}}, parameter::AbstractString)
    # parse data
    data, header = DelimitedFiles.readdlm(file, header=true)
    header[1] = strip(header[1], '#')
    # datetime
    ncol_date = header[5] == "mm" ? 5 : 4
    datevec = string.(Int.(data[:, 1:ncol_date]), pad=2)
    two_digit_year = length(datevec[1, 1]) == 2
    fmt = two_digit_year ? "yy" : "yyyy"
    fmt *= "mmddHH"
    ncol_date == 5 && (fmt *= "MM")
    dates = DateTime[]
    for row in eachrow(datevec)
        push!(dates, DateTime(join(row), fmt))
    end
    # data
    unit_dict = Dict("swden" => m * m / Hz, "swdir" => °, "swdir2" => °, "swr1" => 1, "swr2" => 1)
    data = data[:, ncol_date+1:end] * unit_dict[parameter]
    # frequency
    frequency = parse.(Float64, header[ncol_date+1:end]) * Hz
    # AxisArray
    AxisArray(data; time=dates, frequency=frequency)
end


"""
    _spectrum(data, nDirections)

Converts the 3D AxisArray (time x frequency x parameter) of NDBC directional
spectrum data into a WaveSpectra.jl Spectrum structure. See the NDBC 
documentation on calculating the directional wave spectrum from the 5 parameters:
https://www.ndbc.noaa.gov/faq/measdes.shtml
"""
function _spectrum(data::AxisArray, nDirections::Int)
    # omnidirectional spectrum S(f)
    omniS = data[parameter=:den]

    # spread function D(f, θ)
    # Note: AxisArrays currently drops axes information when multiplying values 
    # together (e.g. 0.01*r1,  theta - a1, omniS*D)
    r₁ = 0.01 * data[parameter=:r1]
    r₂ = 0.01 * data[parameter=:r2]
    α₁ = data[parameter=:dir]
    α₂ = data[parameter=:dir2]
    θ = Vector(0:360.0/nDirections:360-360/nDirections) * °
    D = 1/π * (0.5 .+ r₁ .* cos.(θ' .- α₁) .+ r₂ .* cos.(2*(θ' .- α₂))) 
    
    # Create WaveSpectra.jl Spectrum structure
    # S(f, θ) = S(f) * D(f, θ)
    return Spectra.Spectrum(omniS .* D, data.axes[1].val, θ)
end


"""
    read(file, parameter=nothing)

Read a local NDBC data file (.txt or .txt.gz). If `parameter` is not provided,
it is inferred from the filename.
"""
function read(file::AbstractString, parameter::Union{AbstractString,Nothing}=nothing)
    param_dict = Dict("w" => "swden", "d"  => "swdir", "i" => "swdir2", "j" => "swr1", "k" => "swr2")
    isnothing(parameter) && (parameter = param_dict[string(basename(file)[6])])
    file[end-2:end] == ".gz" && (file = transcode(GzipDecompressor, Base.read(file, String)))
    _read(file, parameter)
end


"""
    available_omnidirectional()
    available_omnidirectional(buoy)

Return a DataFrame listing available omnidirectional wave spectrum files.

When a `buoy` is provided, returns the available years for that buoy.
"""
function available_omnidirectional()
    _available("swden")
end


"""
    available_omnidirectional(buoy)

See `available_omnidirectional()` for details.
"""
function available_omnidirectional(buoy::Union{AbstractString,Int})
    _available("swden", buoy)
end


"""
    request_omnidirectional(buoy, year, b_file=false)

Download and parse the omnidirectional spectrum for a given buoy and year.
Set `b_file=true` to request the alternate \"b\" file when available.
"""
function request_omnidirectional(buoy::Union{AbstractString,Int}, year::Int, b_file::Bool=false)
    data = _request("swden", buoy, year, b_file)
    time = data.axes[1]
    S = Array{Spectra.OmnidirectionalSpectrum}(undef, length(time))
    for it in 1:length(time)
        S[it] = Spectra.OmnidirectionalSpectrum(data[time[it]].data, data.axes[2].val)
    end
    S = AxisArray(S; time=time)
    return S
end


"""
    available()
    available(buoy)

Return a DataFrame listing buoy-year combinations where all five spectral
files exist (swden, swdir, swdir2, swr1, swr2). When a `buoy` is provided,
returns only rows for that buoy.
"""
function available()
    den = _available("swden")
    dir = _available("swdir")
    dir2 = _available("swdir2")
    r1 = _available("swr1")
    r2 = _available("swr2")

    # buoy-year combinations for which all 5 files exist
    innerjoin(den, dir, dir2, r1, r2, on=[:buoy, :year, :b_file])
end


"""
    available(buoy)

See `available()` for details.
"""
function available(buoy::Union{AbstractString,Int})
    data = available()
    _filterbuoy(data, buoy)
end


"""
    request(buoy, year, b_file=false)

Download and parse the full directional wave spectrum for a buoy and year.
Returns an AxisArray indexed by time, frequency, and parameter.
"""
function request(buoy::Union{AbstractString,Int}, year::Int, b_file::Bool=false)
    buoy = string(buoy)
    den = _request("swden", buoy, year, b_file)
    dir = _request("swdir", buoy, year, b_file)
    dir2 = _request("swdir2", buoy, year, b_file)
    r1 = _request("swr1", buoy, year, b_file)
    r2 = _request("swr2", buoy, year, b_file)

    time, frequency = den.axes
    parameter = AxisArrays.Axis{:parameter}([:den, :dir, :dir2, :r1, :r2])
    data = AxisArray(cat(den, dir, dir2, r1, r2; dims=3), time, frequency, parameter)
    S = Array{Spectra.Spectrum}(undef, length(time))
    for it in 1:length(time)
        S[it] = _spectrum(data[time[it]], 360)
    end
    S = AxisArray(S; time=time)
    return S
end


"""
    metadata(buoy)

Fetch station metadata for a buoy, including latitude, longitude, water depth,
and watch circle radius. Values are returned with Unitful units.
"""
function metadata(buoy::Union{AbstractString,Int})
    keys = ["Water depth", "Watch circle radius"]
    url = "https://www.ndbc.noaa.gov/station_page.php?station=" * string(buoy)
    raw = split(String(HTTP.get(url, status_exception=false).body), '\n')
    dict = Dict{String, Union{Nothing, Quantity}}()
    for key in keys
        data = filter(x -> occursin(key, x), raw)
        if length(data) == 0
            value =  NaN * 1m
        else
            value = replace(replace(data[1], "\t" => ""), r"<.+?(?=>)>" => "")
            value = strip(split(value, ':')[2])
            if key == "Water depth"
                @assert split(value)[2] == "m"
                value = parse(Float64, split(value)[1]) * 1m
            elseif key == "Watch circle radius"
                @assert split(value)[2] == "yards"
                value = parse(Float64, split(value)[1]) * 0.9144m
            end
        end
        dict[key] = value
    end
    # coordinates
    pattern = r"([+-]?([0-9]*[.])?[0-9]+) +([NS]) +([+-]?([0-9]*[.])?[0-9]+) +([EW])"
    loc = filter(x -> occursin(pattern, x), raw)
    if length(loc) == 0
        lat = NaN
        lon = NaN
    else
        loc = match(pattern, loc[1])
        lat = parse(Float64, loc[1])
        (loc[3] == "S") && (lat *= -1)
        lon = parse(Float64, loc[4])
        (loc[6] == "W") && (lon *= -1)
    end
    dict["Latitude"] = lat * 1°
    dict["Longitude"] = lon * 1°
    return dict
end

end
