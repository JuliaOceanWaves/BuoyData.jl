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
using WaveSpectra
using NCDatasets

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

function _available(parameter::AbstractString, buoy::Union{AbstractString, Int})
    data = _available(parameter)
    _filterbuoy(data, buoy)
end

function _request(parameter::AbstractString, buoy::Union{AbstractString, Int},
        year::Int, b_file::Bool = false, source::Symbol = :historical)
    # get data
    sep_dict = Dict(
        "swden" => "w", "swdir" => "d", "swdir2" => "i", "swr1" => "j", "swr2" => "k")
    sep = b_file ? sep_dict[parameter] * "b" : sep_dict[parameter]

    if source == :historical
        filename = string(buoy) * sep * string(year) * ".txt.gz"
        url = "https://www.ndbc.noaa.gov/data/historical/" * parameter * "/" * filename
        cache_dir = joinpath(homedir(), ".cache", "JuliaOceanWaves", "BuoyData", "NDBC", "historical")
    elseif source == :thredds
        # THREDDS contains all directional spectra parameters in one .nc file named using the swden seperator, "w"
        filename = string(buoy) * "w" * string(year) * ".nc"
        url = "https://dods.ndbc.noaa.gov/thredds/fileServer/data/" * parameter * "/" * string(buoy) * "/" * filename
        cache_dir = joinpath(homedir(), ".cache", "JuliaOceanWaves", "BuoyData", "NDBC", "thredds")
    else
        throw(ArgumentError("source must be a Symbol with value :historical or :thredds"))
    end

    # If available, use cached data. If not, cache and then use data.
    cache_file = joinpath(cache_dir, filename)
    if !isdir(cache_dir)
        mkpath(cache_dir)
    end
    if !isfile(cache_file) || filesize(cache_file) < 500
        HTTP.download(url, cache_file)
    end

    if source == :historical
        return read(cache_file, parameter)
    elseif source == :thredds
        return read_netcdf(cache_file, parameter)
    end
    
end


function _filterbuoy(data::DataFrame, buoy::Union{AbstractString, Int})
    filter!(row -> row.buoy == string(buoy), data)
    select!(data, Not(:buoy))
end

function _read(file::Union{AbstractString, Vector{UInt8}}, parameter::AbstractString)
    # parse data
    data, header = DelimitedFiles.readdlm(file, header = true)
    header[1] = strip(header[1], '#')
    # datetime
    ncol_date = header[5] == "mm" ? 5 : 4
    datevec = string.(Int.(data[:, 1:ncol_date]), pad = 2)
    two_digit_year = length(datevec[1, 1]) == 2
    fmt = two_digit_year ? "yy" : "yyyy"
    fmt *= "mmddHH"
    ncol_date == 5 && (fmt *= "MM")
    dates = DateTime[]
    for row in eachrow(datevec)
        push!(dates, DateTime(join(row), fmt))
    end
    # data
    unit_dict = Dict(
        "swden" => m * m / Hz, "swdir" => °, "swdir2" => °, "swr1" => 1, "swr2" => 1)
    data = data[:, (ncol_date + 1):end] * unit_dict[parameter]
    # frequency
    frequency = parse.(Float64, header[(ncol_date + 1):end]) * Hz
    # AxisArray
    AxisArray(data; time = dates, frequency = frequency)
end

function _convert_to_spectrum(data::AxisArray, nDirections::Int)
    # See https://www.ndbc.noaa.gov/faq/measdes.shtml on calculating the 
    # directional wave spectrum from the 5 parameters:

    # omnidirectional spectrum S(f)
    omniS = data[parameter = :den]

    # spread function D(f, θ)
    # Note: AxisArrays currently drops axes information when multiplying values 
    # together (e.g. 0.01*r1,  theta - a1, omniS*D)
    r₁ = 0.01 * data[parameter = :r1]
    r₂ = 0.01 * data[parameter = :r2]
    α₁ = data[parameter = :dir]
    α₂ = data[parameter = :dir2]
    θ = Vector(0:(360.0 / nDirections):(360 - 360 / nDirections)) * °
    D = 1 / π * (0.5 .+ r₁ .* cos.(θ' .- α₁) .+ r₂ .* cos.(2 * (θ' .- α₂)))

    # Create WaveSpectra.jl Spectrum structure
    # S(f, θ) = S(f) * D(f, θ)
    return WaveSpectra.Spectrum(omniS .* D, data.axes[1].val, θ)
end

"""
    read(file, parameter=nothing)

Read a local NDBC data file (.txt or .txt.gz). If `parameter` is not provided,
it is inferred from the filename.
"""
function read(file::AbstractString, parameter::Union{AbstractString, Nothing} = nothing)
    param_dict = Dict(
        "w" => "swden", "d" => "swdir", "i" => "swdir2", "j" => "swr1", "k" => "swr2")
    isnothing(parameter) && (parameter = param_dict[string(basename(file)[6])])
    file[(end - 2):end] == ".gz" &&
        (file = transcode(GzipDecompressor, Base.read(file, String)))
    _read(file, parameter)
end

"""
    read_netcdf(file, parameter=nothing)

Read a local NDBC THREDDS data file (.nc). If `parameter` is not provided,
it is inferred from the filename.
"""
function read_netcdf(file::AbstractString, parameter::Union{AbstractString, Nothing} = nothing)
    param_dict = Dict(
        "w" => "swden", "d" => "swdir", "i" => "swdir2", "j" => "swr1", "k" => "swr2")
    isnothing(parameter) && (parameter = param_dict[string(basename(file)[6])])
    
    short_to_long_name_dict = Dict(
        "swden" => "spectral_wave_density", "swdir" => "mean_wave_dir", "swdir2" => "principal_wave_dir", "swr1" => "wave_spectrum_r1", "swr2" => "wave_spectrum_r1"
    )
    nc_parameter = short_to_long_name_dict[parameter]

    unit_dict = Dict(
        "swden" => m * m / Hz, "swdir" => °, "swdir2" => °, "swr1" => 1, "swr2" => 1)
    ds = NCDataset(file, "r")
    data = ds[nc_parameter][1,1,:,:] * unit_dict[parameter]
    dates = ds["time"][:]
    frequency = ds["frequency"][:] * Hz
    close(ds)
    return AxisArray(data'; time = dates, frequency = frequency)
end

"""
    available(:spectrum)
    available(:omnidirectional_spectrum)

Return a DataFrame listing buoy-year combinations where spectral
files exist. For type=:spectrum (default; directional), all five 
parameters are checked (swden, swdir, swdir2, swr1, swr2). 
For type=:omnidirectional_spectrum, only swden is used.
When a `buoy` is provided, returns only rows for that buoy.
"""
function available(type::Symbol = :spectrum)
    if type == :spectrum
        den = _available("swden")
        dir = _available("swdir")
        dir2 = _available("swdir2")
        r1 = _available("swr1")
        r2 = _available("swr2")

        # buoy-year combinations for which all 5 files exist
        return innerjoin(den, dir, dir2, r1, r2, on = [:buoy, :year, :b_file])
    elseif type == :omnidirectional_spectrum
        return _available("swden")
    else
        throw(ArgumentError("type must be a Symbol with value :spectrum or :omnidirectional_spectrum"))
    end
end

"""
    available(buoy, :spectrum)
    available(buoy, :omnidirectional_spectrum)

See `available()` for details.
"""
function available(buoy::Union{AbstractString, Int}, type::Symbol = :spectrum)
    data = available(type)
    _filterbuoy(data, buoy)
end

"""
    request(buoy, year, b_file=false, type=:spectrum)

Download and parse the wave spectrum for a buoy and year.
Returns an AxisArray of WaveSpectra.Spectrum or WaveSpectra.OmnidirectionalSpectrum
structs indexed by time.
type=:spectrum (default) returns the full directional WaveSpectra.Spectrum struct.
type=:omnidirectional_spectrum returns the single direction WaveSpectra.OmnidirectionalSpectrum struct.
source=:historical (default) pulls data from NDBC's historical archive (.txt.gz files)
source=:thredds pulls data from NDBC's THREDDS server, containing both historical and real time data (.nc files)
"""
function request(buoy::Union{AbstractString, Int}, year::Int,
        b_file::Bool = false, type::Symbol = :spectrum, source::Symbol = :historical)
    buoy = string(buoy)
    den = _request("swden", buoy, year, b_file, source)
    if type == :spectrum
        dir = _request("swdir", buoy, year, b_file, source)
        dir2 = _request("swdir2", buoy, year, b_file, source)
        r1 = _request("swr1", buoy, year, b_file, source)
        r2 = _request("swr2", buoy, year, b_file, source)

        time, frequency = den.axes
        parameter = AxisArrays.Axis{:parameter}([:den, :dir, :dir2, :r1, :r2])
        data = AxisArray(cat(den, dir, dir2, r1, r2; dims = 3), time, frequency, parameter)
        S = Array{WaveSpectra.Spectrum}(undef, length(time))
        for it in 1:length(time)
            S[it] = _convert_to_spectrum(data[time[it]], 360)
        end
        return AxisArray(S; time = time)
    elseif type::Symbol == :omnidirectional_spectrum
        time = den.axes[1]
        S = Array{WaveSpectra.OmnidirectionalSpectrum}(undef, length(time))
        for it in 1:length(time)
            S[it] = WaveSpectra.OmnidirectionalSpectrum(
                den[time[it]].data, den.axes[2].val)
        end
        return AxisArray(S; time = time)
    else
        throw(ArgumentError("type must be a Symbol with value :spectrum or :omnidirectional_spectrum"))
    end
end

"""
    metadata(buoy)

Fetch station metadata for a buoy, including latitude, longitude, water depth,
and watch circle radius. Values are returned with Unitful units.
"""
function metadata(buoy::Union{AbstractString, Int})
    keys = ["Water depth", "Watch circle radius"]
    url = "https://www.ndbc.noaa.gov/station_page.php?station=" * string(buoy)
    raw = split(String(HTTP.get(url, status_exception = false).body), '\n')
    dict = Dict{String, Union{Nothing, Quantity}}()
    for key in keys
        data = filter(x -> occursin(key, x), raw)
        if length(data) == 0
            value = NaN * 1m
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
