module PRMaps

import Stripeline as Sl
using Healpix, Plots
export Setup
export makeIdealMap, makeErroredMap, makeMapPlots, makeErroredMap_old
export getPixelIndex

Base.@kwdef struct Setup
    τ_s :: Float64 = 0.0
    times :: StepRangeLen
    NSIDE :: Int32 = 0
end

function add2map!(map, sky_value, pixel_idx, pixel_hits)
    map.pixels[pixel_idx] += sky_value
    pixel_hits[pixel_idx] += 1
end

#= function makeErroredMap(
    cam_ang,
    telescope_ang,
    signal,
    setup
    )

    map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    hits = zeros(Int32, 1, 12*setup.NSIDE*setup.NSIDE)

    for t in setup.times

        (dirs_ideal, _) = Sl.genpointings(cam_ang, t) do time
            return (0.0, deg2rad(20.0), Sl.timetorotang(time, setup.τ_s*60.))
        end  
        pixel_index_ideal = ang2pix(signal, dirs_ideal[1], dirs_ideal[2])
        
        (dirs, _) = Sl.genpointings(cam_ang, t; telescope_ang = telescope_ang) do time
            return (0.0, deg2rad(20.0), Sl.timetorotang(time, setup.τ_s*60.))
        end  

        sky_value = Healpix.interpolate(signal, dirs[1], dirs[2])

        tod2map!(map, sky_value, pixel_index_ideal, hits)

    end
    
    map
end =#

function makeErroredMap(
    cam_ang :: Sl.CameraAngles,
    telescope_ang :: Sl.TelescopeAngles,
    signal :: HealpixMap,
    setup :: Setup
    )

    map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    hits = zeros(Int32, 12*setup.NSIDE*setup.NSIDE)
    dirs = Array{Float64}(undef, 1, 2)
    psi = Array{Float64}(undef, 1)

    for t in setup.times
        
        wheelfunction = x -> (0.0, deg2rad(20.0), Sl.timetorotang(x, setup.τ_s*60.))

        Sl.genpointings!(wheelfunction, cam_ang, t, dirs, psi)
        pixel_index_ideal = ang2pix(signal, dirs[1], dirs[2])
        
        Sl.genpointings!(wheelfunction, cam_ang, t, dirs, psi; telescope_ang = telescope_ang)
        pixel_index = ang2pix(signal, dirs[1], dirs[2])

        sky_value = signal.pixels[pixel_index]

        add2map!(map, sky_value, pixel_index_ideal, hits)

    end
    
    map.pixels = map.pixels ./ hits
    map
end

"""
    getPixelIndex(
        cam_ang :: Stripeline.CameraAngles,
        telescope_ang :: Stripeline.TelescopeAngles,
        signal :: Healpix.HealpixMap,
        setup :: PRMaps.Setup
    )

This function return the indeces of the pixel seeing by the telescope given:

    - `cam_ang :: Stripeline.CameraAngles` : encoding the boresight directions of the detector;
    - `telescope_ang :: Stripeline.TelescopeAngles` :  encoding the non idealities angles of the telescope;
    - `signal :: Healpix.HealpixMap` : the input map that the telescope is going to observe;
    - `setup :: PRMaps.Setup` : encoding the information about period of observation and resolution.

This function is provided in two flavours. The second one accepting `telescope_ang = nothing` is used to
simulate the ideal case.

"""
function getPixelIndex(
    cam_ang :: Sl.CameraAngles,
    telescope_ang :: Sl.TelescopeAngles,
    signal :: HealpixMap,
    setup :: Setup
    )
    
    pixel_index = Array{Int}(undef, length(setup.times))
    (dirs, _) = Sl.genpointings(cam_ang, setup.times; telescope_ang = telescope_ang) do time_s
        return (0.0, deg2rad(20.0), Sl.timetorotang(time_s, setup.τ_s*60.))
    end
    for i in 1:length(setup.times)
        colat, long = dirs[i,:]
        pixel_index[i] = ang2pix(signal, colat, long)
    end
    pixel_index
end

function getPixelIndex(
    cam_ang :: Sl.CameraAngles,
    telescope_ang :: Nothing,
    signal :: HealpixMap,
    setup :: Setup
    )
    
    pixel_index = Array{Int}(undef, length(setup.times))
    (dirs, _) = Sl.genpointings(cam_ang, setup.times) do time_s
        return (0.0, deg2rad(20.0), Sl.timetorotang(time_s, setup.τ_s*60.))
    end
    for i in 1:length(setup.times)
        colat, long = dirs[i,:]
        pixel_index[i] = ang2pix(signal, colat, long)
    end
    pixel_index
end

"""
    makeErroredMap(
        cam_ang :: Sl.CameraAngles, 
        telescope_angles,
        signal :: Healpix.HealpixMap,
        pixel_index_ideal :: Array{Int, 1},
        setup :: Setup
    )

Generate a collection of Healpix maps using a telescope model that takes into account of
the non idealities to generate the pointing direction. The observed values instead are
calculated taking into account the ideal pointing directions. The result is a map affected
by an error due to the non idealities of the system.

This function comes into two flavours: telescope_ang could be both a single Sl.TelescopeAngles
or a collection of them returning an Array of HealpixMap.

"""
function makeErroredMap_old(
    cam_ang :: Sl.CameraAngles, 
    telescope_ang :: Sl.TelescopeAngles,
    signal :: Healpix.HealpixMap,
    setup::Setup
    )
    
    pixel_index_ideal = getPixelIndex(cam_ang, nothing, signal, setup)
    pixel_index = getPixelIndex(cam_ang, telescope_ang, signal, setup)
    
    # Return the tod containing the observed values associated with the directions with error
    sky_tod = signal.pixels[pixel_index]

    # Create a map using the values (with error) associated to the pixel that we belive we observe
    map_values = Sl.tod2map_mpi(pixel_index_ideal, sky_tod, 12*(setup.NSIDE^2))

    map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    map.pixels = map_values
    
    map
end

#= function makeErroredMap(
    cam_ang :: Sl.CameraAngles, 
    telescope_angles :: Vector{Sl.TelescopeAngles},
    signal :: Healpix.HealpixMap,
    pixel_index_ideal :: Array{Int, 1},
    setup :: Setup
    )

    [makeErroredMap(cam_ang, tel, signal, pixel_index_ideal, setup) for tel in telescope_angles]
    
end =#

# Flavour to make ideal maps
function makeIdealMap(
    cam_ang :: Sl.CameraAngles, 
    telescope_ang :: Nothing,
    signal :: Healpix.HealpixMap,
    setup::Setup
    )

    pixel_index = getPixelIndex(cam_ang, telescope_ang, signal, setup)

    sky_tod = signal.pixels[pixel_index]

    map_values = Sl.tod2map_mpi(pixel_index, sky_tod, 12*(setup.NSIDE^2))

    map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    map.pixels = map_values

    map
end

#= function makeMapPlots(
    cam_ang :: Sl.CameraAngles, 
    telescope_angles :: Sl.TelescopeAngles,
    signal :: Healpix.HealpixMap,
    map_ideal :: Healpix.HealpixMap,
    setup::Setup
    )
    
    pixel_index_ideal = getPixelIndex(cam_ang, nothing, signal, setup)
    map = makeErroredMap(cam_ang, telescope_angles, signal, pixel_index_ideal, setup)
    result = (map-map_ideal)/map_ideal
    plot(result)
end

function makeMapPlots(
    cam_ang :: Sl.CameraAngles, 
    telescope_angles :: Vector{Sl.TelescopeAngles},
    signal :: Healpix.HealpixMap,
    map_ideal :: Healpix.HealpixMap,
    setup::Setup
    )
    
    pixel_index_ideal = getPixelIndex(cam_ang, nothing, signal, setup)
    maps = makeErroredMap(cam_ang, telescope_angles, signal, pixel_index_ideal, setup)
    [plot((map-map_ideal)/map_ideal) for map in maps]
end =#

end # module PrmMaps
