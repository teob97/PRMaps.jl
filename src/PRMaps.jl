module PRMaps

import Stripeline as Sl
using Healpix, Plots
export Setup
export makeMap, makeErroredMap, makeMapPlots
export getPixelIndex

Base.@kwdef struct Setup
    τ_s :: Float64 = 0.0
    times :: StepRangeLen
    NSIDE :: Int32 = 0
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

function makeErroredMap(
    cam_ang :: Sl.CameraAngles, 
    telescope_ang :: Sl.TelescopeAngles,
    signal :: Healpix.HealpixMap,
    pixel_index_ideal :: Array{Int, 1},
    setup::Setup
    )
    
    pixel_index = getPixelIndex(cam_ang, telescope_ang, signal, setup)
    
    # Return the tod containing the observed values associated with the directions with error
    sky_tod = signal.pixels[pixel_index]

    # Create a map using the values (with error) associated to the pixel that we belive we observe
    map_values = Sl.tod2map_mpi(pixel_index_ideal, sky_tod, 12*(setup.NSIDE^2))

    map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    map.pixels = map_values
    
    map
end

"""
    makeErroredMap(
        cam_ang :: Sl.CameraAngles, 
        telescope_angles :: Vector{Stripeline.TelescopeAngles},
        signal :: Healpix.HealpixMap,
        pixel_index_ideal,
        setup :: Setup
    )

Generate a collection of Healpix maps.
"""
function makeErroredMap(
    cam_ang :: Sl.CameraAngles, 
    telescope_angles :: Vector{Sl.TelescopeAngles},
    signal :: Healpix.HealpixMap,
    pixel_index_ideal :: Array{Int, 1},
    setup :: Setup
    )

    [makeErroredMap(cam_ang, tel, signal, pixel_index_ideal, setup) for tel in telescope_angles]
    
end

# Flavour to make ideal maps
function makeMap(
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

function makeMapPlots(
    cam_ang :: Sl.CameraAngles, 
    telescope_angles :: Sl.TelescopeAngles,
    signal :: Healpix.HealpixMap,
    map_ideal :: Healpix.HealpixMap,
    setup::Setup
    )
    
    pixel_index_ideal = getPixelIndex(cam_ang, nothing, signal, setup)
    map = makeErroredMap(cam_ang, telescope_angles, signal, pixel_index_ideal, setup)
    plot((map-map_ideal)/map_ideal)  
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
end

end # module PrmMaps
