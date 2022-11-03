module PRMaps

import Stripeline as Sl
using Healpix, Plots
export Setup
export makeMap, makeErroredMap, makeErroredMaps, makeMapPlots
export getPixelIndex

Base.@kwdef struct Setup
    τ_s = 0.0
    times
    NSIDE :: Int = 0
end

function getPixelIndex(
    cam_ang :: Sl.CameraAngles,
    telescope_ang :: Sl.TelescopeAngles,
    signal,
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
    signal,
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

function makeErroredMap(cam_ang :: Sl.CameraAngles, 
                 telescope_ang :: Sl.TelescopeAngles,
                 signal,
                 pixel_index_ideal,
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
    makeErroredMaps(
        cam_ang :: Sl.CameraAngles, 
        telescope_angles,
        signal,
        pixel_index_ideal,
        setup::Setup
        )

Generate a collection of Healpix maps.
"""
function makeErroredMap(cam_ang :: Sl.CameraAngles, 
                  telescope_angles,
                  signal,
                  pixel_index_ideal,
                  setup::Setup)
    maps = [makeErroredMap(cam_ang, tel, signal, pixel_index_ideal, setup) for tel in telescope_angles]
    maps
end

# Flavour to make ideal maps
function makeMap(cam_ang :: Sl.CameraAngles, 
    telescope_ang :: Nothing,
    signal,
    setup::Setup)

    pixel_index = getPixelIndex(cam_ang, telescope_ang, signal, setup)

    sky_tod = signal.pixels[pixel_index]

    map_values = Sl.tod2map_mpi(pixel_index, sky_tod, 12*(setup.NSIDE^2))

    map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    map.pixels = map_values

    map
end

function makeMapPlots(cam_ang :: Sl.CameraAngles, 
                      telescope_angles :: Sl.TelescopeAngles,
                      signal,
                      map_ideal,
                      setup::Setup)
    
    pixel_index_ideal = getPixelIndex(cam_ang, nothing, signal, setup)
    map = makeErroredMap(cam_ang, telescope_angles, signal, pixel_index_ideal, setup)
    plots = plot((map-map_ideal)/map_ideal)
    plots
end

function makeMapPlots(cam_ang :: Sl.CameraAngles, 
                      telescope_angles,
                      signal,
                      map_ideal,
                      setup::Setup)
    
    pixel_index_ideal = getPixelIndex(cam_ang, nothing, signal, setup)
    maps = makeErroredMap(cam_ang, telescope_angles, signal, pixel_index_ideal, setup)
    plots = [plot((map-map_ideal)/map_ideal) for map in maps]
    plots
end

end # module PrmMaps
