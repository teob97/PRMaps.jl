module PrmMaps

import Stripeline as Sl
using Healpix, Plots
export Setup
export makeMap, makeMapPlots

Base.@kwdef struct Setup
    τ_s = 0.0
    times
    NSIDE :: Int = 0
end

function makeMap(cam_ang :: Sl.CameraAngles, 
                 telescope_ang :: Sl.TelescopeAngles,
                 signal,
                 setup::Setup)
    
    pixel_index = Array{Int}(undef, length(setup.times))
    pixel_index_ideal = Array{Int}(undef, length(setup.times))

    # Generate directions with error
    (dirs, _) = Sl.genpointings(cam_ang, setup.times; telescope_ang = telescope_ang) do time_s
        return (0.0, deg2rad(20.0), Sl.timetorotang(time_s, setup.τ_s*60.))
    end

    # Generate ideal directions
    (dirs_ideal, _) = Sl.genpointings(cam_ang, setup.times; telescope_ang = Sl.TelescopeAngles()) do time_s
        return (0.0, deg2rad(20.0), Sl.timetorotang(time_s, setup.τ_s*60.))
    end
    
    # Find the pixel index associated to the directions with error
    for i in 1:length(setup.times)
        colat, long = dirs[i,:]
        pixel_index[i] = ang2pix(signal, colat, long)
    end

    # Find the pixel index associated to the ideal directions
    for i in 1:length(setup.times)
        colat, long = dirs_ideal[i,:]
        pixel_index_ideal[i] = ang2pix(signal, colat, long)
    end
    
    # Return the tod containing the observed values associated with the directions with error
    sky_tod = signal.pixels[pixel_index]

    # Create a map using the values (with error) associated to the pixel that we belive we observe
    map_values = Sl.tod2map_mpi(pixel_index_ideal, sky_tod, 12*(setup.NSIDE^2))

    map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    map.pixels = map_values
    
    map
end

# Flavour to make ideal maps
function makeMap(cam_ang :: Sl.CameraAngles, 
    telescope_ang :: Nothing,
    signal,
    setup::Setup)

pixel_index = Array{Int}(undef, length(setup.times))

# Generate ideal directions
(dirs, _) = Sl.genpointings(cam_ang, setup.times) do time_s
return (0.0, deg2rad(20.0), Sl.timetorotang(time_s, setup.τ_s*60.))
end

# Find the pixel index associated to the ideal directions
for i in 1:length(setup.times)
colat, long = dirs[i,:]
pixel_index[i] = ang2pix(signal, colat, long)
end

sky_tod = signal.pixels[pixel_index]

map_values = Sl.tod2map_mpi(pixel_index, sky_tod, 12*(setup.NSIDE^2))

map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
map.pixels = map_values

map
end

function makeMapPlots(cam_ang :: Sl.CameraAngles, 
                      telescope_angles,
                      signal,
                      map_ideal,
                      setup::Setup)
    maps = [makeMap(cam_ang, tel, signal, setup) for tel in telescope_angles]
    plots = [plot((map-map_ideal)/map) for map in maps]
    plots
end

end # module PrmMaps
