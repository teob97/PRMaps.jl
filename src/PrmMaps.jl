module PrmMaps

import Stripeline as Sl
using Healpix
export makeMap

function makeMap(cam_ang :: Sl.CameraAngles, 
                 telescope_ang :: Sl.TelescopeAngles, 
                 times, τ_s, signal, NSIDE)
    
    pixel_index = Array{Int}(undef, length(times))

    (dirs, Ψ) = Sl.genpointings(cam_ang, times; telescope_ang = telescope_ang) do time_s
        return (0.0, deg2rad(20.0), Sl.timetorotang(time_s, τ_s*60.))
    end
    
    for i in 1:length(times)
        colat, long = dirs[i,:]
        pixel_index[i] = ang2pix(signal, colat, long)
    end
    
    sky_tod = signal.pixels[pixel_index]
    map_values = Sl.tod2map_mpi(pixel_index, sky_tod, 12*(NSIDE^2); unseen = UNSEEN)

    map = HealpixMap{Float64, RingOrder}(NSIDE)
    map.pixels = map_values
    
    map
end

# Importare .fits
# Generare puntamenti
# cilare sui puntamenti e trovare indici associati ai pixel che vedo
# genero TOD
# tod2map
# plottare mappa

end # module PrmMaps
