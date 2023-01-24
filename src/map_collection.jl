#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# Function to create a collection of maps, this case is useful when you have to simulate 
# differents frequencies. The pointing is always the same but there are differents signals.
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------

using PRMaps
import Stripeline as Sl
using Healpix

export makeErroredMapsIQU, makeIdealMapsIQU


function fill_IQU_IdealMaps!(
    wheelfunction,
    maps,
    cam_ang :: Sl.CameraAngles,
    signals :: Vector{Healpix.PolarizedHealpixMap},
    setup :: Setup,
    hits :: HealpixMap,
    )
    
    pixbuf = Array{Int}(undef, 4)
    weightbuf = Array{Float64}(undef, 4)
    dirs = Array{Float64}(undef, 1, 2)
    psi = Array{Float64}(undef, 1)
    times = 0 : 1.0/setup.sampling_freq_Hz : setup.total_time_s
    
    for t in times
        
        Sl.genpointings!(wheelfunction, cam_ang, t, dirs, psi)
        pixel_index_ideal = ang2pix(hits, dirs[1], dirs[2])

        for (indx,signal) in enumerate(signals)

            i_value = Healpix.interpolate(signal.i, dirs[1], dirs[2], pixbuf, weightbuf)
            q_value = Healpix.interpolate(signal.q, dirs[1], dirs[2], pixbuf, weightbuf)
            u_value = Healpix.interpolate(signal.u, dirs[1], dirs[2], pixbuf, weightbuf)

            maps[indx].i[pixel_index_ideal] += i_value
            maps[indx].q[pixel_index_ideal] += q_value
            maps[indx].u[pixel_index_ideal] += u_value
        end

        hits[pixel_index_ideal] += 1
    end
    return nothing
end

function makeIdealMapsIQU(
    cam_ang :: Sl.CameraAngles,
    signals :: Vector{Healpix.PolarizedHealpixMap},
    setup :: Setup
)
    maps = [ PolarizedHealpixMap{Float64, RingOrder}(signals[1].i.resolution.nside) for i in axes(signals,1)]
    hits = HealpixMap{Int64, RingOrder}(signals[1].i.resolution.nside)
    wheelfunction = x -> (0.0, deg2rad(20.0), Sl.timetorotang(x, 1.0))

    fill_IQU_IdealMaps!(wheelfunction, maps, cam_ang, signals, setup, hits)
    
    for indx in axes(maps,1)
        maps[indx].i.pixels ./= hits.pixels
        maps[indx].q.pixels ./= hits.pixels
        maps[indx].u.pixels ./= hits.pixels
    end

    return (maps, hits)
    
end


function fill_IQU_ErroredMaps!(
    wheelfunction,
    maps,
    cam_ang :: Sl.CameraAngles,
    tel_ang :: Sl.TelescopeAngles,
    signals :: Vector{Healpix.PolarizedHealpixMap},
    setup :: Setup,
    hits :: HealpixMap,
    )
    
    pixbuf = Array{Int}(undef, 4)
    weightbuf = Array{Float64}(undef, 4)
    dirs = Array{Float64}(undef, 1, 2)
    psi = Array{Float64}(undef, 1)
    times = 0 : 1.0/setup.sampling_freq_Hz : setup.total_time_s
    
    for t in times
        
        Sl.genpointings!(wheelfunction, cam_ang, t, dirs, psi)
        pixel_index_ideal = ang2pix(hits, dirs[1], dirs[2])

        Sl.genpointings!(wheelfunction, cam_ang, t, dirs, psi; telescope_ang = tel_ang)
        
        for (indx,signal) in enumerate(signals)

            i_value = Healpix.interpolate(signal.i, dirs[1], dirs[2], pixbuf, weightbuf)
            q_value = Healpix.interpolate(signal.q, dirs[1], dirs[2], pixbuf, weightbuf)
            u_value = Healpix.interpolate(signal.u, dirs[1], dirs[2], pixbuf, weightbuf)

            maps[indx].i[pixel_index_ideal] += i_value
            maps[indx].q[pixel_index_ideal] += q_value
            maps[indx].u[pixel_index_ideal] += u_value
        end

        hits[pixel_index_ideal] += 1

    end
    return nothing
end

function makeErroredMapsIQU(
    cam_ang :: Sl.CameraAngles,
    tel_ang :: Sl.TelescopeAngles,
    signals :: Vector{Healpix.PolarizedHealpixMap},
    setup :: Setup
)
    maps = [ PolarizedHealpixMap{Float64, RingOrder}(signals[1].i.resolution.nside) for i in axes(signals,1)]
    hits = HealpixMap{Int32, RingOrder}(signals[1].i.resolution.nside)
    wheelfunction = x -> (0.0, deg2rad(20.0), Sl.timetorotang(x, 1.0))

    fill_IQU_ErroredMaps!(wheelfunction, maps, cam_ang, tel_ang, signals, setup, hits)

    for m in maps
        m.i.pixels ./= hits.pixels
        m.q.pixels ./= hits.pixels
        m.u.pixels ./= hits.pixels
    end

    return (maps, hits)
    
end