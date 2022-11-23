import Stripeline as Sl
using Healpix

export Setup
export makeIdealMap, makeErroredMap
export add2pixel!

Base.@kwdef struct Setup
    τ_s :: Float64 = 0.0
    times :: StepRangeLen
    NSIDE :: Int32 = 0
end

function add2pixel!(map, sky_value, pixel_idx, hits_map)
    map.pixels[pixel_idx] += sky_value
    hits_map.pixels[pixel_idx] += 1
    return nothing
end

function fillMap!(
    wheelfunction,
    map :: HealpixMap,
    cam_ang :: Sl.CameraAngles,
    telescope_ang :: Sl.TelescopeAngles,
    signal :: HealpixMap,
    setup :: Setup,
    hits :: HealpixMap,
    )
    
    pixbuf = Array{Int}(undef, 4)
    weightbuf = Array{Float64}(undef, 4)
    dirs = Array{Float64}(undef, 1, 2)
    psi = Array{Float64}(undef, 1)
    
    for t in setup.times
        
        Sl.genpointings!(wheelfunction, cam_ang, t, dirs, psi)
        pixel_index_ideal = ang2pix(signal, dirs[1], dirs[2])
        
        Sl.genpointings!(wheelfunction, cam_ang, t, dirs, psi; telescope_ang = telescope_ang)

        sky_value = Healpix.interpolate(signal, dirs[1], dirs[2], pixbuf, weightbuf)

        add2pixel!(map, sky_value, pixel_index_ideal, hits)

    end
    return nothing
end

"""
    makeErroredMap(
        cam_ang :: Sl.CameraAngles, 
        telescope_angles,
        signal :: Healpix.HealpixMap,
        setup :: Setup
    )

Generate a Healpix map using a telescope model that takes into account
the non idealities to generate the pointing direction. The observed values instead are
calculated taking into account the ideal pointing directions. 
The result is a map affected by an error due to the non idealities of the system.

"""
function makeErroredMap(
    cam_ang :: Sl.CameraAngles,
    telescope_ang :: Sl.TelescopeAngles,
    signal :: HealpixMap,
    setup :: Setup
    )

    map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    hits = HealpixMap{Int32, RingOrder}(setup.NSIDE)
    wheelfunction = x -> (0.0, deg2rad(20.0), Sl.timetorotang(x, setup.τ_s*60.))

    fillMap!(wheelfunction, map, cam_ang, telescope_ang, signal, setup, hits)

    map.pixels .= map.pixels ./ hits
    return map
end

function fillIdealMap!(
    wheelfunction,
    map :: HealpixMap,
    cam_ang :: Sl.CameraAngles,
    signal :: HealpixMap,
    setup :: Setup,
    hits :: HealpixMap,
    )
    
    pixbuf = Array{Int}(undef, 4)
    weightbuf = Array{Float64}(undef, 4)
    dirs = Array{Float64}(undef, 1, 2)
    psi = Array{Float64}(undef, 1)
    
    for t in setup.times
        
        Sl.genpointings!(wheelfunction, cam_ang, t, dirs, psi)
        
        pixel_index_ideal = ang2pix(signal, dirs[1], dirs[2])

        sky_value = Healpix.interpolate(signal, dirs[1], dirs[2], pixbuf, weightbuf)

        add2pixel!(map, sky_value, pixel_index_ideal, hits)

    end
    return nothing
end

function makeIdealMap(
    cam_ang :: Sl.CameraAngles,
    signal :: HealpixMap,
    setup :: Setup
    )

    map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    hits = HealpixMap{Int32, RingOrder}(setup.NSIDE)
    wheelfunction = x -> (0.0, deg2rad(20.0), Sl.timetorotang(x, setup.τ_s*60.))

    fillIdealMap!(wheelfunction, map, cam_ang, signal, setup, hits)

    map.pixels = map.pixels ./ hits.pixels
    return map
end