# MIT License
# 
# Copyright (c) 2022 Matteo Baratto 
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import Stripeline as Sl
import Dates
using Healpix

export Setup
export makeIdealMap, makeErroredMap
export add2pixel!

"""
    Setup(
        sampling_freq_Hz :: Float64,
        total_time_s :: Float64    
    )

Struct containing some useful data.
"""
Base.@kwdef struct Setup
    elevation_ang_rad = 20.0
    sampling_freq_Hz :: Float64 = 0.0
    total_time_s :: Float64 = 0.0
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
    times = 0 : 1.0/setup.sampling_freq_Hz : setup.total_time_s

    for t in times
        
        Sl.genpointings!(wheelfunction, cam_ang, t, dirs, psi)
        pixel_index_ideal = ang2pix(signal, dirs[1], dirs[2])
        
        Sl.genpointings!(wheelfunction, cam_ang, t, dirs, psi; telescope_ang = telescope_ang)

        sky_value = Healpix.interpolate(signal, dirs[1], dirs[2], pixbuf, weightbuf)

        add2pixel!(map, sky_value, pixel_index_ideal, hits)

    end
    return nothing
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
    times = 0 : 1.0/setup.sampling_freq_Hz : setup.total_time_s
    
    for t in times
        
        Sl.genpointings!(wheelfunction, cam_ang, t, dirs, psi)
        
        pixel_index_ideal = ang2pix(signal, dirs[1], dirs[2])

        sky_value = Healpix.interpolate(signal, dirs[1], dirs[2], pixbuf, weightbuf)

        add2pixel!(map, sky_value, pixel_index_ideal, hits)

    end
    return nothing
end

"""
    makeErroredMap(
        cam_ang :: Sl.CameraAngles, 
        telescope_angles :: Sl.TelescopeAngles,
        signal :: Healpix.HealpixMap,
        setup :: Setup
    )

Generate a Healpix map using a telescope model that takes into account
the non idealities to generate the pointing direction. The observed values instead are
calculated taking into account the ideal pointing directions. 
The result is a map affected by an error due to the non idealities of the system.

Return a tuple `(map, hits)::(HealpixMap, HealpixMap)`:
- `map` contains the obserbed values of signal;
- `hits` contains for every pixels the count of how many times the pixel is seen.
"""
function makeErroredMap(
    cam_ang :: Sl.CameraAngles,
    telescope_ang :: Sl.TelescopeAngles,
    signal :: HealpixMap,
    setup :: Setup
    )

    map = HealpixMap{Float64, RingOrder}(signal.resolution.nside)
    hits = HealpixMap{Int32, RingOrder}(signal.resolution.nside)
    wheelfunction = x -> (0.0, deg2rad(20.0), Sl.timetorotang(x, 1.0))

    fillMap!(wheelfunction, map, cam_ang, telescope_ang, signal, setup, hits)

    map.pixels .= map.pixels ./ hits
    return (map, hits)
end

"""
    makeIdealMap(
        cam_ang :: Sl.CameraAngles, 
        signal :: Healpix.HealpixMap,
        setup :: Setup
    )

Generate a Healpix map using an ideal telescope model.

Return a tuple `(map, hits)::(HealpixMap, HealpixMap)`:
- `map` contains the obserbed values of signal;
- `hits` contains for every pixels the count of how many times the pixel is seen.
"""
function makeIdealMap(
    cam_ang :: Sl.CameraAngles,
    signal :: HealpixMap,
    setup :: Setup
    )

    map = HealpixMap{Float64, RingOrder}(signal.resolution.nside)
    hits = HealpixMap{Int32, RingOrder}(signal.resolution.nside)
    wheelfunction = x -> (0.0, deg2rad(20.0), Sl.timetorotang(x, 1.0))

    fillIdealMap!(wheelfunction, map, cam_ang, signal, setup, hits)

    map.pixels = map.pixels ./ hits.pixels
    return (map, hits)
end

# ---------------------------------------------------------------------------------
# -------------- Functions flavour that accept DateTime object --------------------
# ---------------------------------------------------------------------------------

function fillMap!(
    wheelfunction,
    map :: HealpixMap,
    cam_ang :: Sl.CameraAngles,
    telescope_ang :: Sl.TelescopeAngles,
    signal :: HealpixMap,
    setup :: Setup,
    hits :: HealpixMap,
    t_start :: Dates.DateTime
    )
    
    pixbuf = Array{Int}(undef, 4)
    weightbuf = Array{Float64}(undef, 4)
    dirs = Array{Float64}(undef, 1, 2)
    psi = Array{Float64}(undef, 1)
    times = 0 : 1.0/setup.sampling_freq_Hz : setup.total_time_s

    for t in times
        
        Sl.genpointings!(wheelfunction, cam_ang, t, t_start, dirs, psi)
        pixel_index_ideal = ang2pix(signal, dirs[1], dirs[2])
        
        Sl.genpointings!(wheelfunction, cam_ang, t, t_start, dirs, psi; telescope_ang = telescope_ang)

        sky_value = Healpix.interpolate(signal, dirs[1], dirs[2], pixbuf, weightbuf)

        add2pixel!(map, sky_value, pixel_index_ideal, hits)

    end
    return nothing
end

function fillIdealMap!(
    wheelfunction,
    map :: HealpixMap,
    cam_ang :: Sl.CameraAngles,
    signal :: HealpixMap,
    setup :: Setup,
    hits :: HealpixMap,
    t_start :: Dates.DateTime
    )
    
    pixbuf = Array{Int}(undef, 4)
    weightbuf = Array{Float64}(undef, 4)
    dirs = Array{Float64}(undef, 1, 2)
    psi = Array{Float64}(undef, 1)
    times = 0 : 1.0/setup.sampling_freq_Hz : setup.total_time_s
    
    for t in times
        
        Sl.genpointings!(wheelfunction, cam_ang, t, t_start, dirs, psi)
        
        pixel_index_ideal = ang2pix(signal, dirs[1], dirs[2])

        sky_value = Healpix.interpolate(signal, dirs[1], dirs[2], pixbuf, weightbuf)

        add2pixel!(map, sky_value, pixel_index_ideal, hits)

    end
    return nothing
end

"""
    makeErroredMap(
        cam_ang :: Sl.CameraAngles, 
        telescope_angles :: Sl.TelescopeAngles,
        signal :: Healpix.HealpixMap,
        setup :: Setup,
        t_start :: Dates.DateTime
    )

Generate a Healpix map using a telescope model that takes into account
the non idealities to generate the pointing direction. The observed values instead are
calculated using the ideal pointing directions. 
The result is a map affected by an error due to the non idealities of the system.

Return a tuple `(map, hits)::(HealpixMap, HealpixMap)`:
- `map` contains the obserbed values of signal;
- `hits` contains for every pixels the count of how many times the pixel is seen.
"""
function makeErroredMap(
    cam_ang :: Sl.CameraAngles,
    telescope_ang :: Sl.TelescopeAngles,
    signal :: HealpixMap,
    setup :: Setup,
    t_start :: Dates.DateTime
    )

    map = HealpixMap{Float64, RingOrder}(signal.resolution.nside)
    hits = HealpixMap{Int32, RingOrder}(signal.resolution.nside)
    wheelfunction = x -> (0.0, deg2rad(20.0), Sl.timetorotang(x, 1.0))

    fillMap!(wheelfunction, map, cam_ang, telescope_ang, signal, setup, hits, t_start)

    map.pixels .= map.pixels ./ hits
    return (map, hits)
end

"""
    makeIdealMap(
        cam_ang :: Sl.CameraAngles, 
        signal :: Healpix.HealpixMap,
        setup :: Setup
        t_start :: Dates.DateTime
    )

Generate a Healpix map using an ideal telescope model.

Return a tuple `(map, hits)::(HealpixMap, HealpixMap)`:
- `map` contains the obserbed values of signal;
- `hits` contains for every pixels the count of how many times the pixel is seen.
"""
function makeIdealMap(
    cam_ang :: Sl.CameraAngles,
    signal :: HealpixMap,
    setup :: Setup,
    t_start :: Dates.DateTime
    )

    map = HealpixMap{Float64, RingOrder}(signal.resolution.nside)
    hits = HealpixMap{Int32, RingOrder}(signal.resolution.nside)
    wheelfunction = x -> (0.0, deg2rad(20.0), Sl.timetorotang(x, 1.0))

    fillIdealMap!(wheelfunction, map, cam_ang, signal, setup, hits, t_start)

    map.pixels = map.pixels ./ hits.pixels
    return (map, hits)
end