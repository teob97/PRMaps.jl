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

using PRMaps
import Stripeline as Sl
using Healpix

export makeErroredMapIQU, makeIdealMapIQU

function fill_IQU_ErroredMap!(
    wheelfunction,
    map :: PolarizedHealpixMap,
    cam_ang :: Sl.CameraAngles,
    telescope_ang :: Sl.TelescopeAngles,
    signal :: PolarizedHealpixMap,
    setup :: Setup,
    hits :: PolarizedHealpixMap,
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

        i_value = Healpix.interpolate(signal.i, dirs[1], dirs[2], pixbuf, weightbuf)
        q_value = Healpix.interpolate(signal.q, dirs[1], dirs[2], pixbuf, weightbuf)
        u_value = Healpix.interpolate(signal.u, dirs[1], dirs[2], pixbuf, weightbuf)

        add2pixel!(map.i, i_value, pixel_index_ideal, hits.i)
        add2pixel!(map.q, q_value, pixel_index_ideal, hits.q)
        add2pixel!(map.u, u_value, pixel_index_ideal, hits.u)

    end
    return nothing
end


"""
    makeErroredMap(
        cam_ang :: Sl.CameraAngles, 
        telescope_angles :: Sl.TelescopeAngles,
        signal :: Healpix.PolarizedHealpixMap,
        setup :: Setup
    )

Generate a PolarizedHealpix map (I,Q,U) using a telescope model that takes into account
the non idealities to generate the pointing direction. The observed values instead are
calculated taking into account the ideal pointing directions. 
The result is a map affected by an error due to the non idealities of the system.

Return a tuple `(map, hits)::(PolarizedHealpixMap, PolarizedHealpixMap)`:
- `map` contains the obserbed values of signal;
- `hits` contains for every pixels the count of how many times the pixel is seen.
"""
function makeErroredMapIQU(
    cam_ang :: Sl.CameraAngles,
    tel_ang :: Sl.TelescopeAngles,
    signal :: PolarizedHealpixMap,
    setup :: PRMaps.Setup
)
    map = PolarizedHealpixMap{Float64, RingOrder}(signal.i.resolution.nside)
    hits = PolarizedHealpixMap{Int32, RingOrder}(signal.i.resolution.nside)
    wheelfunction = x -> (0.0, deg2rad(20.0), Sl.timetorotang(x, 1.0))
    
    fill_IQU_ErroredMap!(wheelfunction, map, cam_ang, tel_ang, signal, setup, hits)

    map.i.pixels = map.i.pixels ./ hits.i.pixels
    map.q.pixels = map.q.pixels ./ hits.q.pixels
    map.u.pixels = map.u.pixels ./ hits.u.pixels

    return (map, hits)
end

function fill_IQU_IdealMap!(
    wheelfunction,
    map :: PolarizedHealpixMap,
    cam_ang :: Sl.CameraAngles,
    signal :: PolarizedHealpixMap,
    setup :: Setup,
    hits :: PolarizedHealpixMap,
    )
    
    pixbuf = Array{Int}(undef, 4)
    weightbuf = Array{Float64}(undef, 4)
    dirs = Array{Float64}(undef, 1, 2)
    psi = Array{Float64}(undef, 1)
    times = 0 : 1.0/setup.sampling_freq_Hz : setup.total_time_s
    
    for t in times
        
        Sl.genpointings!(wheelfunction, cam_ang, t, dirs, psi)
        pixel_index_ideal = ang2pix(hits, dirs[1], dirs[2])
        
        i_value = Healpix.interpolate(signal.i, dirs[1], dirs[2], pixbuf, weightbuf)
        q_value = Healpix.interpolate(signal.q, dirs[1], dirs[2], pixbuf, weightbuf)
        u_value = Healpix.interpolate(signal.u, dirs[1], dirs[2], pixbuf, weightbuf)

        add2pixel!(map.i, i_value, pixel_index_ideal, hits.i)
        add2pixel!(map.q, q_value, pixel_index_ideal, hits.q)
        add2pixel!(map.u, u_value, pixel_index_ideal, hits.u)

    end
    return nothing
end


"""
    makeIdealMap(
        cam_ang :: Sl.CameraAngles, 
        signal :: Healpix.PolarizedHealpixMap,
        setup :: Setup
    )

Generate a PolarizedHealpix map (I,Q,U) using an ideal telescope model.

Return a tuple `(map, hits)::(PolarizedHealpixMap, PolarizedHealpixMap)`:
- `map` contains the obserbed values of signal;
- `hits` contains for every pixels the count of how many times the pixel is seen.
"""
function makeIdealMapIQU(
    cam_ang :: Sl.CameraAngles,
    signal :: PolarizedHealpixMap,
    setup :: PRMaps.Setup
)
    map = PolarizedHealpixMap{Float64, RingOrder}(signal.i.resolution.nside)
    hits = PolarizedHealpixMap{Int64, RingOrder}(signal.i.resolution.nside)
    wheelfunction = x -> (0.0, deg2rad(20.0), Sl.timetorotang(x, 1.0))
    
    fill_IQU_IdealMap!(wheelfunction, map, cam_ang, signal, setup, hits)

    map.i.pixels ./= hits.i.pixels
    map.q.pixels ./= hits.q.pixels
    map.u.pixels ./= hits.u.pixels

    return (map, hits)
end

# ---------------------------------------------------------------------------------
# -------------- Functions flavour that accept DateTime object --------------------
# ---------------------------------------------------------------------------------

function fill_IQU_ErroredMap!(
    wheelfunction,
    map :: PolarizedHealpixMap,
    cam_ang :: Sl.CameraAngles,
    telescope_ang :: Sl.TelescopeAngles,
    signal :: PolarizedHealpixMap,
    setup :: Setup,
    hits :: PolarizedHealpixMap,
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

        i_value = Healpix.interpolate(signal.i, dirs[1], dirs[2], pixbuf, weightbuf)
        q_value = Healpix.interpolate(signal.q, dirs[1], dirs[2], pixbuf, weightbuf)
        u_value = Healpix.interpolate(signal.u, dirs[1], dirs[2], pixbuf, weightbuf)

        add2pixel!(map.i, i_value, pixel_index_ideal, hits.i)
        add2pixel!(map.q, q_value, pixel_index_ideal, hits.q)
        add2pixel!(map.u, u_value, pixel_index_ideal, hits.u)

    end
    return nothing
end


"""
    makeErroredMap(
        cam_ang :: Sl.CameraAngles, 
        telescope_angles :: Sl.TelescopeAngles,
        signal :: Healpix.PolarizedHealpixMap,
        setup :: Setup,
        t_start :: Dates.DateTime
    )

Generate a PolarizedHealpix map (I,Q,U) using a telescope model that takes into account
the non idealities to generate the pointing direction. The observed values instead are
calculated taking into account the ideal pointing directions. 
The result is a map affected by an error due to the non idealities of the system.

Return a tuple `(map, hits)::(PolarizedHealpixMap, PolarizedHealpixMap)`:
- `map` contains the obserbed values of signal;
- `hits` contains for every pixels the count of how many times the pixel is seen.
"""
function makeErroredMapIQU(
    cam_ang :: Sl.CameraAngles,
    tel_ang :: Sl.TelescopeAngles,
    signal :: PolarizedHealpixMap,
    setup :: PRMaps.Setup,
    t_start :: Dates.DateTime
)
    map = PolarizedHealpixMap{Float64, RingOrder}(signal.i.resolution.nside)
    hits = PolarizedHealpixMap{Int32, RingOrder}(signal.i.resolution.nside)
    wheelfunction = x -> (0.0, deg2rad(20.0), Sl.timetorotang(x, 1.0))
    
    fill_IQU_ErroredMap!(wheelfunction, map, cam_ang, tel_ang, signal, setup, hits, t_start)

    map.i.pixels = map.i.pixels ./ hits.i.pixels
    map.q.pixels = map.q.pixels ./ hits.q.pixels
    map.u.pixels = map.u.pixels ./ hits.u.pixels

    return (map, hits)
end

function fill_IQU_IdealMap!(
    wheelfunction,
    map :: PolarizedHealpixMap,
    cam_ang :: Sl.CameraAngles,
    signal :: PolarizedHealpixMap,
    setup :: Setup,
    hits :: PolarizedHealpixMap,
    t_start :: Dates.DateTime
    )
    
    pixbuf = Array{Int}(undef, 4)
    weightbuf = Array{Float64}(undef, 4)
    dirs = Array{Float64}(undef, 1, 2)
    psi = Array{Float64}(undef, 1)
    times = 0 : 1.0/setup.sampling_freq_Hz : setup.total_time_s
    
    for t in times
        
        Sl.genpointings!(wheelfunction, cam_ang, t, t_start, dirs, psi)
        pixel_index_ideal = ang2pix(hits, dirs[1], dirs[2])
        
        i_value = Healpix.interpolate(signal.i, dirs[1], dirs[2], pixbuf, weightbuf)
        q_value = Healpix.interpolate(signal.q, dirs[1], dirs[2], pixbuf, weightbuf)
        u_value = Healpix.interpolate(signal.u, dirs[1], dirs[2], pixbuf, weightbuf)

        add2pixel!(map.i, i_value, pixel_index_ideal, hits.i)
        add2pixel!(map.q, q_value, pixel_index_ideal, hits.q)
        add2pixel!(map.u, u_value, pixel_index_ideal, hits.u)

    end
    return nothing
end


"""
    makeIdealMap(
        cam_ang :: Sl.CameraAngles, 
        signal :: Healpix.PolarizedHealpixMap,
        setup :: Setup,
        t_start :: Dates.DateTime
    )

Generate a PolarizedHealpix map (I,Q,U) using an ideal telescope model.

Return a tuple `(map, hits)::(PolarizedHealpixMap, PolarizedHealpixMap)`:
- `map` contains the obserbed values of signal;
- `hits` contains for every pixels the count of how many times the pixel is seen.
"""
function makeIdealMapIQU(
    cam_ang :: Sl.CameraAngles,
    signal :: PolarizedHealpixMap,
    setup :: PRMaps.Setup,
    t_start :: Dates.DateTime
)
    map = PolarizedHealpixMap{Float64, RingOrder}(signal.i.resolution.nside)
    hits = PolarizedHealpixMap{Int64, RingOrder}(signal.i.resolution.nside)
    wheelfunction = x -> (0.0, deg2rad(20.0), Sl.timetorotang(x, 1.0))
    
    fill_IQU_IdealMap!(wheelfunction, map, cam_ang, signal, setup, hits, t_start)

    map.i.pixels ./= hits.i.pixels
    map.q.pixels ./= hits.q.pixels
    map.u.pixels ./= hits.u.pixels

    return (map, hits)
end