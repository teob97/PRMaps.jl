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

using Stripeline
using Healpix
using PRMaps

export makePolDegreeMap, makePolAngMap
export polAngleMap!, polDegreeMap!
export differenceAngMaps

"""
    function polDegreeMap!(
        p_map::HealpixMap,
        i_map::HealpixMap,
        q_map::HealpixMap,
        u_map::HealpixMap
    )

This function calculate the polarization degree map and save it
in `ang_map`. The degree of polarization is defined as:

sqrt(Q^2 + U^2) / I

Used internally by [`makePolAngMap`](@ref).

Input:
- `p_map::HealpixMap` an empty HealpixMap
- `i_map::HealpixMap` I parameter HealpixMap 
- `q_map::HealpixMap` Q parameter HealpixMap 
- `u_map::HealpixMap` U parameter HealpixMap
Output:
- nothing
"""
function polDegreeMap!(
    p_map::HealpixMap,
    i_map::HealpixMap,
    q_map::HealpixMap,
    u_map::HealpixMap
)
    @assert Healpix.conformables(p_map,i_map)
    @assert Healpix.conformables(p_map,q_map)
    @assert Healpix.conformables(p_map,u_map)

    p_map.pixels = sqrt.(q_map[:].^2 + u_map[:].^2) ./ i_map[:]
    return nothing
end

function makePolDegreeMap(
    cam_ang::CameraAngles,
    tel_ang::TelescopeAngles,
    signal::PolarizedHealpixMap,
    setup::Setup
)
    signal_map, _ = makeErroredMapIQU(cam_ang, tel_ang, signal, setup)
    p_map = HealpixMap{Float64, RingOrder}(signal.i.resolution.nside)
    
    polDegreeMap!(p_map, signal_map.i, signal_map.q, signal_map.u)

    return p_map
end

function makePolDegreeMap(
    cam_ang::CameraAngles,
    signal::PolarizedHealpixMap,
    setup::Setup
)
    signal_map, _ = makeIdealMapIQU(cam_ang, signal, setup)
    p_map = HealpixMap{Float64, RingOrder}(signal.i.resolution.nside)

    polDegreeMap!(p_map, signal_map.i, signal_map.q, signal_map.u)

    return p_map
end

"""
    makePolDegreeMap(
        cam_ang::Stripeline.CameraAngles,
        signal::Healpix.PolarizedHealpixMap,
        setup::PRMaps.Setup
    )

    makePolDegreeMap(
        cam_ang::Stripeline.CameraAngles,
        tel_ang::Stripeline.TelescopeAngles,
        signal::Healpix.PolarizedHealpixMap,
        setup::PRMaps.Setup
    )

Return the polarization degree map of the sky observed by a telescope 
whose camera point towards a direction encoded by
the CameraAngles struct.

The second flavour, with a TelescopeAngles as input, produce a
map affected by an error. See [`makeErroredMap`](@ref)

Input:
- `cam_ang :: CameraAngles` encoding the pointing direction of the detector;
- `tel_ang::Stripeline.TelescopeAngles` encoding the non idealities of the telescope;
- `signal::Healpix.PolarizedHealpixMap` the signal (Q,U,I) that the telescope are going to observe;
- `setup::PRMaps.Setup` see [`Setup`](@ref).
Output:
- `p_map` an HealpixMap containing observed degree of polarization. 
"""
makePolDegreeMap

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

"""
    function polAngleMap!(
        ang_map::HealpixMap,
        q_map::HealpixMap,
        u_map::HealpixMap
    )

This function calculate the polarization angle map and save it
in `ang_map`. The angle of polarization is defined as:

0.5 * atan(U, Q)

Used internally by [`makePolAngMap`](@ref).

Input:
- `ang_map::HealpixMap` an empty HealpixMap
- `q_map::HealpixMap` Q parameter HealpixMap 
- `u_map::HealpixMap` U parameter HealpixMap
Output:
- nothing
"""
function polAngleMap!(
    ang_map::HealpixMap,
    q_map::HealpixMap,
    u_map::HealpixMap
)
    @assert Healpix.conformables(ang_map,q_map)
    @assert Healpix.conformables(ang_map,u_map)

    ang_map.pixels = 0.5 .* atan.(u_map[:], q_map[:]) 
    return nothing
end

function makePolAngMap(
    cam_ang :: Sl.CameraAngles,
    tel_ang :: Sl.TelescopeAngles,
    signal :: PolarizedHealpixMap,
    setup :: Setup
    )

    signal_map, _ = makeErroredMapIQU(cam_ang, tel_ang, signal, setup)
    p_map = HealpixMap{Float64, RingOrder}(signal.i.resolution.nside)
    
    polAngleMap!(p_map, signal_map.q, signal_map.u)

    return p_map
end

function makePolAngMap(
    cam_ang :: Sl.CameraAngles,
    signal :: PolarizedHealpixMap,
    setup :: Setup
    )

    signal_map, _ = makeIdealMapIQU(cam_ang, signal, setup)
    p_map = HealpixMap{Float64, RingOrder}(signal.i.resolution.nside)
    
    polAngleMap!(p_map, signal_map.q, signal_map.u)

    return p_map
end

"""
    makePolAngMap(
        cam_ang::Stripeline.CameraAngles,
        signal::Healpix.PolarizedHealpixMap,
        setup::PRMaps.Setup
    )

    makePolAngMap(
        cam_ang::Stripeline.CameraAngles,
        tel_ang::Stripeline.TelescopeAngles,
        signal::Healpix.PolarizedHealpixMap,
        setup::PRMaps.Setup
    )

Return the polarization angle map of the sky observed by a telescope 
whose camera point towards a direction encoded by
the CameraAngles struct.

The second flavour, with a TelescopeAngles as input, produce a
map affected by an error. See [`makeErroredMap`](@ref)

Input:
- `cam_ang :: CameraAngles` encoding the pointing direction of the detector;
- `tel_ang::Stripeline.TelescopeAngles` encoding the non idealities of the telescope;
- `signal::Healpix.PolarizedHealpixMap` the signal (Q,U,I) that the telescope are going to observe;
- `setup::PRMaps.Setup` see [`Setup`](@ref).
Output:
- `p_map` an HealpixMap containing observed angle of polarization. 
"""
makePolAngMap

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

"""
    differenceAngMaps(
        a :: HealpixMap,
        b :: HealpixMap
    )

Calculate the difference between two maps containing the observed
polarization angle.

Return an HealpixMap.

Angles MUST be expressed in RADIANS.
"""
function differenceAngMaps(
    a :: HealpixMap,
    b :: HealpixMap
)
    @assert Healpix.conformables(a,b) ["Input maps haven't same shape and ordering."]
    c = HealpixMap{Float64, RingOrder}(a.resolution.nside)
    r = rem2pi.(a.pixels.-b.pixels, RoundNearest)
    c.pixels = mod.(r[:].+π/2, π) .- π/2
    return c
end