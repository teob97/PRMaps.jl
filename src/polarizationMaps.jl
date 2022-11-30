using Stripeline
using Healpix
using PRMaps

export makePolMap, makePolDegreeMap, makePolAngMap
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

Used internally by [makePolAngMap](@ref).

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
    p_map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    
    polDegreeMap!(p_map, signal_map.i, signal_map.q, signal_map.u)

    return p_map
end

function makePolDegreeMap(
    cam_ang::CameraAngles,
    signal::PolarizedHealpixMap,
    setup::Setup
)
    signal_map, _ = makeIdealMapIQU(cam_ang, signal, setup)
    p_map = HealpixMap{Float64, RingOrder}(setup.NSIDE)

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
map affected by an error. See [makeErroredMap](@ref)

Input:
- `cam_ang :: CameraAngles` encoding the pointing direction of the detector;
- `tel_ang::Stripeline.TelescopeAngles` encoding the non idealities of the telescope;
- `signal::Healpix.PolarizedHealpixMap` the signal (Q,U,I) that the telescope are going to observe;
- `setup::PRMaps.Setup` see [Setup](@ref).
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

Used internally by [makePolAngMap](@ref).

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
    p_map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    
    polAngleMap!(p_map, signal_map.q, signal_map.u)

    return p_map
end

function makePolAngMap(
    cam_ang :: Sl.CameraAngles,
    signal :: PolarizedHealpixMap,
    setup :: Setup
    )

    signal_map, _ = makeIdealMapIQU(cam_ang, signal, setup)
    p_map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    
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
map affected by an error. See [makeErroredMap](@ref)

Input:
- `cam_ang :: CameraAngles` encoding the pointing direction of the detector;
- `tel_ang::Stripeline.TelescopeAngles` encoding the non idealities of the telescope;
- `signal::Healpix.PolarizedHealpixMap` the signal (Q,U,I) that the telescope are going to observe;
- `setup::PRMaps.Setup` see [Setup](@ref).
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