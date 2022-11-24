using Stripeline
using Healpix
using PRMaps

export makePolMap, makePolDegreeMap, makePolAngMap
export polAngleMap!, polDegreeMap!

function makePolMap(
    cam_ang::CameraAngles,
    signal::PolarizedHealpixMap,
    setup::Setup
)
    q_map, _ = makeIdealMap(cam_ang, signal.q, setup)
    u_map, _ = makeIdealMap(cam_ang, signal.u, setup)
    
    return (q_map, u_map)

end

function makePolMap(
    cam_ang::CameraAngles,
    tel_ang::TelescopeAngles,
    signal::PolarizedHealpixMap,
    setup::Setup
)
    q_map, _ = makeErroredMap(cam_ang, tel_ang, signal.q, setup)
    u_map, _ = makeErroredMap(cam_ang, tel_ang, signal.u, setup)

    return (q_map, u_map)

end

function polDegreeMap!(
    p_map::HealpixMap,
    i_map::HealpixMap,
    q_map::HealpixMap,
    u_map::HealpixMap
)
    p_map.pixels = sqrt.(q_map[:].^2 + u_map[:].^2) ./ i_map[:]
    return nothing
end

function makePolDegreeMap(
    cam_ang::CameraAngles,
    tel_ang::TelescopeAngles,
    signal::PolarizedHealpixMap,
    setup::Setup
)
    (q_map, u_map) = makePolMap(cam_ang, tel_ang, signal, setup)
    i_map, _ = makeErroredMap(cam_ang, tel_ang, signal, setup)
    p_map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    
    polDegreeMap!(p_map, i_map, q_map, u_map)

    return p_map
end

function makePolDegreeMap(
    cam_ang::CameraAngles,
    signal::PolarizedHealpixMap,
    setup::Setup
)
    (q_map, u_map) = makePolMap(cam_ang, signal, setup)
    i_map, _ = makeIdealMap(cam_ang, signal, setup)
    p_map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    
    polDegreeMap!(p_map, i_map, q_map, u_map)

    return p_map
end

function polAngleMap!(
    ang_map::HealpixMap,
    q_map::HealpixMap,
    u_map::HealpixMap
)
    ang_map.pixels = 0.5 .* atan.(u_map[:], q_map[:]) 
    return nothing
end

function makePolAngMap(
    cam_ang :: Sl.CameraAngles,
    tel_ang :: Sl.TelescopeAngles,
    signal :: PolarizedHealpixMap,
    setup :: Setup
    )

    (q_map, u_map) = makePolMap(cam_ang, tel_ang, signal, setup)
    p_map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    
    polAngleMap!(p_map, q_map, u_map)

    return p_map
end

function makePolAngMap(
    cam_ang :: Sl.CameraAngles,
    signal :: PolarizedHealpixMap,
    setup :: Setup
    )

    (q_map, u_map) = makePolMap(cam_ang, signal, setup)
    p_map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    
    polAngleMap!(p_map, q_map, u_map)

    return p_map
end