using Stripeline
using Healpix
using PRMaps

export makePolMap, makePolMap_old

function makePolMap(
    signal::PolarizedHealpixMap,
    cam_ang::CameraAngles,
    setup::Setup
)
    i_map = makeIdealMap(cam_ang, signal.i, setup)
    q_map = makeIdealMap(cam_ang, signal.q, setup)
    u_map = makeIdealMap(cam_ang, signal.u, setup)

    p_map = HealpixMap{Float64, RingOrder}(setup.NSIDE) 
    
    p_map.pixels = sqrt.(q_map[:].^2 + u_map[:].^2) ./ i_map[:]
    
    return p_map

end

function makePolMap_old(
    signal::PolarizedHealpixMap,
    cam_ang::CameraAngles,
    tel_ang::TelescopeAngles,
    setup::Setup
)
    i_map = makeErroredMap(cam_ang, tel_ang, signal.i, setup)
    q_map = makeErroredMap(cam_ang, tel_ang, signal.q, setup)
    u_map = makeErroredMap(cam_ang, tel_ang, signal.u, setup)

    p_map = HealpixMap{Float64, RingOrder}(setup.NSIDE) 
    
    p_map.pixels = sqrt.(q_map[:].^2 + u_map[:].^2) ./ i_map[:]
    
    return p_map

end







function fillMap!(
    wheelfunction :: Function,
    map :: HealpixMap,
    cam_ang :: Sl.CameraAngles,
    telescope_ang :: Sl.TelescopeAngles,
    signal :: PolarizedHealpixMap,
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

        i_value = Healpix.interpolate(signal.i, dirs[1], dirs[2], pixbuf, weightbuf)
        q_value = Healpix.interpolate(signal.q, dirs[1], dirs[2], pixbuf, weightbuf)
        u_value = Healpix.interpolate(signal.u, dirs[1], dirs[2], pixbuf, weightbuf)

        add2pixel!(map, sqrt(q_value^2+u_value^2)/i_value, pixel_index_ideal, hits)

    end
end

function makePolMap(
    signal :: PolarizedHealpixMap,
    cam_ang :: Sl.CameraAngles,
    telescope_ang :: Sl.TelescopeAngles,
    setup :: Setup
    )

    map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    hits = HealpixMap{Int32, RingOrder}(setup.NSIDE)
    wheelfunction = x -> (0.0, deg2rad(20.0), Sl.timetorotang(x, setup.τ_s*60.))

    fillMap!(wheelfunction, map, cam_ang, telescope_ang, signal, setup, hits)

    map.pixels .= map.pixels ./ hits
    map
end