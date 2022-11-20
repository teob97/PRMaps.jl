using Stripeline
using Healpix
using PRMaps

export makePolDegreeMap, makePolAngMap,makePolMap_old

function makePolDegreeMap(
    signal :: PolarizedHealpixMap,
    cam_ang :: Sl.CameraAngles,
    telescope_ang :: Sl.TelescopeAngles,
    setup :: Setup
    )

    map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    hits = HealpixMap{Int32, RingOrder}(setup.NSIDE)
    wheelfunction = x -> (0.0, deg2rad(20.0), Sl.timetorotang(x, setup.τ_s*60.))

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

    map.pixels .= map.pixels ./ hits
    map
end


# Make the ideal polarization map
function makePolDegreeMap(
    signal :: PolarizedHealpixMap,
    cam_ang :: Sl.CameraAngles,
    setup :: Setup
    )

    map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    hits = HealpixMap{Int32, RingOrder}(setup.NSIDE)
    wheelfunction = x -> (0.0, deg2rad(20.0), Sl.timetorotang(x, setup.τ_s*60.))

    pixbuf = Array{Int}(undef, 4)
    weightbuf = Array{Float64}(undef, 4)
    dirs = Array{Float64}(undef, 1, 2)
    psi = Array{Float64}(undef, 1)
    
    for t in setup.times
        
        Sl.genpointings!(wheelfunction, cam_ang, t, dirs, psi)
        pixel_index_ideal = ang2pix(signal, dirs[1], dirs[2])

        i_value = Healpix.interpolate(signal.i, dirs[1], dirs[2], pixbuf, weightbuf)
        q_value = Healpix.interpolate(signal.q, dirs[1], dirs[2], pixbuf, weightbuf)
        u_value = Healpix.interpolate(signal.u, dirs[1], dirs[2], pixbuf, weightbuf)

        add2pixel!(map, sqrt(q_value^2+u_value^2)/i_value, pixel_index_ideal, hits)

    end

    map.pixels = map.pixels ./ hits.pixels
    map
end

function makePolAngMap(
    signal :: PolarizedHealpixMap,
    cam_ang :: Sl.CameraAngles,
    telescope_ang :: Sl.TelescopeAngles,
    setup :: Setup
    )

    map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    hits = HealpixMap{Int32, RingOrder}(setup.NSIDE)
    wheelfunction = x -> (0.0, deg2rad(20.0), Sl.timetorotang(x, setup.τ_s*60.))

    pixbuf = Array{Int}(undef, 4)
    weightbuf = Array{Float64}(undef, 4)
    dirs = Array{Float64}(undef, 1, 2)
    psi = Array{Float64}(undef, 1)
    
    for t in setup.times
        
        Sl.genpointings!(wheelfunction, cam_ang, t, dirs, psi)
        pixel_index_ideal = ang2pix(signal, dirs[1], dirs[2])
        
        Sl.genpointings!(wheelfunction, cam_ang, t, dirs, psi; telescope_ang = telescope_ang)

        q_value = Healpix.interpolate(signal.q, dirs[1], dirs[2], pixbuf, weightbuf)
        u_value = Healpix.interpolate(signal.u, dirs[1], dirs[2], pixbuf, weightbuf)

        add2pixel!(map, 0.5 * atan(u_value , q_value), pixel_index_ideal, hits)

    end

    map.pixels .= map.pixels ./ hits
    map
end


# Make the ideal polarization angle map
function makePolAngMap(
    signal :: PolarizedHealpixMap,
    cam_ang :: Sl.CameraAngles,
    setup :: Setup
    )

    map = HealpixMap{Float64, RingOrder}(setup.NSIDE)
    hits = HealpixMap{Int32, RingOrder}(setup.NSIDE)
    wheelfunction = x -> (0.0, deg2rad(20.0), Sl.timetorotang(x, setup.τ_s*60.))

    pixbuf = Array{Int}(undef, 4)
    weightbuf = Array{Float64}(undef, 4)
    dirs = Array{Float64}(undef, 1, 2)
    psi = Array{Float64}(undef, 1)
    
    for t in setup.times
        
        Sl.genpointings!(wheelfunction, cam_ang, t, dirs, psi)
        pixel_index_ideal = ang2pix(signal, dirs[1], dirs[2])

        q_value = Healpix.interpolate(signal.q, dirs[1], dirs[2], pixbuf, weightbuf)
        u_value = Healpix.interpolate(signal.u, dirs[1], dirs[2], pixbuf, weightbuf)

        add2pixel!(map, 0.5 * atan(u_value , q_value), pixel_index_ideal, hits)

    end

    map.pixels = map.pixels ./ hits.pixels
    map
end


# --------------
# Slow funciotns
# --------------

function makePolMap_old(
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