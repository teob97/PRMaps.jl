import Stripeline as Sl
using Healpix
using PrmMaps
using Plots

signal = readMapFromFITS("input_maps/map_40GHz.fits", 1, Float64)

days_of_observation = 2.
fsampl_hz = 59.9
total_time_s = days_of_observation * 24 * 3600

setup = Setup(
    τ_s = 1. / fsampl_hz,
    times = 0 : (1. / fsampl_hz) : total_time_s,
    NSIDE = 512
)

cam_ang = Sl.CameraAngles()

telescope_ang = Sl.TelescopeAngles(
    omegaVAXang_rad = deg2rad(0.5), 
    zVAXang_rad = deg2rad(0.5)
    )

map_ideal = makeMap(cam_ang, nothing, signal, setup)
map_with_error = makeErroredMap(cam_ang, telescope_ang, signal, setup)

final_map = (map_with_error - map_ideal) / map_ideal

savefig(plot(final_map), "demo.png")