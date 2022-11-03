import Stripeline as Sl
using PRMaps
using Healpix
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
    omegaVAXang_rad = deg2rad(0.6), 
    zVAXang_rad = deg2rad(0.5)
    )

map_ideal = makeMap(cam_ang, nothing, signal, setup)

map_plots = makeMapPlots(cam_ang, telescope_ang, signal, map_ideal, setup)

savefig(map_plots, "examples/demo.png")