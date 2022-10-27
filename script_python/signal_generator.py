import pysm3
import pysm3.units as u
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
# Create a map object specifing resolution and model 
sky = pysm3.Sky(nside=512, preset_strings=["s1"])
# Restituisce un vettore numpy (3,npixel) con i parametri di stokes [I,Q,U] espressi in microKelvin_RJ (Rayleigh-Jeans)
map_40GHz = sky.get_emission(40 * u.GHz)
map_40GHz = map_40GHz.to(u.uK_CMB, equivalencies=u.cmb_equivalencies(40*u.GHz))
# Convert from celestial coordinates to equatorial
map_40GHz_equat = pysm3.apply_smoothing_and_coord_transform(map_40GHz, rot=hp.Rotator(coord=("G", "C")))
# Save into .fits file
hp.fitsfunc.write_map("../input_maps/map_40GHz.fits", map_40GHz_equat, coord = "C", overwrite=True)
# Visualize temperature map
# hp.mollview(map_40GHz_equat[0], min=0, max=1e2, title="I map", unit=map_40GHz_equat.unit)
# Polarization map
# hp.mollview(np.sqrt(map_40GHz_equat[1]**2 + map_40GHz_equat[2]**2), title="P map", min=0, max=1e1, unit=map_40GHz_equat.unit)