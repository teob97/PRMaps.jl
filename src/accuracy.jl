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
import AngleBetweenVectors: angle

export accuracy

"""
    accuracy(
        cam_ang :: Sl.CameraAngles,
        tel_ang :: Sl.TelescopeAngles,
        setup :: PRMaps.Setup,
        t_start :: Dates.DateTime
    )

Compute the mean difference between the ideal poining direction and the errored direction.
The rerun value is in ARC_MIN.
"""
function accuracy(
    cam_ang :: Sl.CameraAngles,
    tel_ang :: Sl.TelescopeAngles,
    setup :: PRMaps.Setup,
    t_start :: Dates.DateTime
)
    
    wheelfunction = x -> (0.0, deg2rad(20.0), Sl.timetorotang(x, 1.0))

    dirs_ideal = Array{Float64}(undef, 1, 2)
    dirs_errored = Array{Float64}(undef, 1, 2)

    psi = Array{Float64}(undef, 1)
    
    times = 0 : 1.0/setup.sampling_freq_Hz : setup.total_time_s
    
    count = 0
    mean = 0

    for t in times
        
        Sl.genpointings!(wheelfunction, cam_ang, t, t_start, dirs_ideal, psi)
        Sl.genpointings!(wheelfunction, cam_ang, t, t_start, dirs_errored, psi; telescope_ang = tel_ang)

        dirs_ideal_vec = ang2vec(dirs_ideal[1], dirs_ideal[2])
        dirs_errored_vec = ang2vec(dirs_errored[1], dirs_errored[2])

        angle = AngleBetweenVectors.angle(dirs_ideal_vec, dirs_errored_vec)
        count += 1

        mean = mean + (angle - mean)/count

    end

    return mean * (60 * 180) / π

end