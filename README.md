# Pointing Recontruction (Model) Maps
PRMaps.jl is a set of function to generate maps of the observed signal by the [LSPE/Strip](https://lspe.roma1.infn.it/strip.html) telescope given an input map of the entire sky.
The scanning strategy of the telescope is aviable [here](https://lspestrip.github.io/Stripeline.jl/latest/scanning/). 

However, the main purpose of PRMaps.jl is the introduction of systematic errors on the configuration angles of the telescope (see [here](https://lspestrip.github.io/Stripeline.jl/latest/prm/) for more details).
In this way it is possible to create maps affected by an intrinsic error and compare them with the ideal ones.

All the maps are represented by [HealpixMaps](https://github.com/ziotom78/Healpix.jl).

## Documentation
A simple documentation is aviable [here](https://teob97.github.io/PRMaps.jl/dev/).

## License
Stripeline is released under the MIT license.
