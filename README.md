# libzinc

Zinc is a C++ library for spatial processing.
 - Zinc provides efficient functions for processing Morton codes, and things built from Morton codes.
 - Intervals are a start and end point on the Morton curve.
 - Regions are a finite list of intervals, able to represent arbitrary regions in N-dimensional space.
   - Morton regions a.k.a. linear octree
 - Intervals and Regions can store data, so they can be used as map types, not just set types.
 - AABBs (Axis-Aligned Bounding Boxes) can be used for creating regions
 - This is alpha software, but it may be useful

[Read our blog](LINK HERE) for more detail

## TODO

 - The Morton code type is currently limited to a single word size, so the Dimension and BitsPerDimension are currently limited to 2,32 and 3,21.
   - This is under active [development](https://github.com/paddygord/bitarray)

## Usage

For usage examples, see `test/zinc-test.cc`

## Dependencies

 - A C++17 compiler, we use clang++ by default for better sanitisation

For testing and installing:
 - [Meson](https://mesonbuild.com/) `sudo apt install meson`

## Building

`./build.sh`

## Testing

`ninja test -C out`

## Installing

`ninja install -C out` (needs permissions for `/usr/local/include`)

## Name

Zn comes from Z-order curve in N dimensions, and all the cool kids name their libraries after elements.

## See also

 - [libmorton](https://github.com/Forceflow/libmorton)
 - [google s2](https://github.com/google/s2geometry)
