# scatteringSimulations
A series of matlab scripts for simulating scattering functions.

## Coordinate Conversion Tools

### cart2obl.m

converts cartesian coordiantes to oblate spheroidal coordinates.

usage:

`[zeta, lambda, phi] = cart2obl(x,y,z,c)`

to use this function to convert from 2d cartesian coordinates to elliptical coordinates, use it like so:

`[zeta, lambda, -] = cart2obl(x,0,y,c)`

### obl2cart.m

converts oblate spheroidal coordinates to cartesian coordinates.

usage:

`[x,y,z] = obl2cart(zeta,lambda,phi,c)`

to use this function to convert from elliptical coordiantes to 2d cartesian coordiates, use it like so:

`[x, -, y] = obl2cart(zeta,lambda,0,c)`

## Scattering Functions

These are functions which take a point in cartesian space, a specific wave number, and the size of the shape we are bouncing an acoustic plane wave off of, and returning the pressure of the scattered wave at that point.

### spherical_scatter.m

pressure of scattered wave from a perfectly rigid sphere. Usage:

`[pressure] = spherical_scatter([x,y,z], k, a)`

where `k` is the wave number and `a` is the radius of the sphere the plane wave is hitting. Since the scattered wave ends up being *rotationally symmetrical* about the y axis, it is always more efficient to calculate the scattered pressures on the x-y plane and rotate that to get a 3-d depiction of a scattered wave, rather than computing the scattered pressure for every point in 3D space.

### elliptical_scatter.m

approximation pressure of scattered wave from a perfectly rigid ellipse. This algorithm is something I designed, and thus it's still a bit of a work in progress, but I'll remove this comment once I think it is really accurate. Usage:

`[pressure] = elliptical_scatter([x,y], k, a,c)`

where `k` is the wave number, `a` is the *hyperbolic* radius of the ellipse (the zeta value), and `c` is the distance for each of the ellipse's foci from the origin.

## Visualizations

These are scripts that will take forever to run, but will result in a pretty visualization.

### waveHitsSphere.m

Depiction of a plane wave hitting a perfectly rigid sphere. Aspects of this simulation, like the size of the sphere, the speed of sound, the number of harmonics we simulate, the number of frames we render, are all configurable at the beginning of the script.

### waveHitsEllipse.m

Depiction of a plane wave hitting a perfectly rigid ellipse. Aspects of this simulation, like the size of the ellipse, the speed of sound, the number of harmonics we simulate, the number of frames we render, are all configurable at the beginning of the script.

