# Snowfake
Copyright © 2012 Bart Massey

[This is a work-in-progress. It does run (when given one
argument describing the size of the snowfake to be built),
but I don't think the simulation is quite right yet, and
it's incredibly slow. I'll honestly probably not get to
it until around next Christmas, though. :-)]

This code generates Gravner-Griffeath "Snowfakes", as
described in

> Janko Gravner and David Griffeath, "Modeling Snow Crystal
> Growth II: A mesoscopic lattice map with plausible
> dynamics", Physica D: Nonlinear Phenomena (237)385–404
> 2008, URL http://psoup.math.wisc.edu/papers/h2l.pdf
> accessed 2012/12/17.

Gravner and Griffeath's own C code has been generously made
available at However, I found this code fairly
incomprehensible and hard to adapt to my purposes.

I am also aware of an open-source Python implementation by
Giles Hall and Rachel Holmes, available on
[Github](). However, even though it can use `pypy`, I am
suspicious that I can improve its performance, and also
clean up its output routines.

One of the primary purposes of this code base is to use SVG
rather than a raster format as the underlying output
representation. By rendering SVG hexagons, it is hoped that
one can get a bit better view of what's going on.

This work is available under the "MIT License". Please see
the file `COPYING` in this distribution for license terms.
