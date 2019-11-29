#!/bin/bash
pymol -p <<EOF
load ../xyz/ammonium-nitrate.xyz
set_view (\
    -0.007026420,   -0.089702196,    0.995937169,\
     0.333145678,    0.938859105,    0.086913697,\
    -0.942846119,    0.332405895,    0.023287248,\
     0.000000000,    0.000000000,  -16.265995026,\
     0.000000000,    0.000000000,    1.680233479,\
    12.824234962,   19.707756042,  -20.000000000 )

set stick_radius, 0.06
set sphere_scale, 0.15
show sticks
show spheres
set transparency, 0.9
set opaque_background, off
set ray_opaque_background, off
ray 4000, 4000
png ammonium_nitrate, dpi=1200
EOF

