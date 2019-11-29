#!/bin/bash
pymol -p <<EOF
load ../xyz/water-ammonium.xyz
set_view (\
    -0.308164150,   -0.868318737,   -0.388657391,\
    -0.060132705,   -0.389943898,    0.918871284,\
    -0.949428082,    0.306534618,    0.067953609,\
     0.000000000,    0.000000000,  -14.178204536,\
    -0.000389874,    0.000036240,    1.186756134,\
    11.178204536,   17.178203583,  -20.000000000 )

set stick_radius, 0.06
set sphere_scale, 0.15
show sticks
show spheres
set transparency, 0.9
set opaque_background, off
set ray_opaque_background, off
ray 4000, 4000
png water_ammonium, dpi=1200
EOF

