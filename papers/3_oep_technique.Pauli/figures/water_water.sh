#!/bin/bash
pymol -p <<EOF
load ../xyz/water-water.xyz
set_view (\
     0.323731661,    0.940398216,   -0.104131773,\
     0.006585651,    0.107815512,    0.994144022,\
     0.946119666,   -0.322521687,    0.028708948,\
     0.000000000,    0.000000000,  -14.178204536,\
    -0.009005010,   -0.003893316,    1.248221874,\
    11.178204536,   17.178203583,  -20.000000000 )

set stick_radius, 0.06
set sphere_scale, 0.15
show sticks
show spheres
set transparency, 0.9
set opaque_background, off
set ray_opaque_background, off
ray 4000, 4000
png water_water, dpi=1200
EOF

