#!/bin/bash
pymol -p <<EOF
load ../xyz/water-methanol.xyz
set_view (\
     0.020379238,    0.412595391,   -0.910679519,\
    -0.174068779,   -0.895494044,   -0.409611315,\
    -0.984519064,    0.166866705,    0.053572617,\
     0.000005444,    0.000003511,  -14.221484184,\
    -0.747403860,    0.374430537,    1.037632346,\
    11.221332550,   17.239648819,  -20.000000000 )
set stick_radius, 0.06
set sphere_scale, 0.15
show sticks
show spheres
set transparency, 0.9
set opaque_background, off
set ray_opaque_background, off
ray 4000, 4000
png water_methanol, dpi=1200
EOF

