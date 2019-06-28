#!/bin/bash
pymol -p <<EOF
load water_water.xyz
set_view (\
     0.077019170,    0.975163639,   -0.207657859,\
    -0.356887996,    0.221440658,    0.907520950,\
     0.930967569,    0.004214230,    0.365079701,\
     0.000000000,    0.000000000,  -14.178204536,\
    -0.027711868,   -0.007276356,   -0.205429554,\
    11.178204536,   17.178203583,  -20.000000000 )
set stick_radius, 0.06
set sphere_scale, 0.15
show sticks
show spheres
set transparency, 0.9
set opaque_background, off
set ray_opaque_background, off
ray 900, 900
png water_water, dpi=600
EOF

