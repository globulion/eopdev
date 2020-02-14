#!/bin/bash
for xyz in water-dimer-steven-fink_roo_2.8085.psi  \
           water-dimer-steven-fink_roo_3.0085.psi  \
           water-dimer-steven-fink_roo_3.3085.psi 
do
    for basis in "dzp" "dzp-d"
    do
        python3 ./test.py $basis $xyz
        mv -v test.log test_${basis}_${xyz}.log
    done
done
