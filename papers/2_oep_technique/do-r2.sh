#!/bin/bash
base=paper-r2
si=supplementary-r2

for i in $base $si
do
  latex  ${i}.tex
  bibtex ${i}.aux
  dvipdf ${i}.dvi
done
