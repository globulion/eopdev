#!/bin/bash
base=paper-r3
si=supplementary-r2

for i in $base #$si
do
  pdflatex  ${i}.tex
  bibtex ${i}.aux
  #dvipdf ${i}.dvi
done
