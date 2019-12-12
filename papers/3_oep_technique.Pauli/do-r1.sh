#!/bin/bash
base=paper-r1
si=supplementary-r1

for i in $base $si
do
  pdflatex  ${i}.tex
  bibtex ${i}.aux
  #dvipdf ${i}.dvi
done
