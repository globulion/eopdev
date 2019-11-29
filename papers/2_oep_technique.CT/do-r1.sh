#!/bin/bash
base=paper-r1

latex  ${base}.tex
bibtex ${base}.aux
dvipdf ${base}.dvi
