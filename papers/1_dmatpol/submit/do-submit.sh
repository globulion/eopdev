#!/bin/bash
base=paper-submit

latex  ${base}.tex
bibtex ${base}.aux
latex  ${base}.tex
latex  ${base}.tex
dvipdf ${base}.dvi
