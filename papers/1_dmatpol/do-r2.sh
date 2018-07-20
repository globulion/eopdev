#!/bin/bash
base=paper-r2

latex  ${base}.tex
bibtex ${base}.aux
dvipdf ${base}.dvi
