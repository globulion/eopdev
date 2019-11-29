#!/bin/bash
base=paper-r0

latex  ${base}.tex
bibtex ${base}.aux
dvipdf ${base}.dvi
