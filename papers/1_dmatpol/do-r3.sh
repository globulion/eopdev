#!/bin/bash
base=paper-r3

latex  ${base}.tex
bibtex ${base}.aux
dvipdf ${base}.dvi
