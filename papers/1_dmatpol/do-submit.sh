#!/bin/bash
base=paper-submit

latex  ${base}.tex
bibtex ${base}.aux
dvipdf ${base}.dvi
