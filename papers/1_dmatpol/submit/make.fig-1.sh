#!/bin/bash
figure_dir="./"

gnuplot plot.fig-1.gpl

for plot in fig-1 ; do
    epstopdf ${plot}-inc.eps
    pdflatex ${plot}.tex
    pdftops -eps ${plot}.pdf
    #mv -v ${plot}.pdf ${plot}.eps ${figure_dir}
    rm -v ${plot}.pdf
    rm -v ${plot}.aux ${plot}.tex ${plot}.log ${plot}-inc*
done
