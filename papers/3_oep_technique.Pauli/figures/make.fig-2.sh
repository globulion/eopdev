#!/bin/bash
figure_dir="./"

gnuplot plot.fig-2.gpl

for plot in fig-2 ; do
    epstopdf ${plot}-inc.eps
    pdflatex ${plot}.tex
    pdftops -eps ${plot}.pdf
    #mv -v ${plot}.pdf ${plot}.eps ${figure_dir}
    #rm -v ${plot}.pdf
    rm -v ${plot}.aux ${plot}.tex ${plot}.log ${plot}-inc*
done
