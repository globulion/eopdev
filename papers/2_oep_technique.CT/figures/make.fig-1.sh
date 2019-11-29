#!/bin/bash
figure_dir="./"

gnuplot plot.fig-1.gpl

for plot in fig-1 ; do
    epstopdf ${plot}-inc.eps
    sed -f ${plot}.sed ${plot}.tex > ${plot}-s.tex
    pdflatex ${plot}-s.tex
    pdftops -eps -r 600 ${plot}-s.pdf
    #mv -v ${plot}.pdf ${plot}.eps ${figure_dir}
    mv -v ${plot}-s.pdf ${plot}.pdf
    rm -v ${plot}-s.aux ${plot}-s.tex ${plot}-s.log ${plot}-inc*
    #rm -v ${plot}.aux ${plot}.tex ${plot}.log 
    mv -v ${plot}-s.eps ${plot}.eps
done
