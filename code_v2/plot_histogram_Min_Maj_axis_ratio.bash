#!/bin/bash

topcat -stilts plot2plane \
    xpix=544 \
    xlabel=Min_kpc/Maj_kpc ylabel=N \
    xmin=-1.0E-7 xmax=1.102 ymin=0 ymax=1970 \
    legend=false fontsize=20 \
    layer=Histogram \
        in=datatable_generated_galaxies.dump.MSscatter0.29.nosizescatter.fits \
        x=Min_kpc/Maj_kpc \
        binsize=-59 \
     out='plot_histogram_Min_Maj_axis_ratio.pdf'


