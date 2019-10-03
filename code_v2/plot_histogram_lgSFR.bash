#!/bin/bash

topcat -stilts plot2plane \
   xpix=544 \
   xlabel=lgSFR ylabel=N \
   xmin=-0.53 xmax=1.71 ymin=-6.0E-5 ymax=581 \
   legend=true legpos=0.982,0.986 fontsize=18 \
   x=lgSFR binsize=-221 \
   layer_1=Histogram \
      in_1='datatable_generated_galaxies.dump.MSscatter0.29.fits' \
      icmd_1='select z<1.4' \
      color_1=blue \
      leglabel_1='MS scatter 0.29' \
   layer_2=Histogram \
      in_2='datatable_generated_galaxies.dump.MSscatter0.10.fits' \
      icmd_2='select z<1.4' \
      color_2=magenta \
      leglabel_2='MS scatter 0.10' \
   out='plot_histogram_lgSFR.pdf'


