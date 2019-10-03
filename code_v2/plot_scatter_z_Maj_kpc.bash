#!/bin/bash

topcat -stilts plot2plane \
   xpix=544 \
   ylog=true xlabel=z ylabel=Maj_kpc \
   xmin=1 xmax=2.92 ymin=0 ymax=22.5 \
   legend=true legpos=0.991,0.986 \
   x=z y=Maj_kpc shading=auto \
   layer_1=Mark \
      in_1='datatable_generated_galaxies.dump.MSscatter0.29.fits' \
      leglabel_1='MS, with size scatter' \
   layer_2=Mark \
      in_2='datatable_generated_galaxies.dump.MSscatter0.29.nosizescatter.fits' \
      color_2=blue \
      leglabel_2='MS, no size scatter' \
   out='plot_scatter_z_Maj_kpc.pdf'


