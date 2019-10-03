#!/bin/bash

topcat -stilts plot2plane \
   xpix=544 \
   xlabel=lgMstar ylabel=lgSFR \
   xmin=8.998 xmax=9.102 ymin=-0.61 ymax=1.83 \
   auxmin=1.05 auxmax=2.89 \
   auxvisible=true auxlabel=z \
   legend=true legpos=0.991,0.986 \
   layer_1=Mark \
      in_1='datatable_generated_galaxies.dump.MSscatter0.29.fits' \
      x_1=lgMstar y_1=lgSFR aux_1=z shading_1=aux \
      icmd_1='select z<1.4' \
      shape_1=cross \
      leglabel_1='MS scatter 0.29 z<1.4' \
   layer_2=Mark \
      in_2='datatable_generated_galaxies.dump.MSscatter0.01.fits' \
      x_2=lgMstar y_2=lgSFR aux_2=z shading_2=aux \
      icmd_2='select z<1.4' \
      leglabel_2='MS scatter 0.01 z<1.4' \
   out='plot_scatter_lgMstar_lgSFR.pdf'


