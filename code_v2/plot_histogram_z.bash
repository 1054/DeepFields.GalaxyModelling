#!/bin/bash

topcat -stilts plot2plane \
   xpix=614 ypix=418 \
   ylog=true xlabel=z ylabel=N fontsize=20 \
   xmin=-0.02 xmax=10.18 ymin=2 ymax=191979 \
   legend=true legpos=0.981,0.938 \
   binsize=0.1 thick=1 \
   layer_1=Histogram \
      in_1=~/Work/AlmaCosmos/Catalogs/Galaxy_Modelling_SIDES/Mock_cat_Bethermin2017.fits \
      x_1=REDSHIFT \
      leglabel_1='Bethermin2017' \
   layer_2=Histogram \
      in_2=datatable_generated_galaxies_with_coordinates.fits \
      x_2=z \
      color_2=blue \
      leglabel_2='dzliu model'
   out=plot_histogram_z.pdf


