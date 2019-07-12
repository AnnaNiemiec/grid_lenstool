#/bin/bash

map=kappa.fits		# Input density map
region=SL.reg		# SL region with no grid
th=0.9			# Density threshold for multiscale grid
LL=5			# Max plit number for multiscale grid
gcat=galcat.dat		# Galaxy catalogue
zl=0.35			# Lens redshift
srccat=arclets.dat	# WL sources catalogue

fit_hex.py $map $th -ll $LL
mk_holes.py $map $region
build_hex.py sdens_holes.dat $gcat $zl $srccat




