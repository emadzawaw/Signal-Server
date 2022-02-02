#!/bin/sh
################################################################################
# utils/sdf/convert_sdf_hd.sh
#
# Converts from *.hgt files to *.sdf files
#
# HD quality
#
################################################################################

for file in *.hgt
do
	srtm2sdf-hd -d /dev/null $file
done
