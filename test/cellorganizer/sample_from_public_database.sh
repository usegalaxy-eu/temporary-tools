#!/usr/bin/env bash

export PATH=$PATH:$(dirname $0)

WORKING_DIRECTORY=`pwd`

OMERO_ID=`python -c "from __future__ import print_function; import random; print(random.choice(range(51, 251) + range(304, 425) + range(451, 532)), end='')"`

wget -O output.tif -nc "http://omero.compbio.cs.cmu.edu:8080/webgateway/render_ome_tiff/i/"$OMERO_ID"/"
