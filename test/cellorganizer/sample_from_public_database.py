#!/usr/bin/env python

import random
import urllib

IFRAME_HTML = '<div id="omeroviewport"><iframe width="850" height="600" src="http://omero.compbio.cs.cmu.edu:8080/webclient/img_detail/{}/" id="omeroviewport" name="omeroviewport"></iframe></div>'

ids = range(51, 251) + range(304, 425) + range(451, 532)
omero_id = random.choice(ids)
html = IFRAME_HTML.format(omero_id)
with open('index.html', 'w') as f:
	f.write(html)

urllib.urlretrieve ("http://omero.compbio.cs.cmu.edu:8080/webgateway/render_ome_tiff/i/" + str(omero_id) + "/", "output.tif")

import glob
print glob.glob("*")
