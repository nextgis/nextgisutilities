# NextGIS utilities

Various NextGIS console utilities

![Language](https://img.shields.io/badge/Language-C++-yellow.svg?maxAge=2592000)

# Install

```
sudo apt-get install nextgisutilities-bin
```

# Performance checking

## Normal run

> time ogr2ogr -overwrite -progress -skipfailures -lco ENCODING=UTF-8 -clipsrc ./osm_cut/irk_bnd.shp ./osm_cut/water-polygon-ogr-final_tmp.shp ./osm_cut/water-polygon1.shp
0ERROR 1: Attempt to write non-polygon (LINESTRING) geometry to POLYGON type shapefile.
ERROR 1: Attempt to write non-polygon (LINESTRING) geometry to POLYGON type shapefile.
.ERROR 1: TopologyException: Input geom 0 is invalid: Self-intersection at or near point 106.2761148 56.456976099999999 at 106.2761148 56.456976099999999
..10...20...30...40...50...60...70...80...90...100 - done.

> real	3m15.196s

> user	3m11.837s

> sys	0m2.982s

## Check if cut geometry contains output geometry

> time ./ng_cutter -overwrite -progress -skipfailures -lco ENCODING=UTF-8 -clipsrc /Volumes/Data/tmp/osm_cut/irk_bnd.shp /Volumes/Data/tmp/osm_cut/water-polygon-ogr-final_tmp.shp /Volumes/Data/tmp/osm_cut/water-polygon1.shp
0..10...20...30...40...50...60...70...80...90...100 - done.

> real	0m35.785s

> user	0m35.068s

> sys	0m0.536s

## Multithreaded cut

> time ./ng_cutter -overwrite -progress -skipfailures -lco ENCODING=UTF-8 -clipsrc /Volumes/Data/tmp/osm_cut/irk_bnd.shp /Volumes/Data/tmp/osm_cut/water-polygon-ogr-final_tmp.shp /Volumes/Data/tmp/osm_cut/water-polygon1.shp --config GDAL_NUM_THREADS 8
0...10...20...30...40...50...60...70...80...90...100 - done.

> real	0m10.708s

> user	1m1.093s

> sys	0m15.187s

# License

This program is licensed under GNU GPL v2 or any later version

![License](https://img.shields.io/badge/License-GPL%20v2-blue.svg?maxAge=2592000)

# Version

0.1.0-27

# Commercial support

Need to fix a bug or add a feature to NextGIS Web? We provide custom development
and support for this software. [Contact us](http://nextgis.ru/en/contact/) to
discuss options!

[![http://nextgis.com](http://nextgis.ru/img/nextgis.png)](http://nextgis.com)
