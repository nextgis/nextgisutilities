/******************************************************************************
 * Project:  multithreaded fix and cut geometry utility
 * Purpose:  main source file
 * Author:   Dmitry Baryshnikov, dmitry.baryshnikov@nexgis.com
 *
 ******************************************************************************
 * Copyright (C) 2017 NextGIS
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 2 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ****************************************************************************/

#ifndef LWGEOM_LIB_H_INCLUDED
#define LWGEOM_LIB_H_INCLUDED

#include "geos_c.h"
#include "ogr_geometry.h"

extern GEOSGeom MakeValid(GEOSContextHandle_t geosContext, GEOSGeom invalidGeom, OGRwkbGeometryType eGType);

#endif // LWGEOM_LIB_H_INCLUDED
