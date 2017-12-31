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

#include "lwgeom_lib.h"

extern "C" {
    #include "liblwgeom.h"
}

/* Return a POINTARRAY from a GEOSCoordSeq */
POINTARRAY * ptarray_from_GEOSCoordSeq(GEOSContextHandle_t geosContext, const GEOSCoordSequence *cs, char want3d)
{
    uint32_t dims=2;
    uint32_t size, i;
    POINTARRAY *pa;
    POINT4D point;

    if ( ! GEOSCoordSeq_getSize_r(geosContext, cs, &size) )
//		lwerror("Exception thrown");
        return nullptr;

    if ( want3d )
    {
        if ( ! GEOSCoordSeq_getDimensions_r(geosContext, cs, &dims) )
//			lwerror("Exception thrown");
            return nullptr;

        /* forget higher dimensions (if any) */
        if ( dims > 3 ) dims = 3;
    }

    pa = ptarray_construct((dims==3), 0, size);

    for (i=0; i<size; i++)
    {
        GEOSCoordSeq_getX_r(geosContext, cs, i, &(point.x));
        GEOSCoordSeq_getY_r(geosContext, cs, i, &(point.y));
        if ( dims >= 3 ) GEOSCoordSeq_getZ_r(geosContext, cs, i, &(point.z));
        ptarray_set_point4d(pa,i,&point);
    }

    return pa;
}

GEOSCoordSeq ptarray_to_GEOSCoordSeq(GEOSContextHandle_t geosContext, const POINTARRAY *pa)
{
    uint32_t dims = 2;
    uint32_t i;
    const POINT3DZ *p3d = nullptr;
    const POINT2D *p2d;
    GEOSCoordSeq sq;

    if ( FLAGS_GET_Z(pa->flags) )
        dims = 3;

    if ( ! (sq = GEOSCoordSeq_create_r(geosContext, pa->npoints, dims)) )
//		lwerror("Error creating GEOS Coordinate Sequence");
        return nullptr;

    for ( i=0; i < pa->npoints; i++ )
    {
        if ( dims == 3 )
        {
            p3d = getPoint3dz_cp(pa, i);
            p2d = (const POINT2D *)p3d;
        }
        else
        {
            p2d = getPoint2d_cp(pa, i);
        }

        GEOSCoordSeq_setX_r(geosContext, sq, i, p2d->x);
        GEOSCoordSeq_setY_r(geosContext, sq, i, p2d->y);

        if ( dims == 3 )
            GEOSCoordSeq_setZ_r(geosContext, sq, i, p3d->z);
    }
    return sq;
}

static GEOSGeometry* ptarray_to_GEOSLinearRing(GEOSContextHandle_t geosContext, const POINTARRAY *pa, int autofix)
{
    GEOSCoordSeq sq;
    GEOSGeom g;
    POINTARRAY *npa = 0;

    if ( autofix )
    {
        /* check ring for being closed and fix if not */
        if ( ! ptarray_is_closed_2d(pa) )
        {
            npa = ptarray_addPoint(pa, getPoint_internal(pa, 0), FLAGS_NDIMS(pa->flags), pa->npoints);
            pa = npa;
        }
        /* TODO: check ring for having at least 4 vertices */
    }

    sq = ptarray_to_GEOSCoordSeq(geosContext, pa);
    if ( npa ) ptarray_free(npa);
    g = GEOSGeom_createLinearRing_r(geosContext, sq);
    return g;
}
 
/* Return an LWGEOM from a Geometry */
static LWGEOM * GEOS2LWGEOM(GEOSContextHandle_t geosContext, const GEOSGeometry *geom, char want3d)
{
    int type = GEOSGeomTypeId_r(geosContext, geom) ;
    int hasZ;
    int SRID = GEOSGetSRID_r(geosContext, geom);

    /* GEOS's 0 is equivalent to our unknown as for SRID values */
    if ( SRID == 0 ) SRID = SRID_UNKNOWN;

    if ( want3d )
    {
        hasZ = GEOSHasZ_r(geosContext, geom);
        if ( ! hasZ )
        {
            want3d = 0;
        }
    }
    const GEOSCoordSequence *cs;
    POINTARRAY *pa, **ppaa;
    const GEOSGeometry *g;
    LWGEOM **geoms;
    uint32_t i, ngeoms;

    switch (type)
    {
    case GEOS_POINT:
        cs = GEOSGeom_getCoordSeq_r(geosContext, geom);
        if ( GEOSisEmpty_r(geosContext, geom) )
          return (LWGEOM*)lwpoint_construct_empty(SRID, want3d, 0);
        pa = ptarray_from_GEOSCoordSeq(geosContext,cs, want3d);
        return (LWGEOM *)lwpoint_construct(SRID, nullptr, pa);

    case GEOS_LINESTRING:
    case GEOS_LINEARRING:
        if ( GEOSisEmpty_r(geosContext, geom) )
          return (LWGEOM*)lwline_construct_empty(SRID, want3d, 0);

        cs = GEOSGeom_getCoordSeq_r(geosContext, geom);
        pa = ptarray_from_GEOSCoordSeq(geosContext,cs, want3d);
        return (LWGEOM *)lwline_construct(SRID, nullptr, pa);

    case GEOS_POLYGON:
        if ( GEOSisEmpty_r(geosContext, geom) )
          return (LWGEOM*)lwpoly_construct_empty(SRID, want3d, 0);
        ngeoms = GEOSGetNumInteriorRings_r(geosContext, geom);
        ppaa = (POINTARRAY**)lwalloc(sizeof(POINTARRAY *)*(ngeoms+1));
        g = GEOSGetExteriorRing_r(geosContext, geom);
        cs = GEOSGeom_getCoordSeq_r(geosContext, g);
        ppaa[0] = ptarray_from_GEOSCoordSeq(geosContext, cs, want3d);
        for (i=0; i<ngeoms; i++)
        {
            g = GEOSGetInteriorRingN_r(geosContext, geom, i);
            cs = GEOSGeom_getCoordSeq_r(geosContext, g);
            ppaa[i+1] = ptarray_from_GEOSCoordSeq(geosContext, cs,
                                                  want3d);
        }
        return (LWGEOM *)lwpoly_construct(SRID, nullptr,
                                          ngeoms+1, ppaa);

    case GEOS_MULTIPOINT:
    case GEOS_MULTILINESTRING:
    case GEOS_MULTIPOLYGON:
    case GEOS_GEOMETRYCOLLECTION:
        ngeoms = GEOSGetNumGeometries_r(geosContext, geom);
        geoms = nullptr;
        if ( ngeoms )
        {
            geoms = (LWGEOM **)lwalloc(sizeof(LWGEOM *)*ngeoms);
            for (i=0; i<ngeoms; i++)
            {
                g = GEOSGetGeometryN_r(geosContext, geom, i);
                geoms[i] = GEOS2LWGEOM(geosContext, g, want3d);
            }
        }
        return (LWGEOM *)lwcollection_construct(type,
                                                SRID, nullptr, ngeoms, geoms);

    default:
        return nullptr;

    }
}

GEOSGeometry * LWGEOM2GEOS(GEOSContextHandle_t geosContext, const LWGEOM *lwgeom, int autofix, OGRwkbGeometryType eType)
{
    GEOSCoordSeq sq;
    GEOSGeom g, shell;
    GEOSGeom *geoms = nullptr;

    uint32_t ngeoms, i;
    int geostype;

    if (lwgeom_has_arc(lwgeom))
    {
        LWGEOM *lwgeom_stroked = lwgeom_stroke(lwgeom, 32);
        GEOSGeometry *g = LWGEOM2GEOS(geosContext, lwgeom_stroked, autofix, eType);
        lwgeom_free(lwgeom_stroked);
        return g;
    }

    LWPOINT *lwp = nullptr;
    LWPOLY *lwpoly = nullptr;
    LWLINE *lwl = nullptr;
    LWCOLLECTION *lwc = nullptr;

    switch (lwgeom->type)
    {
    case POINTTYPE:
        if(eType != wkbPoint && eType != wkbMultiPoint && eType != wkbGeometryCollection)
            return nullptr;

        lwp = (LWPOINT *)lwgeom;

        if ( lwgeom_is_empty(lwgeom) )
        {
            g = GEOSGeom_createEmptyPolygon_r(geosContext);
        }
        else
        {
            sq = ptarray_to_GEOSCoordSeq(geosContext, lwp->point);
            g = GEOSGeom_createPoint_r(geosContext, sq);
        }
        if ( ! g )
        {
            /* lwnotice("Exception in LWGEOM2GEOS"); */
            return nullptr;
        }
        break;
    case LINETYPE:
        if(eType != wkbLineString && eType != wkbMultiLineString && eType != wkbGeometryCollection)
            return nullptr;

        lwl = (LWLINE *)lwgeom;
        /* TODO: if (autofix) */
        if ( lwl->points->npoints == 1 ) {
            /* Duplicate point, to make geos-friendly */
            lwl->points = ptarray_addPoint(lwl->points,
                                   getPoint_internal(lwl->points, 0),
                                   FLAGS_NDIMS(lwl->points->flags),
                                   lwl->points->npoints);
        }
        sq = ptarray_to_GEOSCoordSeq(geosContext, lwl->points);
        g = GEOSGeom_createLineString_r(geosContext, sq);
        if ( ! g )
        {
            /* lwnotice("Exception in LWGEOM2GEOS"); */
            return nullptr;
        }
        break;

    case POLYGONTYPE:
        if(eType != wkbPolygon && eType != wkbMultiPolygon && eType != wkbGeometryCollection)
            return nullptr;

        lwpoly = (LWPOLY *)lwgeom;
        if ( lwgeom_is_empty(lwgeom) )
        {
            g = GEOSGeom_createEmptyPolygon_r(geosContext);
        }
        else
        {
            shell = ptarray_to_GEOSLinearRing(geosContext, lwpoly->rings[0], autofix);
            if ( ! shell ) return nullptr;
            /*lwerror("LWGEOM2GEOS: exception during polygon shell conversion"); */
            ngeoms = lwpoly->nrings-1;
            if ( ngeoms > 0 )
                geoms = (GEOSGeom*)malloc(sizeof(GEOSGeom)*ngeoms);

            for (i=1; i<lwpoly->nrings; ++i)
            {
                geoms[i-1] = ptarray_to_GEOSLinearRing(geosContext, lwpoly->rings[i], autofix);
                if ( ! geoms[i-1] )
                {
                    --i;
                    while (i) GEOSGeom_destroy_r(geosContext, geoms[--i]);
                    free(geoms);
                    GEOSGeom_destroy_r(geosContext, shell);
                    return nullptr;
                }
                /*lwerror("LWGEOM2GEOS: exception during polygon hole conversion"); */
            }
            g = GEOSGeom_createPolygon_r(geosContext, shell, geoms, ngeoms);
            if (geoms) free(geoms);
        }
        if ( ! g ) return nullptr;
        break;
    case MULTIPOINTTYPE:
    case MULTILINETYPE:
    case MULTIPOLYGONTYPE:
    case COLLECTIONTYPE:
        if ( lwgeom->type == MULTIPOINTTYPE )
            geostype = GEOS_MULTIPOINT;
        else if ( lwgeom->type == MULTILINETYPE )
            geostype = GEOS_MULTILINESTRING;
        else if ( lwgeom->type == MULTIPOLYGONTYPE )
            geostype = GEOS_MULTIPOLYGON;
        else
            geostype = GEOS_GEOMETRYCOLLECTION;

        lwc = (LWCOLLECTION *)lwgeom;

        ngeoms = lwc->ngeoms;
        if ( ngeoms > 0 )
            geoms = (GEOSGeom*)malloc(sizeof(GEOSGeom)*ngeoms);

        for (i=0; i<ngeoms; ++i)
        {
            GEOSGeometry* g = LWGEOM2GEOS(geosContext, lwc->geoms[i], 0, eType);
            if ( ! g )
            {
                while (i) GEOSGeom_destroy_r(geosContext, geoms[--i]);
                free(geoms);
                return nullptr;
            }
            geoms[i] = g;
        }
        g = GEOSGeom_createCollection_r(geosContext, geostype, geoms, ngeoms);
        if ( geoms ) free(geoms);
        if ( ! g ) return nullptr;
        break;

    default:
//		lwerror("Unknown geometry type: %d - %s", lwgeom->type, lwtype_name(lwgeom->type));
        return nullptr;
    }

    GEOSSetSRID_r(geosContext, g, lwgeom->srid);

    return g;
}

GEOSGeom MakeValid(GEOSContextHandle_t geosContext, GEOSGeom invalidGeom, OGRwkbGeometryType eGType) {
    LWGEOM* in = GEOS2LWGEOM(geosContext, invalidGeom, 0);
    LWGEOM* out = lwgeom_make_valid(in);
    if(nullptr == out) {
        return invalidGeom;
    }

    return LWGEOM2GEOS(geosContext, out, 0, wkbFlatten(eGType));
}
