/******************************************************************************
 * Project:  GDAL Utilities
 * Purpose:  GDAL Utilities Public Declarations.
 * Author:   Faza Mahamood, fazamhd at gmail dot com
 *
 * ****************************************************************************
 * Copyright (c) 1998, Frank Warmerdam
 * Copyright (c) 2007-2015, Even Rouault <even.rouault at spatialys.com>
 * Copyright (c) 2015, Faza Mahamood
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 ****************************************************************************/

#ifndef OGR2OGR_LIB_H_INCLUDED
#define OGR2OGR_LIB_H_INCLUDED

#include "gdal.h"
#include "gdal_priv.h"
#include "gdal_alg.h"
#include "ogr_api.h"
#include "ogr_core.h"
#include "ogr_feature.h"
#include "ogr_featurestyle.h"
#include "ogr_geometry.h"
#include "ogr_p.h"
#include "ogr_spatialref.h"
#include "ogrsf_frmts.h"

/* Access modes */
typedef enum
{
    ACCESS_CREATION,
    ACCESS_UPDATE, /* open existing output datasource in update mode rather than trying to create a new one */
    ACCESS_APPEND, /* append to existing layer instead of creating new */
    ACCESS_OVERWRITE /*  delete the output layer and recreate it empty */
} GDALVectorTranslateAccessMode;

struct GDALVectorTranslateOptionsForBinary
{
    char* pszDataSource;
    char* pszDestDataSource;
    int bQuiet;
    char** papszOpenOptions;
    int bFormatExplicitlySet;
    char* pszFormat;
    GDALVectorTranslateAccessMode eAccessMode;
};

/************************************************************************/
/*                     OGR2OGRSpatialReferenceHolder                    */
/************************************************************************/

class OGR2OGRSpatialReferenceHolder
{
        OGRSpatialReference* m_poSRS;

    public:
        OGR2OGRSpatialReferenceHolder() : m_poSRS(NULL) {}
       ~OGR2OGRSpatialReferenceHolder() { if( m_poSRS) m_poSRS->Release(); }

       void assignNoRefIncrease(OGRSpatialReference* poSRS) {
           CPLAssert(m_poSRS == NULL);
           m_poSRS = poSRS;
       }
       OGRSpatialReference* get() { return m_poSRS; }
};

typedef enum
{
    GEOMOP_NONE,
    GEOMOP_SEGMENTIZE,
    GEOMOP_SIMPLIFY_PRESERVE_TOPOLOGY,
} GeomOperation;

typedef enum
{
    GTC_DEFAULT,
    GTC_PROMOTE_TO_MULTI,
    GTC_CONVERT_TO_LINEAR,
    GTC_CONVERT_TO_CURVE,
} GeomTypeConversion;

#define GEOMTYPE_UNCHANGED  -2

#define COORD_DIM_UNCHANGED -1
#define COORD_DIM_LAYER_DIM -2
#define COORD_DIM_XYM -3

typedef struct
{
    OGRLayer *   poSrcLayer;
    GIntBig      nFeaturesRead;
    bool         bPerFeatureCT;
    OGRLayer    *poDstLayer;
    OGRCoordinateTransformation **papoCT; // size: poDstLayer->GetLayerDefn()->GetFieldCount();
    char       ***papapszTransformOptions; // size: poDstLayer->GetLayerDefn()->GetFieldCount();
    int         *panMap;
    int          iSrcZField;
    int          iSrcFIDField;
    int          iRequestedSrcGeomField;
    bool         bPreserveFID;
    GIntBig      nFeaturesOutOfClip, nFeaturesInsideClip, nFeaturesClipped, nFeaturesScipClip;
} TargetLayerInfo;

typedef struct
{
    OGRLayer         *poSrcLayer;
    TargetLayerInfo  *psInfo;
} AssociatedLayers;


/************************************************************************/
/*                        GDALVectorTranslateOptions                    */
/************************************************************************/

/** Options for use with GDALVectorTranslate(). GDALVectorTranslateOptions* must be allocated and
 * freed with GDALVectorTranslateOptionsNew() and GDALVectorTranslateOptionsFree() respectively.
 */
struct GDALVectorTranslateOptions
{
    /*! try to fix invalid geometry */
    bool bFixGeometries;

    /*! continue after a failure, skipping the failed feature */
    bool bSkipFailures;

    /*! use layer level transaction. If set to FALSE, then it is interpreted as dataset level transaction. */
    int nLayerTransaction;

    /*! force the use of particular transaction type based on GDALVectorTranslate::nLayerTransaction */
    bool bForceTransaction;

    /*! group nGroupTransactions features per transaction (default 20000). Increase the value for better
        performance when writing into DBMS drivers that have transaction support. nGroupTransactions can
        be set to -1 to load the data into a single transaction */
    int nGroupTransactions;

    /*! If provided, only the feature with this feature id will be reported. Operates exclusive of
        the spatial or attribute queries. Note: if you want to select several features based on their
        feature id, you can also use the fact the 'fid' is a special field recognized by OGR SQL.
        So GDALVectorTranslateOptions::pszWHERE = "fid in (1,3,5)" would select features 1, 3 and 5. */
    GIntBig nFIDToFetch;

    /*! allow or suppress progress monitor and other non-error output */
    bool bQuiet;

    /*! output file format name (default is ESRI Shapefile) */
    char *pszFormat;

    /*! list of layers of the source dataset which needs to be selected */
    char **papszLayers;

    /*! dataset creation option (format specific) */
    char **papszDSCO;

    /*! layer creation option (format specific) */
    char **papszLCO;

    /*! access modes */
    GDALVectorTranslateAccessMode eAccessMode;

    /*! It has the effect of adding, to existing target layers, the new fields found in source layers.
        This option is useful when merging files that have non-strictly identical structures. This might
        not work for output formats that don't support adding fields to existing non-empty layers. */
    bool bAddMissingFields;

    /*! It must be set to true to trigger reprojection, otherwise only SRS assignment is done. */
    bool bTransform;

    /*! output SRS. GDALVectorTranslateOptions::bTransform must be set to true to trigger reprojection,
        otherwise only SRS assignment is done. */
    char *pszOutputSRSDef;

    /*! override source SRS */
    char *pszSourceSRSDef;

    bool bNullifyOutputSRS;

    /*! If set to false, then field name matching between source and existing target layer is done
        in a more relaxed way if the target driver has an implementation for it. */
    bool bExactFieldNameMatch;

    /*! an alternate name to the new layer */
    char *pszNewLayerName;

    /*! attribute query (like SQL WHERE) */
    char *pszWHERE;

    /*! name of the geometry field on which the spatial filter operates on. */
    char *pszGeomField;

    /*! list of fields from input layer to copy to the new layer. A field is skipped if
        mentioned previously in the list even if the input layer has duplicate field names.
        (Defaults to all; any field is skipped if a subsequent field with same name is
        found.) Geometry fields can also be specified in the list. */
    char **papszSelFields;

    /*! SQL statement to execute. The resulting table/layer will be saved to the output. */
    char *pszSQLStatement;

    /*! SQL dialect. In some cases can be used to use (unoptimized) OGR SQL instead of the
        native SQL of an RDBMS by using "OGRSQL". The "SQLITE" dialect can also be used with
        any datasource. */
    char *pszDialect;

    /*! the geometry type for the created layer */
    int eGType;

    GeomTypeConversion eGeomTypeConversion;

    /*! Geometric operation to perform */
    GeomOperation eGeomOp;

    /*! the parameter to geometric operation */
    double dfGeomOpParam;

    /*! list of field types to convert to a field of type string in the destination layer.
        Valid types are: Integer, Integer64, Real, String, Date, Time, DateTime, Binary,
        IntegerList, Integer64List, RealList, StringList. Special value "All" can be
        used to convert all fields to strings. This is an alternate way to using the CAST
        operator of OGR SQL, that may avoid typing a long SQL query. Note that this does
        not influence the field types used by the source driver, and is only an afterwards
        conversion. */
    char **papszFieldTypesToString;

    /*! list of field types and the field type after conversion in the destination layer.
        ("srctype1=dsttype1","srctype2=dsttype2",...).
        Valid types are : Integer, Integer64, Real, String, Date, Time, DateTime, Binary,
        IntegerList, Integer64List, RealList, StringList. Types can also include subtype
        between parenthesis, such as Integer(Boolean), Real(Float32), ... Special value
        "All" can be used to convert all fields to another type. This is an alternate way to
        using the CAST operator of OGR SQL, that may avoid typing a long SQL query.
        This is a generalization of GDALVectorTranslateOptions::papszFieldTypeToString. Note that this does not influence
        the field types used by the source driver, and is only an afterwards conversion. */
    char **papszMapFieldType;

    /*! set field width and precision to 0 */
    bool bUnsetFieldWidth;

    /*! display progress on terminal. Only works if input layers have the "fast feature count"
    capability */
    bool bDisplayProgress;

    /*! split geometries crossing the dateline meridian */
    bool bWrapDateline;

    /*! offset from dateline in degrees (default long. = +/- 10deg, geometries
    within 170deg to -170deg will be split) */
    double dfDateLineOffset;

    /*! clip geometries when it is set to true */
    bool bClipSrc;

    OGRGeometryH hClipSrc;

    /*! clip datasource */
    char *pszClipSrcDS;

    /*! select desired geometries using an SQL query */
    char *pszClipSrcSQL;

    /*! selected named layer from the source clip datasource */
    char *pszClipSrcLayer;

    /*! restrict desired geometries based on attribute query */
    char *pszClipSrcWhere;

    OGRGeometryH hClipDst;

    /*! destination clip datasource */
    char *pszClipDstDS;

    /*! select desired geometries using an SQL query */
    char *pszClipDstSQL;

    /*! selected named layer from the destination clip datasource */
    char *pszClipDstLayer;

    /*! restrict desired geometries based on attribute query */
    char *pszClipDstWhere;

    /*! split fields of type StringList, RealList or IntegerList into as many fields
        of type String, Real or Integer as necessary. */
    bool bSplitListFields;

    /*! limit the number of subfields created for each split field. */
    int nMaxSplitListSubFields;

    /*! produce one feature for each geometry in any kind of geometry collection in the
        source file */
    bool bExplodeCollections;

    /*! uses the specified field to fill the Z coordinates of geometries */
    char *pszZField;

    /*! the list of field indexes to be copied from the source to the destination. The (n)th value
        specified in the list is the index of the field in the target layer definition in which the
        n(th) field of the source layer must be copied. Index count starts at zero. There must be
        exactly as many values in the list as the count of the fields in the source layer.
        We can use the "identity" option to specify that the fields should be transferred by using
        the same order. This option should be used along with the
        GDALVectorTranslateOptions::eAccessMode = ACCESS_APPEND option. */
    char **papszFieldMap;

    /*! force the coordinate dimension to nCoordDim (valid values are 2 or 3). This affects both
        the layer geometry type, and feature geometries. */
    int nCoordDim;

    /*! destination dataset open option (format specific), only valid in update mode */
    char **papszDestOpenOptions;

    /*! If set to true, does not propagate not-nullable constraints to target layer if they exist
        in source layer */
    bool bForceNullable;

    /*! If set to true, does not propagate default field values to target layer if they exist in
        source layer */
    bool bUnsetDefault;

    /*! to prevent the new default behaviour that consists in, if the output driver has a FID layer
        creation option and we are not in append mode, to preserve the name of the source FID column
        and source feature IDs */
    bool bUnsetFid;

    /*! use the FID of the source features instead of letting the output driver to automatically
        assign a new one. If not in append mode, this behaviour becomes the default if the output
        driver has a FID layer creation option. In which case the name of the source FID column will
        be used and source feature IDs will be attempted to be preserved. This behaviour can be
        disabled by option GDALVectorTranslateOptions::bUnsetFid */
    bool bPreserveFID;

    /*! set it to false to disable copying of metadata from source dataset and layers into target dataset and
        layers, when supported by output driver. */
    bool bCopyMD;

    /*! list of metadata key and value to set on the output dataset, when supported by output driver.
        ("META-TAG1=VALUE1","META-TAG2=VALUE2") */
    char **papszMetadataOptions;

    /*! override spatial filter SRS */
    char *pszSpatSRSDef;

    /*! size of the list GDALVectorTranslateOptions::pasGCPs */
    int nGCPCount;

    /*! list of ground control points to be added */
    GDAL_GCP *pasGCPs;

    /*! order of polynomial used for warping (1 to 3). The default is to select a polynomial
        order based on the number of GCPs */
    int nTransformOrder;

    /*! spatial query extents, in the SRS of the source layer(s) (or the one specified with
        GDALVectorTranslateOptions::pszSpatSRSDef). Only features whose geometry intersects the extents
        will be selected. The geometries will not be clipped unless GDALVectorTranslateOptions::bClipSrc
        is true. */
    OGRGeometryH hSpatialFilter;

    /*! the progress function to use */
    GDALProgressFunc pfnProgress;

    /*! pointer to the progress data variable */
    void *pProgressData;

    /*! Whether layer and feature native data must be transferred. */
    bool bNativeData;

    /*! Maximum number of features, or -1 if no limit. */
    GIntBig nLimit;

    int nClipStep;
};

class SetupTargetLayer
{
public:
    GDALDataset          *m_poSrcDS;
    GDALDataset          *m_poDstDS;
    char                **m_papszLCO;
    OGRSpatialReference  *m_poOutputSRS;
    bool                  m_bNullifyOutputSRS;
    char                **m_papszSelFields;
    bool                  m_bAppend;
    bool                  m_bAddMissingFields;
    int                   m_eGType;
    GeomTypeConversion    m_eGeomTypeConversion;
    int                   m_nCoordDim;
    bool                  m_bOverwrite;
    char                **m_papszFieldTypesToString;
    char                **m_papszMapFieldType;
    bool                  m_bUnsetFieldWidth;
    bool                  m_bExplodeCollections;
    const char           *m_pszZField;
    char                **m_papszFieldMap;
    const char           *m_pszWHERE;
    bool                  m_bExactFieldNameMatch;
    bool                  m_bQuiet;
    bool                  m_bForceNullable;
    bool                  m_bUnsetDefault;
    bool                  m_bUnsetFid;
    bool                  m_bPreserveFID;
    bool                  m_bCopyMD;
    bool                  m_bNativeData;
    bool                  m_bNewDataSource;

    TargetLayerInfo*            Setup(OGRLayer * poSrcLayer,
                                      const char *pszNewLayerName,
                                      GDALVectorTranslateOptions *psOptions,
                                      GIntBig& nTotalEventsDone);
};

class LayerTranslator
{
public:
    GDALDataset                  *m_poSrcDS;
    GDALDataset                  *m_poODS;
    bool                          m_bTransform;
    bool                          m_bWrapDateline;
    CPLString                     m_osDateLineOffset;
    OGRSpatialReference          *m_poOutputSRS;
    bool                          m_bNullifyOutputSRS;
    OGRSpatialReference          *m_poUserSourceSRS;
    OGRCoordinateTransformation  *m_poGCPCoordTrans;
    int                           m_eGType;
    GeomTypeConversion            m_eGeomTypeConversion;
    int                           m_nCoordDim;
    GeomOperation                 m_eGeomOp;
    double                        m_dfGeomOpParam;
    OGRGeometry                  *m_poClipSrc;
    OGRGeometry                  *m_poClipDst;
    bool                          m_bExplodeCollections;
    bool                          m_bNativeData;
    GIntBig                       m_nLimit;

    int                 Translate(OGRFeature* poFeatureIn,
                                  TargetLayerInfo* psInfo,
                                  GIntBig nCountLayerFeatures,
                                  GIntBig* pnReadFeatureCount,
                                  GIntBig& nTotalEventsDone,
                                  GDALProgressFunc pfnProgress,
                                  void *pProgressArg,
                                  GDALVectorTranslateOptions *psOptions);
};

/************************************************************************/
/*                            GCPCoordTransformation()                  */
/*                                                                      */
/*      Apply GCP Transform to points                                   */
/************************************************************************/

class GCPCoordTransformation : public OGRCoordinateTransformation
{
public:

  void               *hTransformArg;
  bool                 bUseTPS;
  OGRSpatialReference* poSRS;

  GCPCoordTransformation( int nGCPCount,
                          const GDAL_GCP *pasGCPList,
                          int  nReqOrder,
                          OGRSpatialReference* poSRSIn)
  {
      if( nReqOrder < 0 )
      {
          bUseTPS = true;
          hTransformArg =
              GDALCreateTPSTransformer( nGCPCount, pasGCPList, FALSE );
      }
      else
      {
          bUseTPS = false;
          hTransformArg =
              GDALCreateGCPTransformer( nGCPCount, pasGCPList, nReqOrder, FALSE );
      }
      poSRS = poSRSIn;
      if( poSRS)
          poSRS->Reference();
  }

  bool IsValid() const { return hTransformArg != NULL; }

  virtual ~GCPCoordTransformation()
  {
      if( hTransformArg != NULL )
      {
          if( bUseTPS )
              GDALDestroyTPSTransformer(hTransformArg);
          else
              GDALDestroyGCPTransformer(hTransformArg);
      }
      if( poSRS)
          poSRS->Dereference();
  }

  virtual OGRSpatialReference *GetSourceCS() override { return poSRS; }
  virtual OGRSpatialReference *GetTargetCS() override { return poSRS; }

  virtual int Transform( int nCount,
                         double *x, double *y, double *z = NULL ) override
  {
      int *pabSuccess = (int *) CPLMalloc(sizeof(int) * nCount );

      bool bOverallSuccess = CPL_TO_BOOL(TransformEx( nCount, x, y, z, pabSuccess ));

      for( int i = 0; i < nCount; ++i )
      {
          if( !pabSuccess[i] )
          {
              bOverallSuccess = false;
              break;
          }
      }

      CPLFree( pabSuccess );

      return bOverallSuccess;
  }

  virtual int TransformEx( int nCount,
                           double *x, double *y, double *z = NULL,
                           int *pabSuccess = NULL ) override
  {
      if( bUseTPS )
          return GDALTPSTransform( hTransformArg, FALSE,
                               nCount, x, y, z, pabSuccess );
      else
          return GDALGCPTransform( hTransformArg, FALSE,
                               nCount, x, y, z, pabSuccess );
  }
};

/************************************************************************/
/*                     OGRSplitListFieldLayer                           */
/************************************************************************/

typedef struct
{
 int          iSrcIndex;
 OGRFieldType eType;
 int          nMaxOccurrences;
 int          nWidth;
} ListFieldDesc;

class OGRSplitListFieldLayer : public OGRLayer
{
 OGRLayer                    *poSrcLayer;
 OGRFeatureDefn              *poFeatureDefn;
 ListFieldDesc               *pasListFields;
 int                          nListFieldCount;
 int                          nMaxSplitListSubFields;

 OGRFeature                  *TranslateFeature(OGRFeature* poSrcFeature);

public:
                              OGRSplitListFieldLayer(OGRLayer* poSrcLayer,
                                                     int nMaxSplitListSubFields);
                     virtual ~OGRSplitListFieldLayer();

 bool                        BuildLayerDefn(GDALProgressFunc pfnProgress,
                                             void *pProgressArg);

 virtual OGRFeature          *GetNextFeature() override;
 virtual OGRFeature          *GetFeature(GIntBig nFID) override;
 virtual OGRFeatureDefn      *GetLayerDefn() override;

 virtual void                 ResetReading() override { poSrcLayer->ResetReading(); }
 virtual int                  TestCapability(const char*) override { return FALSE; }

 virtual GIntBig              GetFeatureCount( int bForce = TRUE ) override
 {
     return poSrcLayer->GetFeatureCount(bForce);
 }

 virtual OGRSpatialReference *GetSpatialRef() override
 {
     return poSrcLayer->GetSpatialRef();
 }

 virtual OGRGeometry         *GetSpatialFilter() override
 {
     return poSrcLayer->GetSpatialFilter();
 }

 virtual OGRStyleTable       *GetStyleTable() override
 {
     return poSrcLayer->GetStyleTable();
 }

 virtual void                 SetSpatialFilter( OGRGeometry *poGeom ) override
 {
     poSrcLayer->SetSpatialFilter(poGeom);
 }

 virtual void                 SetSpatialFilter( int iGeom, OGRGeometry *poGeom ) override
 {
     poSrcLayer->SetSpatialFilter(iGeom, poGeom);
 }

 virtual void                 SetSpatialFilterRect( double dfMinX, double dfMinY,
                                                    double dfMaxX, double dfMaxY ) override
 {
     poSrcLayer->SetSpatialFilterRect(dfMinX, dfMinY, dfMaxX, dfMaxY);
 }

 virtual void                 SetSpatialFilterRect( int iGeom,
                                                    double dfMinX, double dfMinY,
                                                    double dfMaxX, double dfMaxY ) override
 {
     poSrcLayer->SetSpatialFilterRect(iGeom, dfMinX, dfMinY, dfMaxX, dfMaxY);
 }

 virtual OGRErr               SetAttributeFilter( const char *pszFilter ) override
 {
     return poSrcLayer->SetAttributeFilter(pszFilter);
 }
};

// /*! Options for GDALVectorTranslate(). Opaque type */
// typedef struct GDALVectorTranslateOptions GDALVectorTranslateOptions;
//
// /** Opaque type */
// typedef struct GDALVectorTranslateOptionsForBinary GDALVectorTranslateOptionsForBinary;

GDALVectorTranslateOptions *GDALVectorTranslateOptionsNew(char** papszArgv, GDALVectorTranslateOptionsForBinary* psOptionsForBinary);

void GDALVectorTranslateOptionsFree( GDALVectorTranslateOptions *psOptions );

void GDALVectorTranslateOptionsSetProgress( GDALVectorTranslateOptions *psOptions,
                                              GDALProgressFunc pfnProgress,
                                              void *pProgressData );

GDALDatasetH GDALVectorTranslate( const char *pszDest, GDALDatasetH hDstDS, int nSrcCount,
                               GDALDatasetH *pahSrcDS,
                               const GDALVectorTranslateOptions *psOptions, int *pbUsageError );

GDALVectorTranslateOptions* GDALVectorTranslateOptionsClone(const GDALVectorTranslateOptions *psOptionsIn);
OGRGeometry* LoadGeometry( const char* pszDS, const char* pszSQL, const char* pszLyr, const char* pszWhere);
GDALDataset* GDALVectorTranslateCreateCopy(GDALDriver* poDriver, const char* pszDest, GDALDataset* poDS, const GDALVectorTranslateOptions* psOptions);
OGRLayer* GetLayerAndOverwriteIfNecessary(GDALDataset *poDstDS,
                                                 const char* pszNewLayerName,
                                                 bool bOverwrite,
                                                 bool* pbErrorOccurred,
                                                 bool* pbOverwriteActuallyDone);
void FreeTargetLayerInfo(TargetLayerInfo* psInfo);
void ApplySpatialFilter(OGRLayer* poLayer, OGRGeometry* poSpatialFilter,
                        OGRSpatialReference* poSpatSRS,
                        const char* pszGeomField,
                        OGRSpatialReference* poSourceSRS);

GDALVectorTranslateOptions *GDALVectorTranslateOptionsNew(char** papszArgv,
                      GDALVectorTranslateOptionsForBinary* psOptionsForBinary);

#endif /* OGR2OGR_LIB_H_INCLUDED */
