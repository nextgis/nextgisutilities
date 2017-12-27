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

 #include "ogrsf_frmts.h"
 #include "ogr_p.h"
 #include "commonutils.h"
 #include "ogr2ogr_lib.h"
 // #include <vector>
 // #include <algorithm>

////////////////////////////////////////////////////////////////////////////////
// Usage()
///////////////////////////////////////////////////////////////////////////////
static void Usage(const char* pszAdditionalMsg, int bShort = TRUE) CPL_NO_RETURN;

static bool StringCISortFunction(const CPLString& a, const CPLString& b)
{
    return STRCASECMP(a.c_str(), b.c_str()) < 0;
}

static void Usage(const char* pszAdditionalMsg, int bShort)
{
    GDALDriverManager *poDM = GetGDALDriverManager();

    printf("Usage: ng_cutter [--help-general] [-skipfailures] [-append] [-update]\n"
        "               [-select field_list] [-where restricted_where|@filename]\n"
        "               [-progress] [-sql <sql statement>|@filename] [-dialect dialect]\n"
        "               [-f format_name] [-overwrite] [[-dsco NAME=VALUE] ...]\n"
        "               dst_datasource_name src_datasource_name\n"
        "               [-lco NAME=VALUE] [-nln name] \n"
        "               [-nlt type|PROMOTE_TO_MULTI|CONVERT_TO_LINEAR|CONVERT_TO_CURVE]\n"
        "               [layer [layer ...]]\n"
        "\n"
        "Advanced options :\n"
        "               [-fix] [-clipstep value]\n"
        "               [-gt n] [-ds_transaction]\n"
        "               [[-oo NAME=VALUE] ...] [[-doo NAME=VALUE] ...]\n"
        "               [-clipsrc [xmin ymin xmax ymax]|WKT|datasource|spat_extent]\n"
        "               [-clipsrcsql sql_statement] [-clipsrclayer layer]\n"
        "               [-clipsrcwhere expression]\n"
        "               [-clipdst [xmin ymin xmax ymax]|WKT|datasource]\n"
        "               [-clipdstsql sql_statement] [-clipdstlayer layer]\n"
        "               [-clipdstwhere expression]\n"
        "               [-wrapdateline][-datelineoffset val]\n"
    );

    if (bShort)
    {
        printf("\nNote: ng_cutter --long-usage for full help.\n");
        if (pszAdditionalMsg)
            fprintf(stderr, "\nFAILURE: %s\n", pszAdditionalMsg);
        exit(1);
    }

    printf("\n -f format_name: output file format name, possible values are:\n");

    for (int iDriver = 0; iDriver < poDM->GetDriverCount(); iDriver++)
    {
        GDALDriver *poDriver = poDM->GetDriver(iDriver);

        if(CPLTestBool( CSLFetchNameValueDef(poDriver->GetMetadata(), GDAL_DCAP_VECTOR, "FALSE")) ) {
            continue;
        }

        if( CPLTestBool( CSLFetchNameValueDef(poDriver->GetMetadata(), GDAL_DCAP_CREATE, "FALSE") ) )
            printf("     -f \"%s\"\n", poDriver->GetDescription());
    }

    printf( " -append: Append to existing layer instead of creating new if it exists\n"
            " -overwrite: delete the output layer and recreate it empty\n"
            " -update: Open existing output datasource in update mode\n"
            " -progress: Display progress on terminal. Only works if input layers have the \n"
            "                                          \"fast feature count\" capability\n"
            " -select field_list: Comma-delimited list of fields from input layer to\n"
            "                     copy to the new layer (defaults to all)\n"
            " -where restricted_where: Attribute query (like SQL WHERE)\n"
            " -wrapdateline: split geometries crossing the dateline meridian\n"
            "                (long. = +/- 180deg)\n"
            " -datelineoffset: offset from dateline in degrees\n"
            "                (default long. = +/- 10deg,\n"
            "                geometries within 170deg to -170deg will be split)\n"
            " -sql statement: Execute given SQL statement and save result.\n"
            " -dialect value: select a dialect, usually OGRSQL to avoid native sql.\n"
            " -skipfailures: skip features or layers that fail to convert\n"
            " -fix: try to fix invalid geometry\n"
            " -clipstep value: Clip geometry envelope will split on value x value precheck envelopes.\n"
            " -gt n: group n features per transaction (default 20000). n can be set to unlimited\n"
            " -dsco NAME=VALUE: Dataset creation option (format specific)\n"
            " -lco  NAME=VALUE: Layer creation option (format specific)\n"
            " -oo   NAME=VALUE: Input dataset open option (format specific)\n"
            " -doo  NAME=VALUE: Destination dataset open option (format specific)\n"
            " -nln name: Assign an alternate name to the new layer\n"
            " -nlt type: Force a geometry type for new layer.  One of NONE, GEOMETRY,\n"
            "      POINT, LINESTRING, POLYGON, GEOMETRYCOLLECTION, MULTIPOINT,\n"
            "      MULTIPOLYGON, or MULTILINESTRING, or PROMOTE_TO_MULTI or CONVERT_TO_LINEAR.  Add \"25D\" for 3D layers.\n"
            "      Default is type of source layer.\n"
            " -dim dimension: Force the coordinate dimension to the specified value.\n"
            " -fieldmap index1,index2,...: Specifies the list of field indexes to be\n"
            "      copied from the source to the destination. The (n)th value specified\n"
            "      in the list is the index of the field in the target layer definition\n"
            "      in which the n(th) field of the source layer must be copied. Index count\n"
            "      starts at zero. There must be exactly as many values in the list as\n"
            "      the count of the fields in the source layer. We can use the 'identity'\n"
            "      setting to specify that the fields should be transferred by using the\n"
            "      same order. This setting should be used along with the append setting.");

    if (pszAdditionalMsg)
        fprintf(stderr, "\nFAILURE: %s\n", pszAdditionalMsg);

    exit(1);
}

static void Usage(int bShort = TRUE)
{
    Usage(NULL, bShort);
}

/************************************************************************/
/*                 GDALVectorTranslateOptionsForBinaryNew()             */
/************************************************************************/

static GDALVectorTranslateOptionsForBinary *GDALVectorTranslateOptionsForBinaryNew(void)
{
    return (GDALVectorTranslateOptionsForBinary*) CPLCalloc(  1, sizeof(GDALVectorTranslateOptionsForBinary) );
}

/************************************************************************/
/*                  GDALVectorTranslateOptionsForBinaryFree()           */
/************************************************************************/

static void GDALVectorTranslateOptionsForBinaryFree( GDALVectorTranslateOptionsForBinary* psOptionsForBinary )
{
    if( psOptionsForBinary )
    {
        CPLFree(psOptionsForBinary->pszDataSource);
        CPLFree(psOptionsForBinary->pszDestDataSource);
        CSLDestroy(psOptionsForBinary->papszOpenOptions);
        CPLFree(psOptionsForBinary->pszFormat);
        CPLFree(psOptionsForBinary);
    }
}

/* -------------------------------------------------------------------- */
/*                  CheckDestDataSourceNameConsistency()                */
/* -------------------------------------------------------------------- */

static
void CheckDestDataSourceNameConsistency(const char* pszDestFilename,
                                        const char* pszDriverName)
{
    int i;
    char* pszDestExtension = CPLStrdup(CPLGetExtension(pszDestFilename));

    if( EQUAL(pszDriverName, "GMT") )
        pszDriverName = "OGR_GMT";
    CheckExtensionConsistency(pszDestFilename, pszDriverName);

    static const char* apszBeginName[][2] =  { { "PG:"      , "PostgreSQL" },
                                               { "MySQL:"   , "MySQL" },
                                               { "CouchDB:" , "CouchDB" },
                                               { "GFT:"     , "GFT" },
                                               { "MSSQL:"   , "MSSQLSpatial" },
                                               { "ODBC:"    , "ODBC" },
                                               { "OCI:"     , "OCI" },
                                               { "SDE:"     , "SDE" },
                                               { "WFS:"     , "WFS" },
                                               { NULL, NULL }
                                             };

    for(i=0; apszBeginName[i][0] != NULL; i++)
    {
        if (EQUALN(pszDestFilename, apszBeginName[i][0], strlen(apszBeginName[i][0])) &&
            !EQUAL(pszDriverName, apszBeginName[i][1]))
        {
            CPLError(CE_Warning, CPLE_AppDefined,
                    "The target file has a name which is normally recognized by the %s driver,\n"
                    "but the requested output driver is %s. Is it really what you want ?\n",
                    apszBeginName[i][1],
                    pszDriverName);
            break;
        }
    }

    CPLFree(pszDestExtension);
}

////////////////////////////////////////////////////////////////////////////////
// main()
////////////////////////////////////////////////////////////////////////////////

int main( int nArgc, char ** papszArgv )
{
    GDALDatasetH hDS = NULL;
    GDALDatasetH hODS = NULL;
    int bCloseODS = TRUE;
    int bUsageError = FALSE;
    GDALDatasetH hDstDS;
    int nRetCode = 1;
    GDALVectorTranslateOptionsForBinary* psOptionsForBinary;
    GDALVectorTranslateOptions *psOptions;

    /* Check strict compilation and runtime library version as we use C++ API */
//    if (! GDAL_CHECK_VERSION(papszArgv[0]))
//        exit(1);

    EarlySetConfigOptions(nArgc, papszArgv);

/* -------------------------------------------------------------------- */
/*      Register format(s).                                             */
/* -------------------------------------------------------------------- */
    OGRRegisterAll();

/* -------------------------------------------------------------------- */
/*      Processing command line arguments.                              */
/* -------------------------------------------------------------------- */
    nArgc = OGRGeneralCmdLineProcessor( nArgc, &papszArgv, 0 );

    if( nArgc < 1 )
    {
        papszArgv = NULL;
        nRetCode = -nArgc;
        goto exit;
    }

    for( int iArg = 1; iArg < nArgc; iArg++ )
    {
        if( EQUAL(papszArgv[iArg], "--utility_version") )
        {
            printf("%s was compiled against GDAL %s and is running against GDAL %s\n",
                   papszArgv[0], GDAL_RELEASE_NAME, GDALVersionInfo("RELEASE_NAME"));
            nRetCode = 0;
            goto exit;
        }
        else if( EQUAL(papszArgv[iArg],"--help") )
        {
            Usage();
            goto exit;
        }
        else if ( EQUAL(papszArgv[iArg], "--long-usage") )
        {
            Usage(FALSE);
            goto exit;
        }
    }

    psOptionsForBinary = GDALVectorTranslateOptionsForBinaryNew();
    psOptions = GDALVectorTranslateOptionsNew(papszArgv + 1, psOptionsForBinary);

    if( psOptions == NULL )
    {
        Usage();
        GDALVectorTranslateOptionsForBinaryFree(psOptionsForBinary);
        goto exit;
    }

    if( psOptionsForBinary->pszDataSource == NULL ||
        psOptionsForBinary->pszDestDataSource == NULL )
    {
        if( psOptionsForBinary->pszDestDataSource == NULL )
            Usage("no target datasource provided");
        else
            Usage("no source datasource provided");
        GDALVectorTranslateOptionsFree(psOptions);
        GDALVectorTranslateOptionsForBinaryFree(psOptionsForBinary);
        goto exit;
    }

    if( strcmp(psOptionsForBinary->pszDestDataSource, "/vsistdout/") == 0 )
        psOptionsForBinary->bQuiet = TRUE;

    if (!psOptionsForBinary->bQuiet && !psOptionsForBinary->bFormatExplicitlySet &&
        psOptionsForBinary->eAccessMode == ACCESS_CREATION)
    {
        CheckDestDataSourceNameConsistency(psOptionsForBinary->pszDestDataSource,
                                           psOptionsForBinary->pszFormat);
    }
/* -------------------------------------------------------------------- */
/*      Open data source.                                               */
/* -------------------------------------------------------------------- */

    /* Avoid opening twice the same datasource if it is both the input and output */
    /* Known to cause problems with at least FGdb, SQlite and GPKG drivers. See #4270 */
    if (psOptionsForBinary->eAccessMode != ACCESS_CREATION &&
        strcmp(psOptionsForBinary->pszDestDataSource, psOptionsForBinary->pszDataSource) == 0)
    {
        hODS = GDALOpenEx( psOptionsForBinary->pszDataSource,
                GDAL_OF_UPDATE | GDAL_OF_VECTOR, NULL, psOptionsForBinary->papszOpenOptions, NULL );
        GDALDriverH hDriver = NULL;
        if( hODS != NULL )
            hDriver = GDALGetDatasetDriver(hODS);

        /* Restrict to those 3 drivers. For example it is known to break with */
        /* the PG driver due to the way it manages transactions... */
        if (hDriver && !(EQUAL(GDALGetDescription(hDriver), "FileGDB") ||
                         EQUAL(GDALGetDescription(hDriver), "SQLite") ||
                         EQUAL(GDALGetDescription(hDriver), "GPKG")))
        {
            hDS = GDALOpenEx( psOptionsForBinary->pszDataSource,
                        GDAL_OF_VECTOR, NULL, psOptionsForBinary->papszOpenOptions, NULL );
        }
        else
        {
            hDS = hODS;
            bCloseODS = FALSE;
        }
    }
    else
    {
        hDS = GDALOpenEx( psOptionsForBinary->pszDataSource,
                        GDAL_OF_VECTOR, NULL, psOptionsForBinary->papszOpenOptions, NULL );
    }

/* -------------------------------------------------------------------- */
/*      Report failure                                                  */
/* -------------------------------------------------------------------- */
    if( hDS == NULL )
    {
        GDALDriverManager *poDM = GetGDALDriverManager();

        fprintf( stderr, "FAILURE:\n"
                "Unable to open datasource `%s' with the following drivers.\n",
                psOptionsForBinary->pszDataSource );

        for( int iDriver = 0; iDriver < poDM->GetDriverCount(); iDriver++ )
        {
            GDALDriver* poIter = poDM->GetDriver(iDriver);
            char** papszDriverMD = poIter->GetMetadata();
            if( CPLTestBool( CSLFetchNameValueDef(papszDriverMD, GDAL_DCAP_VECTOR, "FALSE") ) )
            {
                fprintf( stderr,  "  -> `%s'\n", poIter->GetDescription() );
            }
        }

        GDALVectorTranslateOptionsFree(psOptions);
        GDALVectorTranslateOptionsForBinaryFree(psOptionsForBinary);
        goto exit;
    }

    if( hODS != NULL )
    {
        GDALDriverManager *poDM = GetGDALDriverManager();

        GDALDriver* poDriver = poDM->GetDriverByName(psOptionsForBinary->pszFormat);
        if( poDriver == NULL )
        {
            fprintf( stderr,  "Unable to find driver `%s'.\n", psOptionsForBinary->pszFormat );
            fprintf( stderr,  "The following drivers are available:\n" );

            for( int iDriver = 0; iDriver < poDM->GetDriverCount(); iDriver++ )
            {
                GDALDriver* poIter = poDM->GetDriver(iDriver);
                char** papszDriverMD = poIter->GetMetadata();
                if( CPLTestBool( CSLFetchNameValueDef(papszDriverMD, GDAL_DCAP_VECTOR, "FALSE") ) &&
                    (CPLTestBool( CSLFetchNameValueDef(papszDriverMD, GDAL_DCAP_CREATE, "FALSE") ) ||
                     CPLTestBool( CSLFetchNameValueDef(papszDriverMD, GDAL_DCAP_CREATECOPY, "FALSE") )) )
                {
                    fprintf( stderr,  "  -> `%s'\n", poIter->GetDescription() );
                }
            }
            GDALVectorTranslateOptionsFree(psOptions);
            GDALVectorTranslateOptionsForBinaryFree(psOptionsForBinary);
            goto exit;
        }
    }

    if( !(psOptionsForBinary->bQuiet) )
    {
        GDALVectorTranslateOptionsSetProgress(psOptions, GDALTermProgress, NULL);
    }

    hDstDS = GDALVectorTranslate(psOptionsForBinary->pszDestDataSource, hODS,
                                              1, &hDS, psOptions, &bUsageError);
    if( bUsageError )
        Usage();
    else
        nRetCode = (hDstDS) ? 0 : 1;

    GDALVectorTranslateOptionsFree(psOptions);
    GDALVectorTranslateOptionsForBinaryFree(psOptionsForBinary);

    if(hDS)
        GDALClose(hDS);
    if(bCloseODS)
        GDALClose(hDstDS);

exit:
    CSLDestroy( papszArgv );
    OGRCleanupAll();

    return nRetCode;
}
