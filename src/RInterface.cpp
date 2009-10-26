
#include <stdlib.h>
#include <stdio.h>

#include "astonGeostats.h"
#include "R.h"
#include "Rmath.h"
#include "Rinternals.h"

#define NUM_VARIOGRAM_PARAMETERS 5

using namespace std;

extern "C" {
    SEXP predict(SEXP xData, SEXP yData, SEXP xPred, SEXP psgpPar, 
                 SEXP errorIdx, SEXP sensorIdx, SEXP metaData)
    {

        SEXP meanResult;
        SEXP varResult;
        SEXP ans;
        
        // These are not used by PSGP anymore. Instead, we use the parameters 
        // stored in psgpParameters 
        // double range, sill, nugget, bias;
        int xDataLen, yDataLen, xPredLen;
        double *xDataPtr, *yDataPtr, *xPredPtr;
        
        // PSGP covariance function parameters
        double  *psgpParameters = REAL(psgpPar);
        
        int     *errorPtr, *sensorPtr;
        char    **metaDataTable;

        int metadataSize = length(metaData);
        // The metadata table, if provided, is terminated by an empty line
        // Need to remove it form the size
        if (metadataSize > 0) metadataSize -= 1;

        // there must be a more obvious way to get the dimensions of a matrix?
        xDataLen = length(xData) / 2;
        yDataLen = length(yData);
        xPredLen = length(xPred) / 2;

        PROTECT(meanResult = allocVector(REALSXP, xPredLen));
        PROTECT(varResult = allocVector(REALSXP, xPredLen));
        PROTECT(ans = allocVector(VECSXP, 2));

        xPredPtr = REAL(xPred);
        xDataPtr = REAL(xData);
        yDataPtr = REAL(yData);
        errorPtr = INTEGER(errorIdx);
        sensorPtr= INTEGER(sensorIdx);

        // we need to make a table of pointers to the sensor model strings
        metaDataTable = (char**)calloc(metadataSize, sizeof(char *));

        // casting from a const pointer to a non-const pointer is a bit
        // of nightmare, but I'm not really sure what the other options are?
        for (int i = 0; i < metadataSize; i++ )
        {
            metaDataTable[i] = const_cast<char*>(CHAR(STRING_ELT(VECTOR_ELT(metaData, i),0)));
        }

        /*
        makePredictions(xDataLen, xPredLen, eDataLen, xDataPtr, yDataPtr, xPredPtr, eDataPtr, 
                metadataSize, errorPtr, sensorPtr, metaDataTable,
                REAL(meanResult), REAL(varResult),
                &range, &sill, &nugget, &bias, model);
        */
        printf("-- Making predictions\n");
        makePredictions(xDataLen, xPredLen, xDataPtr, yDataPtr, xPredPtr, 
                        metadataSize, errorPtr, sensorPtr, metaDataTable,
                        REAL(meanResult), REAL(varResult), psgpParameters);
        
        SET_VECTOR_ELT(ans, 0, meanResult);
        SET_VECTOR_ELT(ans, 1, varResult);

        UNPROTECT(3);
        return ans;
    }
}

extern "C" {
    SEXP estParam(SEXP xData, SEXP yData, SEXP vario, SEXP errorIdx, SEXP sensorIdx, SEXP metaData)
    {
        // SEXP meanResult;
        // SEXP varResult;
        SEXP params;
        int xDataLen;
        double *xDataPtr, *yDataPtr, *varioPtr;
        int *errorPtr, *sensorPtr;
        char **metaDataTable;
        

        int metadataSize = length(metaData);
        // The metadata table, if provided, is terminated by an empty line
        // The empty line will be for sensor model - unused at the moment
        // so we simply discard (i.e. remove) it
        if (metadataSize > 0) metadataSize -= 1;
        
        // there must be a more obvious way to get the dimensions of a matrix?
        xDataLen = length(xData) / 2;
        xDataPtr = REAL(xData);
        yDataPtr = REAL(yData);
        varioPtr = REAL(vario);
        errorPtr = INTEGER(errorIdx);
        sensorPtr= INTEGER(sensorIdx);

        // Check that all data arrays are coherent 
        assert(xDataLen == length(yData));
                
        // we need to make a table of pointers to the sensor model strings
        metaDataTable = (char**)calloc(metadataSize, sizeof(char *));

        // casting from a const pointer to a non-const pointer is a bit
        // of nightmare, but I'm not really sure what the other options are?
        for (int i = 0; i < metadataSize; i++ )
        {
            metaDataTable[i] = const_cast<char*>(CHAR(STRING_ELT(VECTOR_ELT(metaData, i),0)));
        }

        // Allocate parameter array
        PROTECT(params = allocVector(REALSXP, NUM_PSGP_PARAMETERS));
        double* psgpParameters = REAL(params);
        UNPROTECT(1);
        
        // Copy current variogram parameters to parameter array
        memcpy(psgpParameters, varioPtr, NUM_VARIOGRAM_PARAMETERS * sizeof(double));
        
        // Estimate parameters.
        // This also updates the parameter values in psgpParameters 
        // and in variogramParameters (if the initial parameters were not
        // valid, i.e. negative...)
        printf("-- Estimating parameters\n");
        learnParameters(xDataLen, xDataPtr, yDataPtr,
                        metadataSize, errorPtr, sensorPtr, 
                        metaDataTable, psgpParameters);
        
        return params;
    }
}

