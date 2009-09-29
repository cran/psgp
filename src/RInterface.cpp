
#include <stdlib.h>
#include <stdio.h>

#include "astonGeostats.h"
#include "R.h"
#include "Rmath.h"
#include "Rinternals.h"

using namespace std;

extern "C" {
    SEXP predict(SEXP xData, SEXP yData, SEXP eData, SEXP xPred, SEXP vario, SEXP errorIdx, SEXP sensorIdx, SEXP metaData)
    {

        SEXP meanResult;
        SEXP varResult;
        SEXP ans;
        double range, sill, nugget, bias;
        int xDataLen, yDataLen, eDataLen, xPredLen;
        //int* xs = INTEGER(paramA);
        //int* ys = INTEGER(paramB);
        double *xDataPtr, *yDataPtr, *xPredPtr, *eDataPtr;
        double *v = REAL(vario);
        int *errorPtr, *sensorPtr;
        int model = 0;
        char **metaDataTable;

        int metadataSize = length(metaData);



        // there must be a more obvious way to get the dimensions of a matrix?
        xDataLen = length(xData) / 2;
        yDataLen = length(yData);
        xPredLen = length(xPred) / 2;
        eDataLen = length(eData);

        // we need some default values here, or maybe we should add something in
        // the learnParamters function to get a good start

        model  = int(v[0]);
        range  = v[1];
        sill   = v[2];
        nugget = v[3];
        bias   = v[4];


        PROTECT(meanResult = allocVector(REALSXP, xPredLen));
        PROTECT(varResult = allocVector(REALSXP, xPredLen));
        PROTECT(ans = allocVector(VECSXP, 2));

        xPredPtr = REAL(xPred);
        xDataPtr = REAL(xData);
        yDataPtr = REAL(yData);
        eDataPtr = REAL(eData);
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

        makePredictions(xDataLen, xPredLen, eDataLen, xDataPtr, yDataPtr, xPredPtr, eDataPtr, 
                metadataSize, errorPtr, sensorPtr, metaDataTable,
                REAL(meanResult), REAL(varResult),
                &range, &sill, &nugget, &bias, model);

        SET_VECTOR_ELT(ans, 0, meanResult);
        SET_VECTOR_ELT(ans, 1, varResult);

        UNPROTECT(3);
        return ans;
    }
}

extern "C" {
    SEXP estParam(SEXP xData, SEXP yData, SEXP eData, SEXP vario, SEXP errorIdx, SEXP sensorIdx, SEXP metaData)
    {
        SEXP meanResult;
        SEXP varResult;
        SEXP ans;
        double range, sill, nugget, bias;
        int xDataLen, eDataLen;
        double *xDataPtr, *yDataPtr, *eDataPtr, *rans, *varioPtr;
        int *errorPtr, *sensorPtr;
        int model = 0;
        char **metaDataTable;

        int metadataSize = length(metaData);

        // there must be a more obvious way to get the dimensions of a matrix?
        xDataLen = length(xData) / 2;
        eDataLen = length(eData);
        xDataPtr = REAL(xData);
        yDataPtr = REAL(yData);
        eDataPtr = REAL(eData);
        varioPtr = REAL(vario);
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

        // we need some default values here, or maybe we should add something in
        // the learnParamters function to get a good start
        model  = int(varioPtr[0]);
        range  = varioPtr[1];
        sill   = varioPtr[2];
        nugget = varioPtr[3];	
        bias	 = varioPtr[4];

        learnParameters(xDataLen, eDataLen, xDataPtr, yDataPtr, eDataPtr,
                metadataSize, errorPtr, sensorPtr, 
                metaDataTable, &range, &sill, &nugget, &bias, model);

        PROTECT(ans = allocVector(REALSXP, 5));

        rans = REAL(ans);
        rans[0] = model;
        rans[1] = range;
        rans[2] = sill;
        rans[3] = nugget;
        rans[4] = bias;
        UNPROTECT(1);
        return ans;
    }
}

