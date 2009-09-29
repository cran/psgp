

void learnParameters(int numObs, int numError, double *xData, double *yData, double *errorData, int numMetadata, int *errPtr, int *sensorPtr, char **metaDataTable, double *range, double *sill, double *nugget, double *bias, int model);

void makePredictions(int numObs, int numPred, int numError, double *xData, double *yData, double *errorData, double *xPred, int numMetadata, int *errPtr, int *sensorPtr, char **metaDataTable, double *meanPred, double *varPred, double *range, double *sill, double *nugget, double *bias, int model);

