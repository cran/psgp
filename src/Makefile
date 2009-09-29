
LIBS = `itpp-config --libs` -L/usr/local/lib/R/lib
CXX =	g++
CXXFLAGS = -I. -I/usr/share/R/include `itpp-config --cflags` -DUSING_R=1 -O3

SRC =	astonGeostats.cpp \
	RInterface.cpp \
	gaussianProcesses/GaussianProcess.cpp \
	gaussianProcesses/SequentialGP.cpp \
	gaussianProcesses/ForwardModel.cpp \
	covarianceFunctions/GaussianCF.cpp \
	covarianceFunctions/ConstantCF.cpp \
	covarianceFunctions/CovarianceFunction.cpp \
	covarianceFunctions/ExponentialCF.cpp \
	covarianceFunctions/WhiteNoiseCF.cpp \
	covarianceFunctions/SumCovarianceFunction.cpp \
	covarianceFunctions/Transform.cpp \
	covarianceFunctions/LogTransform.cpp \
	covarianceFunctions/IdentityTransform.cpp \
	likelihoodModels/AnalyticLikelihood.cpp \
	likelihoodModels/LikelihoodType.cpp \
	likelihoodModels/SamplingLikelihood.cpp \
	likelihoodModels/GaussianLikelihood.cpp \
	likelihoodModels/ExponentialSampLikelihood.cpp \
	likelihoodModels/GaussianSampLikelihood.cpp \
	optimisation/Optimisable.cpp \
	optimisation/ModelTrainer.cpp \
	optimisation/SCGModelTrainer.cpp \


OBJS = $(SRC:.cpp=.o)

TARGET = psgp.so

$(TARGET):
#	R CMD COMPILE $(SRC)
	MAKEFLAGS='CXXFLAGS=\ -O3\ -DUSING_R=1' $(R_HOME)/bin/R CMD SHLIB $(OBJS) -o $(TARGET) $(LIBS)

all: $(SRC) $(TARGET)

clean:
	rm -f $(OBJS) $(TARGET) psgp.so

##.cpp.o:
##	echo "FRED"
##	$(CXX) $(CXXFLAGS) -c $< -o $@

