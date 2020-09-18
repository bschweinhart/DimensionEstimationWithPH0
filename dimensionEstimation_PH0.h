/*
Benjamin Schweinhart (2020)

This and the associated .cpp file implement the estimation of the PH0 dimension, as defined in FRACTAL DIMENSION ESTIMATION WITH PERSISTENT
HOMOLOGY: A COMPARATIVE STUDY by Jaquette and Schweinhart. This header file includes functions to compute the minimum spanning tree of a point
set using the implemenation of the double tree Boruvka algorithm in mlpack (intervalLengths), the alpha-weighted edge sum of a minimum spanning
tree (alphaWeight), power law fit to data (fitPowerLaw), and the PH0 dimension of a data set (estimateDimension).

Please see the readme for additional documentation, and estimateDimension.cpp for an example.

*/


#include <mlpack/core.hpp>
#include<mlpack/methods/emst/dtb.hpp>
#include <mlpack/methods/linear_regression/linear_regression.hpp>


//Computes numberComputed logarithmically spaced points between minPoints and maxPoints.
arma::vec computeNVect(int minPoints,int maxPoints,int numberComputed);


//Computes the lengths of the edges of the MST on a point set, using mlpack's implementation of the Double Tree Boruvka algorithm.
arma::vec intervalLengths(arma::mat pts);

double alphaWeight(arma::vec lengths, double alpha);


//Fits a power law to the data (X,Y) between the indices fitMin and fitMax. Returns a and b where the fit is of the form Y=exp(b*log(X)+a).
std::pair<double,double> fitPowerLaw(arma::vec X, arma::vec Y, int fitMin, int fitMax);



/* 
Estimates the dimension of a data set. Note that if you are using multiple values of alpha, it is faster to compute the intervalLengths once for
each value of n, as in estimateDimension.cpp. The function computes numberComputed MST alpha-weights using the first n_j points for 
1000=n1< n2<...<n_{numberComputed}=nMax. If logspacing=true, the n_j's are spaced logarithmically. Otherwise, they are spaced linearly. A power 
law Y=exp(a*log(N)+b) is then fitted between n_{numberComputed/2} and n_{numberComputed} (if fitRangeAuto=true) or between n_{fitMin} and n_{fitMax} 
( if fitRangeAuto-false). The function returns d, where a=(d-alpha)/d and saves data in two files. outName_fit contains a single line "a b" where a 
and b are the parameters estimated for the power law.   outName_data contains three columns. In the j^th row, the first element is n_j, the second is 
the alpha-weight of the MST on n_j points, and the third is exp(a*n_j+b).
*/


double estimateDimension(arma::mat pts, double alpha=1.0,std::string outName="PH0Estimate",int numberComputed=100, int nMax=std::numeric_limits<int>::max(), bool logSpacing=true,bool fitRangeAuto=true, int fitMin=1, int fitMax=1000);


