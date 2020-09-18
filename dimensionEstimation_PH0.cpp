//Please see the readme and header file for documentation.

#include "dimensionEstimation_PH0.h"



//Computes numberComputed logarithmically spaced points between minPoints and maxPoints.
arma::vec computeNVect(int minPoints,int maxPoints,int numberComputed){
	arma::vec toReturn(numberComputed);
	double logBase=std::pow(((double)maxPoints)/minPoints,1.0/((double) numberComputed-1));
	

	for (int k=1;k<=numberComputed;k++){
		if (k==numberComputed){toReturn(k-1)=maxPoints;} //at the last step, include all points
		else{toReturn(k-1)=floor(minPoints*pow(logBase,k));} //logarithmically spaced points
	}
		
	return toReturn;
}


//Computes the lengths of the edges of the MST on a point set, using mlpack's implementation of the Double Tree Boruvka algorithm.
arma::vec intervalLengths(arma::mat pts){



	//compute the MST using mlpack's implementation of the Double Tree Boruvka algorithm
	arma::mat mstResults;
	mlpack::emst::DualTreeBoruvka<> dtb(pts);
	dtb.ComputeMST(mstResults);
	return mstResults.row(2).t();
	
}

//Given the lengths of edges of an MST, computes the alpha-weighted sum (the sum of the lengths of the edges to the alpha power).
double alphaWeight(arma::vec lengths, double alpha){
	double toReturn=0;
	
	for (int i=0;i<lengths.size();i++){
		toReturn=toReturn+pow(lengths(i),alpha);
	}
	return toReturn;
}


//Fits a power law to the data (X,Y) between the indices fitMin and fitMax. Returns a and b where the fit is of the form Y=exp(b*log(X)+a).
std::pair<double,double> fitPowerLaw(arma::vec X, arma::vec Y, int fitMin, int fitMax){
	if (((fitMin<0) or (fitMin>=fitMax)) or (fitMax>=X.size())){
		std::cout<<"Fit bounds in fitPowerLaw incorrect."<<std::endl;
		return  std::make_pair(0,0);
	}
	fitMax=std::min((int)X.size()-1,fitMax);

	arma::mat logY(fitMax-fitMin+1,1);
	arma::rowvec logX(fitMax-fitMin+1);

	for (int i=0;i<fitMax-fitMin+1;i++){
		logX(i)=log(X(fitMin+i));
		logY(i,0)=log(Y(fitMin+i));
	}

	mlpack::regression::LinearRegression lr(logX, logY.t());
	arma::vec parameters = lr.Parameters();

	return std::make_pair(parameters(0),parameters(1));

}



/* Estimates the dimension of a data set. Note that if you are using multiple values of alpha, it is faster to compute the intervalLengths once for each value of n, as in estimateDimension.cpp.
The function computes numberComputed MST alpha-weights using the first n_j points for 1000=n1< n2<...<n_{numberComputed}=nMax. If logspacing=true, the n_j's are spaced logarithmically. Otherwise, they
are spaced linearly. A power law Y=exp(a*log(N)+b) is then fitted between n_{numberComputed/2} and n_{numberComputed} (if fitRangeAuto=true) or between n_{fitMin} and n_{fitMax} ( if fitRangeAuto-false). The function returns d, where a=(d-alpha)/d and saves data in two files. outName_fit contains a single line "a b" where a and b are the parameters estimated for the power law.   outName_data contains three columns. In the j^th row, the first element is n_j, the second is the alpha-weight of the MST on n_j points, and the third is exp(a*n_j+b).*/


double estimateDimension(arma::mat pts, double alpha,std::string outName,int numberComputed, int nMax, bool logSpacing,bool fitRangeAuto, int fitMin, int fitMax){
	
	nMax=std::min(nMax,(int) pts.n_cols);

	if (fitRangeAuto){
		fitMax=floor(numberComputed-1);
		fitMin=floor((numberComputed-1)/2);
	}

	arma::mat alphaWeights(numberComputed,1);
	arma::vec nVect(numberComputed);





	//if the points are linearly spaced spaced, this computes the spacing required to obtained the desired number of points between 1000 and the size
	// of the data set
	double freq=((double) nMax-1000.0)/(numberComputed-1.0);

	//if the points are logarithmically spaced, this computes the log base required to obtained the desired number of points between 1000 and the size
	// of the data set
	double logBase=std::pow(((double) nMax)/1000.0,1.0/((double) numberComputed-1));

	std::cout<<"Computing MST data."<<std::endl;

	for (int k=0;k<numberComputed;k++){
		int n;
		if (k==numberComputed-1){n=nMax;} //at the last step, include all points
		else if (logSpacing){n=floor(1000.0*pow(logBase,k+1));} //logarithmically spaced points
		else {n=floor(1000.0+freq*(k)) ;}//linearly spaced points
		
		nVect(k)=n;

		//compute the MST using mlpack's implementation of the Double Tree Boruvka algorithm
		arma::mat mstResults;
		mlpack::emst::DualTreeBoruvka<> dtb(pts.cols(0,n-1));
		dtb.ComputeMST(mstResults);
		alphaWeights(k,0)=alphaWeight(mstResults.row(2).t(),alpha);
	}


	//fit the power law
	std::pair<double,double> parameters=fitPowerLaw(nVect,alphaWeights,fitMin,fitMax); 
	double a=parameters.second;
	double b=parameters.first;
	double d=alpha/(1.0-a);
	std::cout<<"The power law fit is: Y="<<exp(b)<<"X^"<<a<<". The dimension estimate is d= "<<d<<"."<<std::endl;

	//output results
	std::string outName1=outName+"_fit";
	std::ofstream fs1(outName1);
	fs1<<a<<" "<<b<<std::endl;
	fs1.close();

	std::string outName2=outName+"_data";
	std::ofstream fs2(outName2);
	for (int i=0;i<nVect.size();i++){
		fs2<<nVect[i]<<" "<<alphaWeights(i,0)<<" "<<exp(a*log(nVect(i))+b)<<std::endl;
	}
	fs2.close();

	return d;

}



