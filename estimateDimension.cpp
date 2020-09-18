

/*
Benjamin Schweinhart (2020)

A command-line program for dimension estimation using the PH0 dimension. Please see the readme for documentation.

Compile with
 g++ estimateDimension.cpp dimensionEstimation_PH0.cpp -larmadillo -lmlpack  -o outPH0 -std=c++11 -O2 -fopenmp
*/



#include "dimensionEstimation_PH0.h"



//Parses a string of doubles separated by commas. 
std::vector<double> parseStringDouble(std::string toParse){
	std::vector<double> toReturn;

	std::stringstream ss(toParse);

	while(ss.good()){
		std::string substring;
		getline( ss, substring, ',' );
		toReturn.push_back(atof(substring.c_str()));
	}
	return toReturn;
}


		
int main(int argc, char** argv) {

	std::string infile="";
	int minPoints=1000;
	int maxPoints=std::numeric_limits<int>::max();
	int numComputed=100;
	int maxTrials=std::numeric_limits<int>::max();

	int opt;

	std::vector<double> alphaVect={1.0};
	while ((opt = getopt(argc,argv,"f:a:n:m:t:")) != EOF)
	switch(opt)
	{
		case 'f': infile=optarg; break;
		case 'a': alphaVect=parseStringDouble(optarg) ; break;
		case 'n': numComputed=atoi(optarg); break;
		case 'm': maxPoints=atoi(optarg); break;
		case 't': maxTrials=atoi(optarg); break;
		case '?': fprintf(stderr, "Usage is \n -f : for name of the file to load \n -a: input the values of alpha for which to compute the dimension, separated by commas (default is alpha=1.0) \n -m: for the maximum number of points to use (default is to use the entire data set)\n Please see the readme for more details. \n -t: for the maximum number of trials to use (default is to use them all) \n -n: for how many values of n to compute the alpha-weighted sum (default is 100) ");
	}

	if (infile==""){
		std::cout<<"Please provide a data file. See the readme for documentation."<<std::endl;
		return 0;
	}

	if (alphaVect.size()==0){alphaVect={1.0};}

	
	std::ifstream input_stream(infile);
	std::string line;
	clock_t start=clock();
	int count=0;



	arma::vec nVect;//the values of n at which to compute the alpha-weights. Determined once the size of the first trial is known.
	std::vector<arma::mat > weightData={};//records the alpha weights for each value of n, each trial, and each value of alpha.
	//weightData[i] - data for the i-th trial, an m x p matrix A so A(j,k)=alpha_k weight for n_j points




	std::getline(input_stream,line);
	while (!line.empty()){
		if (count>maxTrials){break;}

	
		//load the data from a trial
		std::cout<<"Loading data from trial "<<count<<"."<<std::endl;

		if ((line[0]!='>') and (line[0]!='<')){
			std::cout<<"Input file is in the wrong format. Please consult the readme."<<std::endl;
			std::cout<<line[0]<<std::endl;
			return 0;
		}
		
		
		std::vector<std::vector<double> > pts={};
		
		bool continue0=true;

		std::getline(input_stream,line);
		while (continue0){


		
			std::vector<double> point;

			double value; 
			std::stringstream s(line);
	
		
			while (s >> value) {
				if (s.fail()){
					std::cout<<"Input file is in the wrong format. Please consult the readme."<<std::endl;
					return 0;
				}
				point.push_back(value);
				s.ignore();
			}
			if (!point.empty()) pts.push_back(point);
			else{
				std::cout<<"Input file is in the wrong format. Please consult the readme."<<std::endl;
				return 0;
			}
			if (point.size() != pts.front().size()){
				std::cout<<"Input file is in the wrong format. Please consult the readme."<<std::endl;
				return 0;

			}
			
			std::getline(input_stream,line);
			
			if (line.empty()){continue0=false;}
			else if (line[0]=='>'){continue0=false;}
		}



		
		int ambientDimension=pts[0].size();

	
		if (count==0){maxPoints=std::min(maxPoints,(int) pts.size());}
		else if (pts.size()<maxPoints){
			std::cout<<"Point count mismatch. This trial contains "<<pts.size()<<" points but we are collecting data for "<<maxPoints<<" points."<<std::endl;
			return 0;
		}

		if (maxPoints==0){
			std::cout<<"No points loaded. Check the filename and the format of the input file. Please consult the readme."<<std::endl;
			return 0;
		}
		std::cout<<"Loaded "<<pts.size()<<" points, computing MST data for "<<maxPoints<<" points."<<std::endl;


		arma::mat pointMatrix(ambientDimension,pts.size());
		for (int j=0;j<pts.size();j++){for (int i=0;i<ambientDimension;i++){pointMatrix(i,j)=pts[j][i];}}


		if (count==0){//determine the values of n at which to compute the MST data
			nVect=computeNVect(minPoints,maxPoints,numComputed);		
		}
		
		//compute the MSTs and the alpha-weights
		arma::mat trialData(numComputed,alphaVect.size());//trialData(i,j)=alpha_j-weight computed with n_i points
		for (int i=0;i<nVect.size();i++){
			arma::vec mstLengths=intervalLengths(pointMatrix.cols(0,nVect(i)-1));
			for (int j=0;j<alphaVect.size();j++){
				trialData(i,j)=alphaWeight(mstLengths,alphaVect[j]);
			}
		}
		weightData.push_back(trialData);
		std::cout<<"MST data computed for trial "<<count<<"."<<std::endl<<std::endl;;

		count++;
	}




	std::cout<<"MST data computed for "<<count<<" trials. Computing dimension estimates."<<std::endl<<std::endl<<std::endl;


	//Compute a dimension estimate for each value of alpha and each value of n_j. Write the resulting estimates to a file.
	
	

	for (int i=0;i<alphaVect.size();i++){


		std::cout<<"Computing dimension estimate with alpha= "<<alphaVect[i]<<"."<<std::endl<<std::endl;;
		std::string outfile=infile+"_PH0Dimension_"+std::to_string(alphaVect[i]);
		std::ofstream fs(outfile);
		double estimateMean;
		double estimateStdDev;
		
		for (int j=4;j<nVect.size();j++){
			arma::vec estimates(weightData.size());
			//compute a dimension estimate for each trial with n_j points
			for (int k=0;k<weightData.size();k++){
				double exponent=fitPowerLaw(nVect,weightData[k].col(i),floor(j/2),j).second;
				estimates[k]=alphaVect[i]/(1.0-exponent);
			}
			estimateMean=arma::mean(estimates);
			estimateStdDev=arma::stddev(estimates);
			fs<<nVect[j]<<" "<<estimateMean<<" "<<estimateStdDev<<std::endl;
			
		}


		std::cout<<"Dimension estimate for alpha= "<<alphaVect[i]<<": "<<estimateMean<<"."<<std::endl;
		std::cout<<"Remember to check for oscillations using the data in the file "<<outfile<<"."<<std::endl<<std::endl;

	
		fs.close();

	}
		
	std::cout<<"Computation complete. Time elapsed: "<<((double) clock()-start) / ((double)CLOCKS_PER_SEC)<<" seconds."<<std::endl;		


}
	
