Dimension Estimation with PH0

DEVELOPER: Benjamin Schweinhart (bschweinhart@albany.edu)
DATE: September 17, 2020
LICENSE: MPL 1.1/GPL 2.0/LGPL 2.1 (see license.txt)


*****AND*****

Dimension Estimation with the Correlation Dimension
DEVELOPER: Jonathan Jaquette (jaquette@bu.edu)
DATE: September 17, 2020
LICENSE: MPL 1.1/GPL 2.0/LGPL 2.1 (see license.txt)



OVERVIEW:

This software estimates the 0-dimensional persistent homology dimension of a data set using minimum spanning trees (MSTs), based on FRACTAL DIMENSION ESTIMATION WITH PERSISTENT
HOMOLOGY: A COMPARATIVE STUDY by Jaquette and Schweinhart. It computes the MST of a point set in Euclidean space using the implementation of the double-tree Boruvka algorithm in mlpack. Then, it computes sum of the lengths of the edges of the MST to the alpha power (the "alpha-weight.") It then fits a power law Y=aX^b to the data (n,alpha weight for an MST on n points). The dimension estimate is then d=alpha/(1.0-a). Please see the article for more information. 

The software can be run from the command line as described below, or by including the header file "dimensionEstimation_PH0.h." The header file includes additional documentation for functions. 

In addition, the MATLAB code file "Corr_Dim_File.m" may be used to estimate the correlation dimension of a data set. For its documentation, see the comments in the file.

This is a "beta" version. Please email the author at bschweinhart@albany.edu if you experience any errors, or have any questions/comments/requests for additional functionality.

This material is  based upon work supported by a NSF Mathematical Sciences Postdoctoral Research Fellowship under award number DMS-1606259.

CITATIONS:

If you use this code, cite the paper, the mlpack library, and the paper introducing the double-tree Boruvka algorithm (March, 2010).

@article{jaquette2020fractal,
  title={Fractal dimension estimation with persistent homology: A comparative study},
  author={Jaquette, Jonathan and Schweinhart, Benjamin},
  journal={Communications in Nonlinear Science and Numerical Simulation},
  volume={84},
  pages={105163},
  year={2020},
  publisher={Elsevier}
}



@Article{2018curtin,
  author  = {R.R. Curtin AND M. Edel AND M. Lozhnikov AND Y. Mentekidis AND S. Ghaisas AND S. Zhang},
  title   = {mlpack 3: a fast, flexible machine learning library},
  journal = {Journal of Open Source Software},
  year    = {2018},
}

@Article{2010march,
  author  = {William B. March AND Parikshit Ram AND Alexander G. Gray},
  title   = {Fast {E}uclidean Minimum Spanning Tree: Algorithm, Analysis, and Applications},
  journal = {KDD '10 Proceedings of the 16th ACM SIGKDD international conference on Knowledge discovery and data mining},
  year    = {2010},
}


DEPENDENCIES:

C++.

mlpack https://www.mlpack.org/.

Armadillo: http://arma.sourceforge.net/.



INSTALLATION:

After installing the dependencies, compile the command line program "estimateDimension.cpp" as follows:

 g++ estimateDimension.cpp dimensionEstimation_PH0.cpp -larmadillo -lmlpack  -o estimatePH0 -std=c++11 -O2 -fopenmp


INPUT FORMAT:

Point cloud data is inputted as a series of trials, each assumed to include the same number of points. The first line corresponding to the j^th trial is ">>Trialj" and the following lines contain the coordinates of the points of the trial. For example:

	Trial0
	0.853183 0.027619
	0.258144 0.447079
	0.701377 0.456282
	0.603918 0.0175923
	0.580954 0.721831
	0.793012 0.032166
	0.265038 0.406934
	0.417551 0.114177
	0.475683 0.657667
	Trial1
	0.728806 0.388454
	0.860009 0.216639
	0.666785 0.462113
	0.470014 0.0476032
	0.739897 0.220043
	0.332839 0.00507469
	0.378793 0.655559
	0.091009 0.13559
	0.922048 0.074543

A larger example containing 5 trials of 50,000 random points sampled from the natural measure on the Sierpinski Triangle is included in sierpinski_50K.


COMMAND LINE, BASIC USAGE: 

Usage getopt -f filename 

To use the command line option, make sure you have compiled "estimateDimension" as described in the installation section. 

Takes as input data in the format described above. Each trial is assumed to have the same number of points (n). Computes PH0 dimension estimates at 100 logarithmically spaced values of n_j between 1000 and n (using alpha=1.0), and outputs the data to filename_PH0Dimension_1.0. Each line of the output file is in the form "n_j, dimension estimate using n_j points" for each of the lograthmically spaced values of n_j after the fifth one. If the input file contains more than one trial, the dimension estimate is averaged over the trials and each line of the output file is of the form "n_j, average dimension estimate using n_j points, standard deviation of dimension estimate using n_j points". For example, running
	./estimatePH0 -f sierpinski_50K
will produce a file  sierpinski_50K_PH0Dimension_1.0 with the first three lines

	1218 1.61575 0.111565
	1267 1.60143 0.0733027
	1318 1.58239 0.0850917


The first line indicates that the dimension estimate using 1218 points is 1.61575, with a standard deviation of 0.111565 between trials. This can be used to generate plots similar to those in the paper FRACTAL DIMENSION ESTIMATION WITH PERSISTENT HOMOLOGY: A COMPARATIVE STUDY. It is recommended that you check this data for oscillations before using the dimension estimate.


COMMAND LINE, ADVANCED USAGE: 

Different options can be selecting by using the following flags.


Usage: getopt -f filename [-a alpha1,alpha2,...] [-n numComputed] [-m maxPoints] [-t maxTrials]


-f: The filename of the input file, in the format described above.

-a: A list of the values of alpha to estimate the dimension with, separated by commas. The default is alpha=1.0. 

-n: The number of values of n_j at which to compute the alpha-weights, which will be logarithmically spaced between 1000 and maxPoints. The default is numComputed=100.

-m: The maximum number of points to use from each trial. Make sure to choose a value smaller than the minimum size of the trials. The default is to use the number of points in the first trial.
 
-t: The maximum number of trials to load. The default is to use all of them.

Running the program will produce a file filename_PH0Dimension_alpha for each value of alpha specified using the -a flag. The output will be in the same format described in the basic usage section.
