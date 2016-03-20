#ifndef __HELPER_FUNCTIONS_H__
#define __HELPER_FUNCTIONS_H__

#include <math.h>
#include <vector>
#include <algorithm>
#include <set>

using namespace std;

// Whether a given integer is present in a vector of integers
bool findValVec( const vector<int> &Labels , int i)		
{
	unsigned int n;
	bool found = false;
	for(n=0;n < Labels.size();n++)
		if(Labels[n] == i)
			found = true;

	return found;
}

// returns L2-distance between two 3-vectors
double diffVec( const vector<double> &V1, const vector<double> &V2)	
{
	return pow((V1[0]-V2[0])*(V1[0]-V2[0]) + (V1[1]-V2[1])*(V1[1]-V2[1]) + (V1[2]-V2[2])*(V1[2]-V2[2]), 0.5);
}

// Given a 3-vector, this function finds its closest 3-vector from among a vector of 3-vectors in the L2-distance norm
double minVecD( const vector< vector<double> > &CClusters, const vector<double> &MeanColors)	
{
	double minVal = diffVec(CClusters[0],MeanColors);
	unsigned int i;
	double tmpVal;
	for(i=0;i<CClusters.size();i++)
	{
		tmpVal = diffVec(CClusters[i],MeanColors);
		minVal = (minVal<tmpVal)? minVal : tmpVal;
	}
	return minVal;
}
	


#endif