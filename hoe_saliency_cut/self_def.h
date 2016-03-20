#ifndef __SELF_DEF_H__
#define __SELF_DEF_H__

#include <math.h>
#include <vector>
#include <algorithm>

using namespace std;

//get the mean value of the pixels in a local window
void get_win_mean( const vector< vector<double>>& Colors , int y , int x , int r , int c , vector<double>& mean ) 
{
	//printf("*\n");
	int loc_y[9] = { 0 , 1, -1, -1 , -1,  0,  1 , 1 , 1 } ;
	int loc_x[9] = { 0 , 1,  1, 0 ,  -1, -1, -1 , 0 , 1 } ;

	int count =0 ;
	double ch1=0 , ch2=0 , ch3=0 ;
	for( int i=0 ; i<9 ; i++ )
	{
		if( y+loc_y[i]<= r-1 && y+loc_y[i]>=0 && x+loc_x[i]<= c-1 && x+loc_x[i]>=0 )
		{
			int this_ind ;
			this_ind = (x+loc_x[i])*r+ y+loc_y[i] ;
			count++ ;
			ch1 += Colors[this_ind][0] ;
			ch2 += Colors[this_ind][1] ;
			ch3 += Colors[this_ind][2] ;
		}
	}

	mean[0] = ch1/count ;
	mean[1] = ch2/count ;
	mean[2] = ch3/count ;
}

//get the variance of the pixels in the high-order clique
float get_var( const vector< vector<double>>& Colors_L , const vector< vector<double>>& Colors_R , 
			   int h, int c, int cr , int numRows , int numCols ,
			   vector<double>& mean_l , vector<double>& mean_r )
{
	int loc_y[9] = { 0 , 1, -1, -1 , -1,  0,  1 , 1 , 1 } ;
	int loc_x[9] = { 0 , 1,  1, 0 ,  -1, -1, -1 , 0 , 1 } ;

	int count =0 ;
	float var_1=0 , var_2=0 , var_3=0 ;
	float ave_1=0 , ave_2=0 , ave_3=0 ;

	ave_1 = (mean_l[0] + mean_r[0]) /2.0 ;
	ave_2 = (mean_l[1] + mean_r[1]) /2.0 ;
	ave_3 = (mean_l[2] + mean_r[2]) /2.0 ;

	for( int i=0 ; i<9 ; i++ )
	{
		int ind_l , ind_r ;
		ind_l = (c+loc_x[i]) *numRows + h+loc_y[i] ;
		ind_r = (cr+loc_x[i])*numRows + h+loc_y[i] ;

		var_1 += (Colors_L[ind_l][0] -var_1)*(Colors_L[ind_l][0] -var_1) + (Colors_R[ind_r][0] -var_1)*(Colors_R[ind_r][0] -var_1) ;
		var_2 += (Colors_L[ind_l][1] -var_2)*(Colors_L[ind_l][1] -var_2) + (Colors_R[ind_r][1] -var_2)*(Colors_R[ind_r][1] -var_2) ;
		var_3 += (Colors_L[ind_l][2] -var_3)*(Colors_L[ind_l][2] -var_3) + (Colors_R[ind_r][2] -var_3)*(Colors_R[ind_r][2] -var_3) ;
	}

	var_1 = var_1/17 ;
	var_2 = var_2/17 ;
	var_3 = var_3/17 ;

	return (var_1+var_2+var_3)/3.0 ;
}

//get the number of high-order cliques
int get_num_corres( const vector<int> &Disparity , int numRows , int numCols ) 
{
	int count=0 ;

	for( int c=1 ; c<numCols-1 ; c++ )
	    for( int h=1 ; h<numRows-1 ; h++ )
		{
			int this_pix = c*numRows + h ;
			int this_disp = Disparity[this_pix] ;

			if( c-this_disp >= 1 )
				count ++ ;
	    }
	return count ;
}



#endif