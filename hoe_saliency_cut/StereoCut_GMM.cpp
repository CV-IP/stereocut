/**************************************************************
GRAPHCUTSEGMENT.CPP - Main graph cuts segmentation C++ function
Authors - Mohit Gupta, Krishnan Ramnath
Affiliation - Robotics Institute, CMU, Pittsburgh
2006-05-15
***************************************************************/

#include <stdio.h>
#include "mex.h"
#include "energy.h"
#include "helper_functions.h"
#include "self_def.h"

 void MakeAdjacencyList( int numLabels, int numCols, int numRows, vector< set<int> > &neighbors )
{
	int i,j;
	
	neighbors.clear() ;

	// Declaring the Adjacency list
	for( i=0 ; i<2*numLabels ; i++ )
	{
		set<int> tmp;
		tmp.clear();
		neighbors.push_back(tmp);
	}
	printf("Adjacency code starts... \n");

	//Filling up the adjacency list 
	for( j=0 ; j<numCols ; j++ )
	    for( i=0 ; i<numRows ; i++ )
		{
			int this_pix = j*numRows + i ;

			if ( i>0 )              //up
			    neighbors[this_pix].insert( j*numRows + (i-1) );
			if ( i<(numRows-1) )    //down
			    neighbors[this_pix].insert( j*numRows + (i+1) );
			if ( j>0 )              //left
			    neighbors[this_pix].insert( (j-1)*numRows + i );
			if ( j<(numCols-1) )    //right
			    neighbors[this_pix].insert( (j+1)*numRows + i );
	    }

	for( j=0 ; j<numCols ; j++ )
	    for( i=0 ; i<numRows ; i++ )
		{
			int this_pix = j*numRows + i ;

			if ( i>0 )              //up
			    neighbors[numLabels+ this_pix].insert( numLabels+ j*numRows + (i-1) );
			if ( i<(numRows-1) )    //down
			    neighbors[numLabels+ this_pix].insert( numLabels+ j*numRows + (i+1) );
			if ( j>0 )              //left
			    neighbors[numLabels+ this_pix].insert( numLabels+ (j-1)*numRows + i );
			if ( j<(numCols-1) )    //right
			    neighbors[numLabels+ this_pix].insert( numLabels+ (j+1)*numRows + i );
	    }
}

void MakeGraph( Energy<float,float,float> *e , int *vars , vector< set<int> > &neighbors,      vector<int> &Disparity,
				vector< vector<double> > Colors_L ,
				vector<double> FDistVec_Uc_L , vector<double> BDistVec_Uc_L ,
				vector<double> FDistVec_Us_L , vector<double> BDistVec_Us_L ,
				vector< vector<double> > Colors_R, 
				vector<double> FDistVec_Uc_R , vector<double> BDistVec_Uc_R ,
				vector<double> FDistVec_Us_R , vector<double> BDistVec_Us_R ,		  
				int numLabels , int num_high , int numRows , int numCols , 
				double lambda_s , double lambda_p , double lambda_c )
{
	int i;
	set<int>::iterator pIter;

	double beta;

	double countCand = 0;
	double SumColors = 0;

	for( i=0 ; i<numLabels ; i++ )
	{
		for ( pIter = neighbors[i].begin() ; pIter != neighbors[i].end() ; pIter++ )
		{
			SumColors += pow( diffVec(Colors_L[i],Colors_L[*pIter]) , 2 );
			countCand++;
		}
	}
	for( i=numLabels ; i<2*numLabels ; i++ )
	{
		for ( pIter = neighbors[i].begin() ; pIter != neighbors[i].end() ; pIter++ )
		{
			SumColors += pow( diffVec(Colors_R[i-numLabels],Colors_R[*pIter-numLabels]) , 2 );
			countCand++;
		}
	}
	beta = countCand/(2*SumColors);


	//when reducing the order of high-order term, every high-order clique needs 5 new auxiliary variables 
	//calulate the factors for data term
	vector<double> ForeEdges, BackEdges;
	ForeEdges.clear();
	BackEdges.clear();
	
	printf("Graph Making code starts... \n");

	/****** Making Terminal Edge Weights ******/
	for( i=0; i<numLabels; i++)
	{
		ForeEdges.push_back(BDistVec_Uc_L[i] + lambda_s*BDistVec_Us_L[i]);
		BackEdges.push_back(FDistVec_Uc_L[i] + lambda_s*FDistVec_Us_L[i]);
	}

	for( i=0; i<numLabels; i++)
	{
		ForeEdges.push_back(BDistVec_Uc_R[i] + lambda_s*BDistVec_Us_R[i]);
		BackEdges.push_back(FDistVec_Uc_R[i] + lambda_s*FDistVec_Us_R[i]);
	}

	printf("Terminal Edge Weights Made... \n");

	///////////////////////////////////////////////
	/******* Start making the graph **************/
	///////////////////////////////////////////////	

	// Add Nodes
	for( i=0 ; i<2*numLabels+4*num_high*5 ; i++ )
		vars[i] = e -> add_variable();

	printf("Nodes Added to Graph... \n");

	//data term
	for( i=0 ; i<2*numLabels ; i++ )
	{
		e -> add_term1(vars[i], ForeEdges[i], BackEdges[i]);
	}

	printf("Terminal Edge Weights Set... \n");

	//smooth term
	for( i=0 ; i<numLabels ; i++ )
		for (pIter = neighbors[i].begin(); pIter != neighbors[i].end(); pIter++)
		{
			int tmpN = *pIter;
			double Energ = exp( -beta * pow(diffVec( Colors_L[i], Colors_L[tmpN]) ,2 ) );
			e -> add_term2( vars[i], vars[tmpN] , 0 , lambda_p*Energ, lambda_p*Energ , 0);
		}
	for( i=numLabels ; i<2*numLabels ; i++ )
		for (pIter = neighbors[i].begin(); pIter != neighbors[i].end(); pIter++)
		{
			int tmpN = *pIter;
			double Energ = exp( -beta * pow(diffVec( Colors_R[i-numLabels], Colors_R[tmpN-numLabels]) ,2 ) );
			e -> add_term2( vars[i], vars[tmpN] , 0 , lambda_p*Energ, lambda_p*Energ , 0);
		}

	//correspondence term: the low-order terms decomposed from the high-order term
	int count_high =0 ;
	int loc_y[4] = { -1 , 1,  0,  0 } ;
	int loc_x[4] = { 0 ,  0, -1,  1 } ;  //up down left right

	float f1 , f2 , a ;
	f2 = 2 ;
	float theta1 , theta2 , gamma1, gamma2 ;

	theta1=0.5 ;
	theta2=1.5 ;  ////////////////////////////////
	////////////////////////////////////////
	gamma1 = 200 ;
	gamma2 = 1000 ;

	for( int c=1 ; c<numCols-1 ; c++ )
	    for( int h=1 ; h<numRows-1 ; h++ )
		{
		    int this_ind , this_disp , cr ;
		    this_ind = c*numRows + h ;
			this_disp = Disparity[this_ind] ;
			cr = c- this_disp ;

		    if( cr >= 1 )
		    {
		    	vector<double> Meanwin_l(3) ;
				vector<double> Meanwin_r(3) ;
				get_win_mean( Colors_L , h , c , numRows , numCols , Meanwin_l ) ;
				get_win_mean( Colors_R , h , cr, numRows , numCols , Meanwin_r ) ;

				double C_rl = exp(    -sqrt(  1/(2*16*16) * pow(diffVec( Meanwin_l, Meanwin_r) ,2)  )    );
				C_rl *= lambda_c ; 

				float local_var=0;
				local_var = get_var( Colors_L , Colors_R , h, c, cr , numRows , numCols , Meanwin_l , Meanwin_r);

				if ( local_var<theta1 )
	                f1 = gamma1 ;
				else if ( local_var>theta2 )
					f1 = gamma2 ;
				else
				{
					f1 = -(gamma2-gamma1)/(theta2-theta1)*local_var + gamma2*theta2-gamma1*theta1 ;
				}

				a = f2 - 2*f1 ;

				if(a>0)
				{
					//one order terms transformed from the high-order term
					e -> add_term1(vars[2*numLabels + count_high*5]   , 0, -C_rl*2*a);
					e -> add_term1(vars[2*numLabels + count_high*5 +1], 0, -C_rl*2*a);
					e -> add_term1(vars[2*numLabels + count_high*5 +2], 0, -C_rl*2*a);
					e -> add_term1(vars[2*numLabels + count_high*5 +3], 0, -C_rl*2*a);
					e -> add_term1(vars[2*numLabels + count_high*5 +4], 0, -C_rl*12*a);    // w1 w2 w3 w4 w

					//two order terms transformed from the high-order term
					for( int u=0 ; u<4 ; u++ )
					{
						//find the corresponding pixels and their neighbours
						int  neigh_y , neigh_xl , neigh_xr , ind_l , ind_r , neigh_ind_l , neigh_ind_r ;
						neigh_y   = h  + loc_y[u] ;
						neigh_xl  = c  + loc_x[u] ;
						neigh_xr  = cr + loc_x[u] ;

						ind_l = c  * numRows + h ;
						ind_r = cr * numRows + h ;
						neigh_ind_l = neigh_xl * numRows + neigh_y ;
						neigh_ind_r = neigh_xr * numRows + neigh_y ;
	
						e -> add_term1(vars[ind_l]                 , 0, C_rl*(f1+6*a) );
						e -> add_term1(vars[ind_r+numLabels]       , 0, C_rl*(f1+6*a) );
						e -> add_term1(vars[neigh_ind_l]           , 0, C_rl*(f1+6*a) );
						e -> add_term1(vars[neigh_ind_r+numLabels] , 0, C_rl*(f1+6*a) );         //x0 y0 x1 y1
	
						e -> add_term2( vars[ind_l]          , vars[neigh_ind_l]           , 0 , 0 , 0 , -C_rl*3*a);   //x0x1 y0y1 x0y1 x1y0
						e -> add_term2( vars[ind_r+numLabels], vars[neigh_ind_r+numLabels] , 0 , 0 , 0 , -C_rl*3*a);
						e -> add_term2( vars[ind_l]          , vars[neigh_ind_r+numLabels] , 0 , 0 , 0 , -C_rl*3*a);
						e -> add_term2( vars[neigh_ind_l]    , vars[ind_r+numLabels]       , 0 , 0 , 0 , -C_rl*3*a);

						e -> add_term2( vars[ind_l]      , vars[ind_r+numLabels]       , 0 , 0 , 0 , C_rl*(-2*f1-4*a) );
						e -> add_term2( vars[neigh_ind_l], vars[neigh_ind_r+numLabels] , 0 , 0 , 0 , C_rl*(-2*f1-4*a) );   //x0y0 x1y1

						e -> add_term2( vars[2*numLabels + count_high*5] , vars[ind_l]           , 0 , 0 , 0 , -C_rl*2*a );  //w1 *( x0 + x1 + y0 )
						e -> add_term2( vars[2*numLabels + count_high*5] , vars[neigh_ind_l]     , 0 , 0 , 0 , -C_rl*2*a );
						e -> add_term2( vars[2*numLabels + count_high*5] , vars[ind_r+numLabels] , 0 , 0 , 0 , -C_rl*2*a );

						e -> add_term2( vars[2*numLabels + count_high*5+1] , vars[ind_l]                 , 0 , 0 , 0 , -C_rl*2*a );  //w2 *( x0 + x1 + y1 )
						e -> add_term2( vars[2*numLabels + count_high*5+1] , vars[neigh_ind_l]           , 0 , 0 , 0 , -C_rl*2*a );
						e -> add_term2( vars[2*numLabels + count_high*5+1] , vars[neigh_ind_r+numLabels] , 0 , 0 , 0 , -C_rl*2*a );

						e -> add_term2( vars[2*numLabels + count_high*5+2] , vars[ind_l]                 , 0 , 0 , 0 , -C_rl*2*a );  //w3 *( x0 + y0 + y1 )
						e -> add_term2( vars[2*numLabels + count_high*5+2] , vars[ind_r+numLabels]       , 0 , 0 , 0 , -C_rl*2*a );
						e -> add_term2( vars[2*numLabels + count_high*5+2] , vars[neigh_ind_r+numLabels] , 0 , 0 , 0 , -C_rl*2*a );

						e -> add_term2( vars[2*numLabels + count_high*5+3] , vars[neigh_ind_l]           , 0 , 0 , 0 , -C_rl*2*a );  //w4 *( x1 + y0 + y1 )
						e -> add_term2( vars[2*numLabels + count_high*5+3] , vars[ind_r+numLabels]       , 0 , 0 , 0 , -C_rl*2*a );
						e -> add_term2( vars[2*numLabels + count_high*5+3] , vars[neigh_ind_r+numLabels] , 0 , 0 , 0 , -C_rl*2*a );

						e -> add_term2( vars[2*numLabels + count_high*5+4] , vars[ind_l]                 , 0 , 0 , 0 , C_rl*4*a );  //w*(x0+x1+y0+y1)
						e -> add_term2( vars[2*numLabels + count_high*5+4] , vars[neigh_ind_l]           , 0 , 0 , 0 , C_rl*4*a );
						e -> add_term2( vars[2*numLabels + count_high*5+4] , vars[ind_r+numLabels]       , 0 , 0 , 0 , C_rl*4*a );
						e -> add_term2( vars[2*numLabels + count_high*5+4] , vars[neigh_ind_r+numLabels] , 0 , 0 , 0 , C_rl*4*a );
					}
				}
				else
				{
					e -> add_term1(vars[2*numLabels + count_high*5]   , 0, C_rl*4*a);
					e -> add_term1(vars[2*numLabels + count_high*5 +1], 0, C_rl*4*a);
					e -> add_term1(vars[2*numLabels + count_high*5 +2], 0, C_rl*4*a);
					e -> add_term1(vars[2*numLabels + count_high*5 +3], 0, C_rl*4*a);
					e -> add_term1(vars[2*numLabels + count_high*5 +4], 0, C_rl*12*a);    // w1 w2 w3 w4 w

					for( int u=0 ; u<4 ; u++ )
					{
						int  neigh_y , neigh_xl , neigh_xr , ind_l , ind_r , neigh_ind_l , neigh_ind_r ;
						neigh_y   = h  + loc_y[u] ;
						neigh_xl  = c  + loc_x[u] ;
						neigh_xr  = cr + loc_x[u] ;

						ind_l = c  * numRows + h ;
						ind_r = cr * numRows + h ;
						neigh_ind_l = neigh_xl * numRows + neigh_y ;
						neigh_ind_r = neigh_xr * numRows + neigh_y ;
		
						e -> add_term1(vars[ind_l]                 , 0, C_rl*f1 );
						e -> add_term1(vars[ind_r+numLabels]       , 0, C_rl*f1 );
						e -> add_term1(vars[neigh_ind_l]           , 0, C_rl*f1 );
						e -> add_term1(vars[neigh_ind_r+numLabels] , 0, C_rl*f1 );         //x0 y0 x1 y1

						e -> add_term2( vars[ind_l]          , vars[neigh_ind_l]           , 0 , 0 , 0 , C_rl*5*a);   //x0x1 y0y1 x0y1 x1y0
						e -> add_term2( vars[ind_r+numLabels], vars[neigh_ind_r+numLabels] , 0 , 0 , 0 , C_rl*5*a);
						e -> add_term2( vars[ind_l]          , vars[neigh_ind_r+numLabels] , 0 , 0 , 0 , C_rl*5*a);
						e -> add_term2( vars[neigh_ind_l]    , vars[ind_r+numLabels]       , 0 , 0 , 0 , C_rl*5*a);

						e -> add_term2( vars[ind_l]      , vars[ind_r+numLabels]       , 0 , 0 , 0 , C_rl*(-2*f1+4*a) );
						e -> add_term2( vars[neigh_ind_l], vars[neigh_ind_r+numLabels] , 0 , 0 , 0 , C_rl*(-2*f1+4*a) );   //x0y0 x1y1

						e -> add_term2( vars[2*numLabels + count_high*5] , vars[ind_l]           , 0 , 0 , 0 , -C_rl*2*a );  //w1 *( x0 + x1 + y0 )
						e -> add_term2( vars[2*numLabels + count_high*5] , vars[neigh_ind_l]     , 0 , 0 , 0 , -C_rl*2*a );
						e -> add_term2( vars[2*numLabels + count_high*5] , vars[ind_r+numLabels] , 0 , 0 , 0 , -C_rl*2*a );

						e -> add_term2( vars[2*numLabels + count_high*5+1] , vars[ind_l]                 , 0 , 0 , 0 , -C_rl*2*a );  //w2 *( x0 + x1 + y1 )
						e -> add_term2( vars[2*numLabels + count_high*5+1] , vars[neigh_ind_l]           , 0 , 0 , 0 , -C_rl*2*a );
						e -> add_term2( vars[2*numLabels + count_high*5+1] , vars[neigh_ind_r+numLabels] , 0 , 0 , 0 , -C_rl*2*a );

						e -> add_term2( vars[2*numLabels + count_high*5+2] , vars[ind_l]                 , 0 , 0 , 0 , -C_rl*2*a );  //w3 *( x0 + y0 + y1 )
						e -> add_term2( vars[2*numLabels + count_high*5+2] , vars[ind_r+numLabels]       , 0 , 0 , 0 , -C_rl*2*a );
						e -> add_term2( vars[2*numLabels + count_high*5+2] , vars[neigh_ind_r+numLabels] , 0 , 0 , 0 , -C_rl*2*a );

						e -> add_term2( vars[2*numLabels + count_high*5+3] , vars[neigh_ind_l]           , 0 , 0 , 0 , -C_rl*2*a );  //w4 *( x1 + y0 + y1 )
						e -> add_term2( vars[2*numLabels + count_high*5+3] , vars[ind_r+numLabels]       , 0 , 0 , 0 , -C_rl*2*a );
						e -> add_term2( vars[2*numLabels + count_high*5+3] , vars[neigh_ind_r+numLabels] , 0 , 0 , 0 , -C_rl*2*a );

						e -> add_term2( vars[2*numLabels + count_high*5+4] , vars[ind_l]                 , 0 , 0 , 0 , -C_rl*8*a );  //w*(x0+x1+y0+y1)
						e -> add_term2( vars[2*numLabels + count_high*5+4] , vars[neigh_ind_l]           , 0 , 0 , 0 , -C_rl*8*a );
						e -> add_term2( vars[2*numLabels + count_high*5+4] , vars[ind_r+numLabels]       , 0 , 0 , 0 , -C_rl*8*a );
						e -> add_term2( vars[2*numLabels + count_high*5+4] , vars[neigh_ind_r+numLabels] , 0 , 0 , 0 , -C_rl*8*a );
					}

				}

		    	count_high ++ ;

		    }

    	}


	printf("Graph Made... \n");
}

void SegmentImage( Energy<float,float,float> *e , int *vars , int numLabels ,
				   int numRows, int numCols, double *SegImage_L, double *SegImage_R )
{
	int i,j;
	int Emin = e -> minimize();
    printf("Energy = %d\n", Emin);

	for( i=0 ; i<numRows ; i++ )
		for( j=0 ; j<numCols ; j++ )
		{
			int this_pix = j*numRows + i;

			if ( e->get_var(vars[this_pix]) == 1 )		// Do the classification...
				SegImage_L[j*numRows + i] = 1.0;
			else
				SegImage_L[j*numRows + i] = 0.0;						
		}

	for( i=0 ; i<numRows ; i++ )
		for( j=0 ; j<numCols ; j++ )
		{
			int this_pix = numLabels+ j*numRows + i;

			if ( e->get_var(vars[this_pix]) == 1 )		// Do the classification...
				SegImage_R[j*numRows + i] = 1.0;
			else
				SegImage_R[j*numRows + i] = 0.0;						
		}
}

///////////////////////////////////////////////////////////////////////////////

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	//the index of the variables from the left image are 0~numLabels-1 ; right:numLabels~2*numLabels-1
	double *r,*c, *Colors_L, *Colors_R, 
		   *FDist_Uc_L , *BDist_Uc_L , *FDist_Uc_R , *BDist_Uc_R , 
		   *FDist_Us_L , *BDist_Us_L , *FDist_Us_R , *BDist_Us_R , 
		   *disp , *RelateParam_s , *RelateParam_p , *RelateParam_c ;	
	// Input Arguments

	double *SegImage_L , *SegImage_R ;	// Output Arguments
	int numRows, numCols, numLabels ;
	int i;
	  
	if(nrhs!=16) 
	{
		mexErrMsgTxt("Incorrect No. of inputs");
	} 
	else if(nlhs!=2) 
	{
		mexErrMsgTxt("Incorrect No. of outputs");
	}
	  
   	r          = mxGetPr(prhs[0]);
	c          = mxGetPr(prhs[1]);
	Colors_L   = mxGetPr(prhs[2]);
	Colors_R   = mxGetPr(prhs[3]);
	FDist_Uc_L    = mxGetPr(prhs[4]);
	BDist_Uc_L    = mxGetPr(prhs[5]);
	FDist_Uc_R    = mxGetPr(prhs[6]);
	BDist_Uc_R    = mxGetPr(prhs[7]);
	FDist_Us_L    = mxGetPr(prhs[8]);
	BDist_Us_L    = mxGetPr(prhs[9]);
	FDist_Us_R    = mxGetPr(prhs[10]);
	BDist_Us_R    = mxGetPr(prhs[11]);
	disp  	      = mxGetPr(prhs[12]);	
	RelateParam_s  = mxGetPr(prhs[13]); 
	RelateParam_p  = mxGetPr(prhs[14]);  
 	RelateParam_c  = mxGetPr(prhs[15]);   

	double lambda_s = (*RelateParam_s) ;
    double lambda_p = (*RelateParam_p);
    double lambda_c = (*RelateParam_c);
    
	numCols = (int)(*c);		// Image Size
	numRows = (int)(*r);
	numLabels = numRows*numCols;	// Number of labels

	plhs[0] = mxCreateDoubleMatrix(numRows,numCols, mxREAL);	// Memory Allocated for output SegImage
	plhs[1] = mxCreateDoubleMatrix(numRows,numCols, mxREAL);	// Memory Allocated for output SegImage

	SegImage_L = mxGetPr(plhs[0]);	// variable assigned to the output SegImage
	SegImage_R = mxGetPr(plhs[1]);	// variable assigned to the output SegImage

	printf("Starting MEX Code \n");

	/**** Cast these data structures into correct types 
	and convert into STL Data Structures ****/

	printf("Vectorization 1: colors of pixels \n");
	vector< vector<double> > ColorsVec_L;
	ColorsVec_L.clear();
    for( i=0 ; i<numLabels ; i++ )
	{
		vector<double> tmp;
		tmp.clear();
		tmp.push_back( Colors_L[0*numLabels + i] );
		tmp.push_back( Colors_L[1*numLabels + i] );
		tmp.push_back( Colors_L[2*numLabels + i] );

		ColorsVec_L.push_back(tmp);
	}
	vector< vector<double> > ColorsVec_R;
	ColorsVec_R.clear();
    for( i=0 ; i<numLabels ; i++ )
	{
		vector<double> tmp;
		tmp.clear();
		tmp.push_back( Colors_R[0*numLabels + i] );
		tmp.push_back( Colors_R[1*numLabels + i] );
		tmp.push_back( Colors_R[2*numLabels + i] );

		ColorsVec_R.push_back(tmp);
	}


	printf("Vectorization 2 3 4 5: COLOR distance from foreground and background GMM \n");
	vector<double> FDistVec_Uc_L;
	FDistVec_Uc_L.clear();
    for(i=0;i<numLabels;i++)
		FDistVec_Uc_L.push_back(FDist_Uc_L[i]);

	vector<double> BDistVec_Uc_L;
	BDistVec_Uc_L.clear();
    for(i=0;i<numLabels;i++)
		BDistVec_Uc_L.push_back(BDist_Uc_L[i]);

	vector<double> FDistVec_Uc_R;
	FDistVec_Uc_R.clear();
    for(i=0;i<numLabels;i++)
		FDistVec_Uc_R.push_back(FDist_Uc_R[i]);

	vector<double> BDistVec_Uc_R;
	BDistVec_Uc_R.clear();
    for(i=0;i<numLabels;i++)
		BDistVec_Uc_R.push_back(BDist_Uc_R[i]);
	
	printf("Vectorization 6 7 8 9: SALIENCY distance from foreground and background GMM \n");
	vector<double> FDistVec_Us_L;
	FDistVec_Us_L.clear();
    for(i=0;i<numLabels;i++)
		FDistVec_Us_L.push_back(FDist_Us_L[i]);

	vector<double> BDistVec_Us_L;
	BDistVec_Us_L.clear();
    for(i=0;i<numLabels;i++)
		BDistVec_Us_L.push_back(BDist_Us_L[i]);

	vector<double> FDistVec_Us_R;
	FDistVec_Us_R.clear();
    for(i=0;i<numLabels;i++)
		FDistVec_Us_R.push_back(FDist_Us_R[i]);

	vector<double> BDistVec_Us_R;
	BDistVec_Us_R.clear();
    for(i=0;i<numLabels;i++)
		BDistVec_Us_R.push_back(BDist_Us_R[i]);


	printf("Vectorization 10: disparity map \n");
	vector<int> Disparity;
	Disparity.clear();
    for( int c=0 ; c<numCols ; c++ )
        for( int h=0 ; h<numRows ; h++ )
	    {
	    	Disparity.push_back( (int)disp[c*numRows+h] );
	    }


	/**************************************************************/
	/********** Data Structures Converted into Vector Form*********/
	/**************************************************************/



	/**** Actual Image Segmentation Code *********/
    // Adjacency list declaration
	vector< set<int> > neighbors;					

    // Adjacency list made
	MakeAdjacencyList( numLabels, numCols, numRows , neighbors );	

	printf("Adjacency lists made \n");
	
	//Calculate how many high-order cliques

	int num_high= 0 ;  //the number of high-order cliques
	num_high = get_num_corres( Disparity , numRows , numCols ) ;

	Energy<float,float,float> *e = new Energy<float,float,float>( 2*numLabels+4*num_high*5 , (2*numLabels+4*num_high*5)*3 ) ;
	int *vars = new int[2*numLabels+4*num_high*5];   

	MakeGraph( e , vars , neighbors, Disparity ,
		       ColorsVec_L , FDistVec_Uc_L, BDistVec_Uc_L , FDistVec_Us_L, BDistVec_Us_L,
	   		   ColorsVec_R , FDistVec_Uc_R, BDistVec_Uc_R , FDistVec_Us_R, BDistVec_Us_R,
			   numLabels,  num_high , numRows , numCols , lambda_s , lambda_p , lambda_c );
	printf("Graph constructed \n");
	
	SegmentImage( e , vars , numLabels , numRows, numCols,  SegImage_L , SegImage_R );
	printf("Image Segmented \n");
}
