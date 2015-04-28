/*
clustering_sc_dp.cpp -- Performs clustering multidimensional data by a dynamic programming
                     approach that is guaranteed to be optimal if only subsequent items may form a cluster
Based on Ckmeans_1d_dp.cpp -- The original algorithm works only on one dimensional data and returns clustering for fixed number of clusters k. It is adapted to multidemensional data and returns the optimum for any number of clusters less than or equal to k.

  Tibor Szkaliczki
  eLearning Department
  Institute for Automation and Control
  szkaliczki.tibor@sztaki.mta.hu

Created: September 2, 2013
*/

#include "clustering.sc.dp.h"
#include <iostream>
#include <algorithm>
#include <string>
#include <math.h>

#include <assert.h>

using namespace std;

//all vectors in this program is considered starting at position 1, position 0 is not used.
void clustering_sc_dp(vector<vector<double> > &x, int K, vector<double> &dist, vector<vector<int> > &B)
{
  // Input:
  //  x -- arrays of numbers containing the feature vectors of subsequent items
  //  K -- the maximum number of clusters expected

	int N = x.size()-1;  //N: is the size of input vector
	int feature_size = x[1].size();
	vector<vector<double> > D(K + 1, vector<double>(N + 1));

	for(int i=1;i<=K;i++)
    {
        D[i][1] = 0;
        B[i][1] = 1;
    }
    
    vector <double> mean_x1(feature_size);
	vector <double> mean_xj(feature_size);
	vector <double> d_j(feature_size);
	vector <double> d_1(feature_size);
	double d;

	for(int i = 0; i < feature_size; i++ ) {
		mean_x1[i] = x[1][i];
		d_1[i] = 0.0;
	}

   for(int k=1;k<=K;k++)
    {
      for(int i=max(2,k);i<=N;i++)		// the number of clusters is less than or equal to the number of items
      {
        D[k][i] = -1;
        if(k == 1) 
        {
			d = 0.0;
			for(int l = 0; l < feature_size; l++ ) {
				d_1[l] += (i-1)/(double)i * (x[i][l] - mean_x1[l]) * (x[i][l] - mean_x1[l]);
				mean_x1[l] = ((i-1) * mean_x1[l] + x[i][l])/(double)i;
				d += d_1[l];  
			}
			D[1][i] = d;
			B[1][i] = 1;
        }
        else 
        {
			for(int l = 0; l < feature_size; l++ ) {
				d_j[l] = 0;
				mean_xj[l] = 0;
			}
          
          for(int j=i;j>=k;j--)				// the first item of the kth cluster is at least the kth item
          {
            d = 0.0;
			for(int l = 0; l < feature_size; l++ ) {
              d_j[l] += (i-j)/(double)(i-j+1) * (x[j][l] - mean_xj[l]) * (x[j][l] - mean_xj[l]);
              mean_xj[l] = (x[j][l] + (i-j)*mean_xj[l])/(double)(i-j+1);
			  d += d_j[l];
			}
            if(D[k][i] == -1) //initialization of D[k,i]
            { 
            
              if(j == 1) 
              {        
                D[k][i] = d;
                B[k][i] = j;
              } 
              else 
              { 
				if(D[k-1][j-1] == -1) {
					continue;
				}
				D[k][i] = d + D[k-1][j-1];
				B[k][i] = j;
              }
            } 
            else 
            {
              if(j == 1) 
              {
                if(d <= D[k][i]) 
                {
                  D[k][i] = d;
                  B[k][i] = j;
                }
              } 
              else 
              {
				if(D[k-1][j-1] == -1) {
					continue;
				}
				if(d + D[k-1][j-1] < D[k][i]) 
				{
				  D[k][i] = d + D[k-1][j-1];
				  B[k][i] = j;
				}
              }
            }
          }
        }
      }
   }
   dist.resize(K+1);
   for(int i = 0; i<K; i++) {
	dist[i] = D[i+1][N];
   }
}

data backtracking(vector<vector<double> > &x, int K, vector<vector<int> > &B) {

	data result;
	int N = x.size()-1;  //N: is the size of input vector
	int feature_size = x[1].size();

   //Backtrack to find the clusters of the data points
    int cluster_right = N;
    int cluster_left;
    result.cluster.resize(N+1);
    result.centers.resize(K+1);
    for( vector<vector<double> >::iterator it = result.centers.begin(); it != result.centers.end(); ++it) {
        it->resize(feature_size);
    }
    result.withinss.resize(K+1);
	result.size.resize(K+1);

	/*Forming final result*/
    for(int k=K;k>=1;k--)
    {
      cluster_left = B[k][cluster_right];
      
      for(int i=(int)cluster_left;i<=cluster_right;i++)
      	result.cluster[i] = k;

      vector<double> sum(feature_size);
	  for(int l = 0; l < feature_size; l++ ) {
		sum[l] = 0.0;
	  }
		  
      for(int a=(int)cluster_left;a<=cluster_right;a++)
		for(int l = 0; l < feature_size; l++ ) {
			sum[l] += x[a][l];
		}
		
		for(int l = 0; l < feature_size; l++ ) {
			result.centers[k][l] = sum[l]/(cluster_right-cluster_left+1);
		}

	  double error;
	  result.withinss[k] = 0.0;
      for(int a=(int)cluster_left;a<=cluster_right;a++) {
		  error = 0.0;
		for(int l = 0; l < feature_size; l++ ) {
	      	 error += pow( (x[a][l] - result.centers[k][l]), 2 );
		}
		result.withinss[k] += error;
	  }

       // Instead of avarage, the center contains the item closest to the average
		
      result.size[k] = cluster_right - cluster_left + 1;
      
      if(k > 1) {
        cluster_right = cluster_left - 1;
      }
    }

  return result;
}  


