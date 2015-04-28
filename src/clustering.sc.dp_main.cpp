/*
clustering_dp_main.cpp --- wrap function for 
						   "clustering_sc_dp()"
Based on Ckmeans.1d.dp_main.cpp

  Tibor Szkaliczki
  eLearning Department
  Institute for Automation and Control
  szkaliczki.tibor@sztaki.mta.hu
				  
Created: March 14, 2015
*/

#include "clustering.sc.dp.h"
#include <vector>   

using namespace std;

/* Wrap functions */
extern "C" {
	
	void Cclustering_sc_dp(double *x, int *s, int *fs, int *k, int* cluster, double* centers, double* withinss, int* size) {
		int length = s[0];
		int feature_size = fs[0];
		int level = k[0];
	
		/* input data array */
		vector<vector<double> > input(length +1, vector<double>(feature_size));

		for(int j = 0; j < feature_size; j++) {
			for(int i=1;i<(length+1);i++) {
				input[i][j] = x[j * length + i - 1];
			}
		}

		/* vector of total distances */
		vector<double> d(level + 1);

		/* array for backtracking */
		vector<vector<int> > B (level + 1, vector<int>(length+1));
				
		/* Call C++ version of dp clustering algorithm */
		clustering_sc_dp(input, level, d, B);

		data result;  /*object of class data, used to convey the final result*/

		/*Call C++ version of backtracking algorithm*/
		result = backtracking( input, level, B);

		/*Since R doesn't allow return value from C/C++ function, using pointers to give back result*/
		for(int i=1;i<(length+1);i++)
			cluster[i-1] = result.cluster[i];
		for (int i = 1; i < (level + 1); i++) {
			for (int j = 0; j < feature_size; j++) {
				centers[(i - 1) * feature_size + j] = result.centers[i][j];
			}
			withinss[i - 1] = result.withinss[i];
			size[i - 1] = result.size[i];
		}
	}

	void Cpreclustering_sc_dp(double *x, int *s, int *fs, int *k, double *d, int *backtrack) {
		int length = s[0];
		int feature_size = fs[0];
		int level = k[0];

		/* input data array */
		vector<vector<double> > input(length +1, vector<double>(feature_size));

		for(int j = 0; j < feature_size; j++) {
			for(int i=1;i<(length+1);i++) {
				input[i][j] = x[j * length + i - 1];
			}
		}
		vector<double> dist(level + 1);
		vector<vector<int> > B (level + 1, vector<int>(length+1));
		
		/* Call C++ version of dp clustering algorithm */
		clustering_sc_dp(input, level, dist, B);

		for(int i=0;i<level;i++) {
			d[i] = dist[i];
			for(int j = 0; j < length; j++) {
				backtrack[i * length + j] = B[i+1][j+1];
			}
		}
	}

	void Cbacktracking(double *x, int *s, int *fs, int *k, int *backtrack, int *hbacktrack, int* cluster, double* centers, double* withinss, int* size) {
		
		int length = s[0];
		int feature_size = fs[0];
		int level = k[0];
		int maxlevel = hbacktrack[0];
	
		/* input data array */
		vector<vector<double> > input(length +1, vector<double> (feature_size));

		for(int j = 0; j < feature_size; j++) {
			for(int i=1;i<(length+1);i++) {
				input[i][j] = x[j * length + i - 1];
			}
		}

		vector<vector<int> > B (maxlevel + 1, vector<int>(length+1));

		for(int j=0;j<length;j++) {
			for(int i = 0; i < maxlevel; i++) {
				B[i+1][j+1] = backtrack[j * maxlevel + i];
			}
		}
		
		data result;  /*object of class data, used to convey the final result*/

		/*Call C++ version of backtracking algorithm*/
		result = backtracking( input, level, B);
		
		/*Since R doesn't allow return value from C/C++ function, using pointers to give back result*/
		for(int i=1;i<(length+1);i++)
			cluster[i-1] = result.cluster[i];
		for (int i = 1; i < (level + 1); i++) {
			for (int j = 0; j < feature_size; j++) {
				centers[(i - 1) * feature_size + j] = result.centers[i][j];
			}
			withinss[i - 1] = result.withinss[i];
			size[i - 1] = result.size[i];
		}
	}
}	
