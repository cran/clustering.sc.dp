/*
clustering_sc_dp.h --- Head file for clustering_sc_dp
                    Declare a class "data" and 
					wrap functions related to "clustering_sc_dp"
Based on Ckmeans.1d.dp.h
              
  Tibor Szkaliczki
  eLearning Department
  Institute for Automation and Control
  szkaliczki.tibor@sztaki.mta.hu
				  
Created: April 10, 2015
*/

#include <vector>
using namespace std;

/* data that return by clustering.sc.dp()*/
class data{
public:
    vector<int> cluster;  	/*record which cluster each point belongs to*/
    vector<vector<double> > centers;	/*record the center of each cluster*/
    vector<double> withinss;/*within sum of distance square of each cluster*/
    vector<int> size;		/*size of each cluster*/
};

/*one-dimensional cluster algorithm implemented in C*/
/*x is input one-dimensional vector and K stands for the cluster level*/
void clustering_sc_dp(vector<vector<double> > &x, int K, vector<double> &dist, vector<vector<int> > &B);
//vector<double> getDistances(vector<vector<double>> &D);
data backtracking(vector<vector<double> > &x, int K, vector<vector<int> > &B);
vector<data> getResults(vector<vector<double> > &x, int k1, int k2, vector<vector<int> > &B);
