#ifndef PCA_H
#define PCA_H

#include <iostream>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include <map>
#include <assert.h>

#include <stdio.h>
#include <stdlib.h>
#include <Eigen/Eigenvalues>

#include <cmath>

class pca
{
public:
    pca();
};

double mean(Eigen::MatrixXf& origin, int i);
double variance(Eigen::MatrixXf&  origin, double& mean,int i);
void adjust_data(Eigen::MatrixXf& origin, Eigen::VectorXf& means);
float compute_covariance(const Eigen::MatrixXf& d, int i, int j);
void compute_covariance_matrix( const Eigen::MatrixXf& d, Eigen::MatrixXf& covar_mat);
Eigen::EigenSolver<Eigen::MatrixXf> eigen(Eigen::MatrixXf& origin);
void transpose(Eigen::MatrixXf origin, Eigen::MatrixXf& transp);
Eigen::MatrixXf multiply(const Eigen::MatrixXf A ,const Eigen::MatrixXf B);
int decreaseDim(Eigen::VectorXf& eigenVal, float parameter);
void writeToText(std::string name, Eigen::MatrixXf matrix);
template<typename M> M load_csv (const std::string & path);
void principalComponentAnalysis(const std::string& pathData, const std::string& pathTest, float parameter);


#endif // PCA_H
