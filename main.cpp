#include <iostream>
#include "pca.cpp"
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

#include "pca.h"


int main(int argc, char *argv[]){

        MatrixXf origin(rows,col);
        for(int i=0; i<rows;i++){
            for(int j=0;j<col;j++){
                origin(i,j)=2.0+j+i;
            }
        }
        cout<< mean(origin, 1)<<endl;
        cout<<origin<<endl;
        VectorXf means(col);
        adjust_data(origin,means);
        cout<< "data ajustee" <<endl;

        MatrixXf covar_matrix(col,col);

        compute_covariance_matrix(origin,covar_matrix);
        cout<< "covariance calculee" <<endl;
        cout<< covar_matrix<<endl;

        //int dim = col;


        EigenSolver<MatrixXf> eig = eigen(covar_matrix);
        cout<<"The eigenvalues"<<endl;
        cout<<eig.eigenvalues().transpose()<<endl;
        cout<<"The eigenvectors"<<endl;

       VectorXf eigenVal(col);
       MatrixXf final_data(rows,col);
       MatrixXf eigenVec(col,col);
       for(int i=0; i<col;i++){
           eigenVal(i)= eig.eigenvalues()[i].real();
           for(int j=0; j<col;j++){
               eigenVec(i,j)= eig.eigenvectors().col(j)[i].real();
           }
       }
       cout<<eigenVec<<endl;

       MatrixXf finalEigenVec = decreaseDim(eigenVec, eigenVal, 0.95);
       int nbVec= finalEigenVec.cols();
       cout<<"vecteurs propres réduits"<<endl;
       cout<<finalEigenVec<<endl;

       MatrixXf baseChange(rows,col);
       MatrixXf transposeEigen(col,col);
       transpose(eigenVec,transposeEigen);
       multiply(origin,transposeEigen,baseChange);

       MatrixXf final(rows,nbVec);
       for(int i=0; i<rows;i++){
           for(int j=0; j<nbVec;j++){
               final(i,j)=baseChange(i,j);
           }
       }
       cout<<"donnees finales "<<endl;
       cout<<final<<endl;

       writeToText("test", final);


return 0;


}
