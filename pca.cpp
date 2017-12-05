#include "pca.h"

using namespace std;
using namespace Eigen;

pca::pca()
{

}

const int rows=10;
const int col = 10;

double mean(MatrixXf& origin, int i){
    double mean =0;
    for(int j=0; j<rows; j++){
        mean+= origin(j,i);
    }
    mean= (double)mean/rows;
    return mean;
}


double variance(MatrixXf&  origin, double& mean,int i){
    double variance=0;
    for(int j=0; j<rows; j++){
        variance+=pow(origin(i,j),2);
    }
    variance= (double)variance/rows;
    variance = variance - pow(mean,2);
    return variance;
}

//Dans le cas ou on centre et reduit la matrice d'origine
/*double* transformation(double* origin){
    newMat= new double* [rows];
    for(int i=0; i<col; i++){
        double m = mean(origin,i);
        double v = variance(origin,m,i);
        for(int j=0; j<rows; j++){
            newMat[j] = new double*[col];
            newMat[j][i]= double(origin[i][j]-m)/v;
        }
    }
    return newMat;
}*/

//le tableau mean sert à stocker les moyennes
//on calcule les moyennes des variables aleatoires et on les centre
void adjust_data(MatrixXf&origin, VectorXf& means){
    for(int i=0; i<col; i++){
        double mean =0;
        for(int j=0; j<rows;j++){
            mean+= origin(i,j);
        }
        mean = (double) mean/rows;
        means(i) = mean;
        for(int j=0;j<rows;j++){
            origin(j,i)-=mean;
        }
    }
}



double compute_covariance(const MatrixXf& d, int i, int j){
    double cov=0;
    int n = d.rows();
    for(int k=0; k<n;k++){
        cov+= d(k,i)*d(k,j);
    }
    return( (double)cov/(d.rows()));
}

void compute_covariance_matrix( const MatrixXf& d, MatrixXf& covar_mat){
    int dim= d.cols();
    if(dim== covar_mat.rows() && dim==covar_mat.cols()){
        for(int i=0; i<dim;i++){
            for(int j=i; j<dim; j++){
                covar_mat(i,j)= compute_covariance(d,i,j);
            }
        }
        //on complete la deuxieme moitiee de la matrice
        for(int i=1; i<dim;i++){
            for(int j=0; j<i; j++){
                covar_mat(i,j)=covar_mat(j,i);
            }
        }
    }
}

//Calcul des valeurs et vecteurs propres de la matrice de covariance
/*void eigen(MatrixXf& covar_mat, MatrixXf& eigenvector, VectorXf& eigenvalues){
    //Eigenvalue<double> eig(covar_matrix);

    eig.getV(eigenvector);
    eig.getD(eigenvalues);
}*/


EigenSolver<MatrixXf> eigen( MatrixXf& origin){
    EigenSolver<MatrixXf> es;
    es.compute(origin, true);
    return es;
}

void transpose(MatrixXf& origin, MatrixXf& transp){
    for(int i=0; i<col;i++){
        for(int j=0; j<rows; j++){
            transp(i,j)=origin(j,i);
        }
    }
}


//retourne A*B dans C
void multiply(MatrixXf& A ,MatrixXf& B, MatrixXf& C){
    int a = A.rows();
    int b = B.cols();
    int n = A.cols();
    if(n!=(int)B.rows()){cout<<"erreur de dimension"<<endl;}
    else{
        for(int i=0; i<a; i++){
            for(int j=0; j<b; j++){
                for(int k=0; k<n;k++){
                    C(i,j)+= A(i,k)*B(k,j);
                }
            }
        }
    }
}

MatrixXf decreaseDim(MatrixXf& eigenVec, VectorXf& eigenVal, float parameter){
    float s=0;
    for(int i=0; i<eigenVal.size(); i++){
        s+=eigenVal(i);
    }
    s*=parameter;
    float sum=0;
    int indice=0;

    while(sum<s){
        sum+=eigenVal(indice);
        indice++;
    }
    int r= eigenVec.rows();
    MatrixXf final(r,indice) ;

    for(int i=0;i<r;i++){
        for(int j=0;j<indice;j++){
            final(i,j)=eigenVec(i,j);
        }
    }
    return final;


}

/*void writeToCSVfile(string name, MatrixXf matrix){
    ofstream file(name.c_str());
    file<<matrix.format(CSVFormat);
}*/

void writeToText(string name, MatrixXf matrix){
    std::ofstream file(name+".csv");
    if(file.is_open()){
        for(int i=0; i<matrix.rows();i++){
            for(int j=0; j<matrix.cols();j++){
                if( j<matrix.cols()-1){file<<matrix(i,j)<<",";}
                else{ file<<matrix(i,j);}
            }
            file<<endl;
        }
       // file<<matrix;
    }
}


/*//calcul de la matrice des correlations a partir d'une matrice centrée réduite
double* correlation( double* reduced){
    int k = reduced.length;
    double*  treduced= transpose(reduced);
    double* C = multiply(treduced,reduced);
    for(int i=0; i<C.length;i++){
        for(int j=0; j<C[0].length; j++){
            C[i][j]=C[i][j] * (double)(1/k);
        }
    }
    return C;
}*/

