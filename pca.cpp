#include "pca.h"
#include <vector>

using namespace std;
using namespace Eigen;

pca::pca()
{

}

int rows=20;
int col=10;
float parameter = 0.98;

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


//le tableau mean sert à stocker les moyennes
//on calcule les moyennes des variables aleatoires et on les centre
void adjust_data(MatrixXf& origin, VectorXf& means){
    for(int i=0; i<col; i++){
        float mean =0;
        for(int j=0; j<rows;j++){
            mean+= origin(j,i);
        }

        mean = (float) mean/rows;
        means(i) = mean;

        double norm=0;
        for(int j=0;j<rows;j++){
            origin(j,i)=origin(j,i) - mean;
            norm+= pow(origin(j,i),2);
        }

        for(int j=0; j<rows;j++){
            if(norm==0){origin(j,i)=0;}
            else{ origin(j,i)= (float) (origin(j,i)/sqrt(norm));}
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
EigenSolver<MatrixXf> eigen( MatrixXf& origin){
    EigenSolver<MatrixXf> es;
    es.compute(origin, true);
    return es;
}

void transpose(MatrixXf& origin, MatrixXf& transp){
    for(int i=0; i<origin.cols();i++){
        for(int j=0; j<origin.rows(); j++){
            transp(i,j)=origin(j,i);
        }
    }
}

//retourne A*B dans C
void multiply(MatrixXf& A ,MatrixXf& B, MatrixXf& C){
    int a = A.rows();
    int n = A.cols();
    int b = B.cols();
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

int decreaseDim(VectorXf& eigenVal){
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
    return indice;
}


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
    }
}

template<typename M>
M load_csv (const std::string & path) {
    std::ifstream indata;
    indata.open(path);
    std::string line;
    std::vector<float> values;
    int r = 0;

    while (std::getline(indata, line)) {
        int colonnes=0;
        std::stringstream lineStream(line);
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            values.push_back(std::stod(cell));
            colonnes+=1;

        }
        col=colonnes-1;
        ++r;
    }

    rows=r;
    return Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, RowMajor>>(values.data(), r, values.size()/r);
}

void principalComponentAnalysis(const std::string& pathData){
    //MatrixXf origin = load_csv<MatrixXf>(pathData);

     MatrixXf origin2 = load_csv<MatrixXf>(pathData);
     MatrixXf origin(rows,col);

     cout<< "lignes origin " << origin2.rows() << "  "<< rows<<endl;
     cout<<"colonnes origin "<<origin2.cols()<<"  "<< col <<endl;
     for(int i=0; i<rows;i++){
     for(int j=0;j<col;j++){
     origin(i,j)=origin2(i,j+1);
     //origin(i,j)=rand();
     }
     }
    VectorXf means(col);
    adjust_data(origin,means);
    MatrixXf covar_matrix(col,col);
    compute_covariance_matrix(origin,covar_matrix);

    cout<< "covariance calculee" <<endl;
    cout<< covar_matrix<<endl;
    cout<<" covar rows "<< covar_matrix.rows()<<endl;
    cout<<"covar cols "<<covar_matrix.cols()<<endl;


    EigenSolver<MatrixXf> eig = eigen(covar_matrix);
    VectorXf eigenVal(col);
    MatrixXf final_data(rows,col);
    MatrixXf eigenVec(col,col);

    cout<<"The eigenvalues"<<endl;
    cout<<eig.eigenvalues().transpose()<<endl;
    cout<<"The eigenvectors"<<endl;

    for(int i=0; i<col;i++){
        eigenVal(i)= eig.eigenvalues()[i].real();
        for(int j=0; j<col;j++){
            eigenVec(i,j)= eig.eigenvectors().col(j)[i].real();
        }
    }


    cout<<eigenVec<<endl;

    int nbVec= decreaseDim(eigenVal);
    cout<<"nbVec "<<nbVec<<endl;
    cout<<"vecteurs propres réduits"<<endl;

    MatrixXf baseChange(rows,col);
    MatrixXf transposeEigen(col,col);
    transpose(eigenVec,transposeEigen);
    cout<<"transpose Eigen "<< transposeEigen<<endl;
    cout<< "origin rows cols "<< origin.rows() <<" "<<origin.cols()<<endl;
    cout<< "transEig rows cols "<< transposeEigen.rows() <<" "<<transposeEigen.cols()<<endl;
    cout<< "basChaneg rows cols "<< baseChange.rows() <<" "<<baseChange.cols()<<endl;
    multiply(origin,transposeEigen,baseChange);
    cout<<"baseChange "<<endl;
    cout<<baseChange<<endl;

    MatrixXf final(rows,nbVec);
    for(int i=0; i<rows;i++){
        for(int j=0; j<nbVec;j++){
            final(i,j)=baseChange(i,j);
        }
    }
    cout<<"donnees finales "<<endl;
    cout<<final<<endl;

    cout<<"lignes finales "<< final.rows()<<endl;
    cout<<"colonnes finales "<< final.cols()<<endl;

    writeToText("resultDimReduction", final);
}
