#include "pca.h"
#include <vector>

using namespace std;
using namespace Eigen;

pca::pca()
{

}

float parameter = 0.98;

double mean(MatrixXf& origin, int i){
    int rows = origin.rows();
    double mean =0;
    for(int j=0; j<rows; j++){
        mean+= origin(j,i);
    }
    mean= (double)mean/rows;
    return mean;
}


double variance(MatrixXf&  origin, double& mean,int i){
    int rows = origin.rows();
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
    int col = origin.cols();
    int rows = origin.rows();
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

float compute_covariance(const MatrixXf& d, int i, int j){
    float cov=0;
    int n = d.rows();
    for(int k=0; k<n;k++){
        cov+= d(k,i)*d(k,j);
    }
    return( (float)cov/(d.rows()-1));
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
    else{
        cout<<"pb de dimension calcul de covariance "<<endl;
    }
}

//Calcul des valeurs et vecteurs propres de la matrice de covariance
EigenSolver<MatrixXf> eigen( MatrixXf& origin){
    EigenSolver<MatrixXf> es;
    es.compute(origin, true);
    return es;
}

void transpose(MatrixXf origin, MatrixXf& transp){
    for(int i=0; i<origin.cols();i++){
        for(int j=0; j<origin.rows(); j++){
            transp(i,j)=origin(j,i);
        }
    }
}

//retourne A*B dans C
MatrixXf multiply(const MatrixXf A ,const MatrixXf B){
    if(A.cols()!=(int)B.rows()){cout<<"erreur de dimension"<<endl;}
    else{
        MatrixXf C(A.rows(),B.cols());
        for(int i=0; i<A.rows(); i++){
            for(int j=0; j<B.cols(); j++){
                for(int k=0; k<A.cols();k++){
                    C(i,j)+= A(i,k)*B(k,j);
                }
            }
        }
          return C;
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
   /* std::ofstream file(name+".csv");
    if(file.is_open()){
        for(int i=0; i<matrix.rows();i++){
            for(int j=0; j<matrix.cols();j++){
                if( j<matrix.cols()-1){file<<matrix(i,j)<<",";}
                else{ file<<matrix(i,j);}
            }
            file<<endl;
        }
    }*/
    std::string a = "../resultaPCA.csv";

            std::ofstream myfile(a.c_str(), std::ofstream::out | std::ofstream::app);
            if (myfile.is_open()){
                // myfile.seekp(0, ios::end); // On se déplace à la fin du fichier
                for(int i=0;i<matrix.rows();i++){
                    for(int j=0; j<matrix.cols();j++){
                        if( j<matrix.cols()-1){myfile<<matrix(i,j)<<",";}
                        else{ myfile<<matrix(i,j);}
                    }
                    myfile <<endl;
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
        //col=colonnes-1;
        ++r;
    }

    //rows=r;
    return Map<const Matrix<typename M::Scalar, M::RowsAtCompileTime, M::ColsAtCompileTime, RowMajor>>(values.data(), r, values.size()/r);
}

void principalComponentAnalysis(const std::string& pathData, const std::string& pathTest){

     MatrixXf origin2 = load_csv<MatrixXf>(pathData);
     cout <<origin2<<endl;
     int rows = origin2.rows();
     int col = origin2.cols()-1;
    //int rows=20;
    //int col=5;
    //int rowsT=20;
     MatrixXf origin(rows,col);
     MatrixXf test2 = load_csv<MatrixXf>(pathTest);
     int rowsT = test2.rows();
     MatrixXf test (rowsT,col);
     for(int i=0; i<rows;i++){
        for(int j=0;j<col;j++){
            origin(i,j)=origin2(i,j+1);
          // origin(i,j)=i*(j+2)/3;//(float)rand()/RAND_MAX;
        }
     }
     for(int i=0; i<rowsT;i++){
        for(int j=0;j<col;j++){
            test(i,j)=test2(i,j+1);
            // test(i,j)=rand();
        }
     }

    VectorXf means(col);
    VectorXf meansTest(col);
    adjust_data(origin,means);
    adjust_data(test,meansTest);
    cout<<"origin "<<origin<<endl;
    cout<<"test "<<test<<endl;
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

    MatrixXf transposeEigen(col,col);

    transpose(eigenVec,transposeEigen);
    cout<<"transpose Eigen"<<endl;
    cout<<transposeEigen<<endl;

    MatrixXf baseChange = multiply(origin,transposeEigen);
    cout<<"base chaneg"<<endl;
    cout<<baseChange<<endl;
    MatrixXf baseChangeTest = multiply(test,transposeEigen);

    cout<<"base changeT"<<endl;
    cout<<baseChangeTest<<endl;


    cout<<"origin "<<origin<<endl;
    cout<<"test "<<test<<endl;

    MatrixXf final(rows,nbVec+1);
    MatrixXf finalTest(rowsT,nbVec+1);

    for(int i=0; i<rows;i++){
        for(int j=0; j<nbVec;j++){
            if(j==0){
               final(i,0)=origin2(i,0);
            }
                final(i,j+1)=baseChange(i,j);

        }
    }
    for(int i=0; i<rowsT;i++){
        for(int j=0; j<nbVec;j++){
            if(j==0){
                finalTest(i,j)=test2(i,0);
            }
                finalTest(i,j+1)=baseChangeTest(i,j);
        }
    }
    cout<<"donnees finales "<<endl;
    cout<<final<<endl;
    cout<<"donnees finales test "<<endl;
    cout<<finalTest<<endl;




    //writeToText("resultDimReduction", final);
}
