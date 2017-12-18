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
        float parameter = 0.98;
        std::string c = "../apprentissage_10.txt";
        std::string test = "../testResult.txt";
        principalComponentAnalysis(c,test, parameter);
        return 0;
}
