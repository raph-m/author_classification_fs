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
        std::string c = "/Users/mariannedevriendt/Desktop/results.txt";
        principalComponentAnalysis(c);
        return 0;
}
