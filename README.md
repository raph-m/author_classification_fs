# author_classification_fs
The principal function is ```principalComponentAnalysis(const std::string& pathData, const std::string& pathTest, float parameter)``` in pca.pp.
This function computes the base change required for the data from pathData, the change base enables to reduce the number of features analyzed in the machine learning program.
The same base change is applied to the data from the texts to test, so that the machine learning algorithm runs on homogeneous data.
The float parameter influences the number of dimensions kept : ```sum of the eigenvalues selected = parameter * the total sum of the eigenvalues```.
