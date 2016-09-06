main_SVM.m implements a support vector machine (SVM) for classification of any dataset.

Two types of SVMs can be trained one with a hard margin that can only classify linearly separable data whereas Soft margin SVM can also separate noisy data as well as quadratically separable data.

To run the code, you need Matlab installed or it can be easily translated to python and we can use numpy and pandas for matrix manipulation.

linear_new.mat,quadratic_new.mat and noisy-linear_new.mat are sample data sets used to test the algorithm. But any data set can be used such that each column is a sample of data whereas and number of rows determine dimension of feature vector.

compute_error.m, classify.m are used for estimating test error and plotSVM_Hplane.m plots the support vectors.


SoftMarg_gammaopt.m implements non-linear kernelled SVM for non-separable data.