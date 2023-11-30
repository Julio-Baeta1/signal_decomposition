#ifndef ICA_H
#define ICA_H

#include <memory>
#include <random>
#include <iostream>
#include <fstream>
#include<Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <functional>
#include <utility>

using Mat = Eigen::MatrixXd;

class Ica{
    private:
        /* XX: Provided mixed signals
            X: ICA unmixed signals  
            W: Unmixing matrix
        */
        std::unique_ptr<Mat> X,W,XX;

        void gsDecorr(Mat *w, int col_num);
        void symDecorr();

        void serialFastICA(int n_sigs, double tol, int max_iter, std::function<double(double)> g, std::function<double(double)> g_der);
        void parallelFastICA(int n_sigs, double tol, int max_iter, std::function<double(double)> g, std::function<double(double)> g_der);

    public:

        Ica(Mat ini_X = Mat::Zero(1,1));
        Ica(Mat* ini_ptr);

        void setSource(Mat new_X);
        void setSource(Mat* new_ptr);
        void setSourceFromFile(std::string filename);

        Mat* getResultPtr() {return X.get();}
        Mat* getMixedPtr() {return XX.get();}
        Mat* getMixingPtr() {return W.get();}
        Mat getResult() {return *X;}
        Mat getMixed() {return *XX;}
        Mat getMixingMat() {return *W;}


        void randW(int rows, int cols, int seed);
        void setW(int rows, int cols, int seed, bool is_rand);

        void sphering();
        void decompose(int n_sigs=2, bool rand_W=false, int seed=1);
        void fastIca(int n_sigs=2, std::string func_type="cosh", bool paral = false, int seed=1, double tol=1e-6, int max_iter=200);
};

#endif